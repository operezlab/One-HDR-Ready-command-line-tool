#!/usr/bin/env python3

import requests, sys, subprocess, csv,time, os, tempfile, json
import pandas as pd
from Bio.Seq import Seq
from Bio.Data import CodonTable     
import primer3
from pathlib import Path

# Function to get the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Function to initiate ID mapping
def initiate_id_mapping(ids, from_db, to_db):
    url = "https://rest.uniprot.org/idmapping/run"
    response = requests.post(url, data={'ids': ids, 'from': from_db, 'to': to_db})
    response.raise_for_status()
    return response.json()["jobId"]
    time.sleep(1)

# Function to check ID mapping results
def check_id_mapping_results(job_id):
    url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

# Function to convert a single Ensembl ID using a GET request
def convert_single_id(ensembl_id):
    url = f"https://biotools.fr/human/ensembl_symbol_converter/?api=1&id={ensembl_id}"
    response = requests.get(url)
    return response.json()

# Function to convert multiple Ensembl IDs using a POST request
def convert_multiple_ids(ensembl_ids):
    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    ids_json = json.dumps(ensembl_ids)
    body = {'api': 1, 'ids': ids_json}
    response = requests.post(url, data=body)
    return response.json()

# Get the Gene IDs from user input
def main(gene_ids):
    print(f"Received Gene IDs: {gene_ids}", flush=True)
    print("Processing CRISPR data...", flush=True)

    job_id = initiate_id_mapping(gene_ids, "GeneCards", "UniProtKB")
    print(f"Job ID: {job_id}")
    time.sleep(1)

    while True:
        results = check_id_mapping_results(job_id)
        if "results" in results:
            break
        print("Waiting for results...")
        time.sleep(1)

    uniprot_accession_codes = [result["to"] for result in results["results"]]
    print("UniProt Accession Codes:", uniprot_accession_codes)

    if uniprot_accession_codes:
        requestURL = f"https://www.ebi.ac.uk/proteins/api/coordinates/{uniprot_accession_codes[0]}"
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        response_data = r.json()
        relevant_info = []

        # Debugging: Print entire API response
        #print("API Response Data:", response_data)
        valid_human_chromosomes = {str(i) for i in range(1, 23)} | {"X", "Y"}

        for gene in response_data.get("gnCoordinate", []):
            ensembl_gene_id = gene.get("ensemblGeneId")
            genomic_location = gene.get("genomicLocation", {})

            if "exon" in genomic_location:
                for exon in genomic_location["exon"]:
                    exon_id = exon.get("id")
                    chromosome = str(genomic_location.get("chromosome"))  # Convert to string for comparison
                    start = exon.get("genomeLocation", {}).get("begin", {}).get("position")
                    end = exon.get("genomeLocation", {}).get("end", {}).get("position")

                    # Ensure chromosome is valid
                    if chromosome in valid_human_chromosomes:
                        #print(f"Valid chromosome found: {chromosome}")
                        if start is not None and end is not None:
                            relevant_info.append({
                                "ensembl_gene_id": ensembl_gene_id,
                                "exon_id": exon_id,
                                "chromosome": chromosome,
                                "start": start,
                                "end": end
                            })
                    else:
                        print(f"Skipping invalid chromosome: {chromosome}")

        if relevant_info:
            # Determine the last exon based on strand orientation
            reverse_strand = genomic_location.get("reverseStrand", False)
            if reverse_strand:
                # For reverse strand, the exon with the smallest start is the last
                last_exon_info = min(relevant_info, key=lambda x: x['start'])
            else:
                # For positive strand, the exon with the largest end is the last
                last_exon_info = max(relevant_info, key=lambda x: x['end'])

            print("Last exon information:")
            print(last_exon_info)

            # Define chromosome, exon_start, and exon_end explicitly
            chromosome = last_exon_info.get("chromosome")
            exon_start = last_exon_info.get("start")
            exon_end = last_exon_info.get("end")

            print(f"Chromosome: {chromosome}, Start: {exon_start}, End: {exon_end}")

            return last_exon_info
        else:
            print("No relevant exon information found for valid gene IDs.")
            return None
    else:
        print("No UniProt accession codes found.")
        return None

# ---------- small utils ----------
def rc(seq: str) -> str:
    t = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(t)[::-1]

def _run(cmd):
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if cp.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}\n{cp.stderr}")
    return cp.stdout

def nw_global(q: str, s: str, match=1, mismatch=-1, gap=-1):
    n, m = len(q), len(s)
    dp = [[0]*(m+1) for _ in range(n+1)]
    tb = [[None]*(m+1) for _ in range(n+1)]
    for i in range(1, n+1): dp[i][0]=i*gap; tb[i][0]='u'
    for j in range(1, m+1): dp[0][j]=j*gap; tb[0][j]='l'
    for i in range(1, n+1):
        for j in range(1, m+1):
            sc = match if q[i-1]==s[j-1] else mismatch
            diag=dp[i-1][j-1]+sc; up=dp[i-1][j]+gap; left=dp[i][j-1]+gap
            best = max(diag, up, left)
            dp[i][j]=best; tb[i][j]='d' if best==diag else ('u' if best==up else 'l')
    i, j = n, m
    aq, as_ = [], []
    while i>0 or j>0:
        mv = tb[i][j]
        if mv=='d': aq.append(q[i-1]); as_.append(s[j-1]); i-=1; j-=1
        elif mv=='u': aq.append(q[i-1]); as_.append('-'); i-=1
        else: aq.append('-'); as_.append(s[j-1]); j-=1
    aq.reverse(); as_.reverse()
    return dp[n][m], ''.join(aq), ''.join(as_)

def build_db_from_template(template_seq: str, workdir: Path):
    fa = workdir / "template.fa"
    fa.write_text(">TEMPLATE\n" + template_seq + "\n")
    db_prefix = str(workdir / "template_db")
    _run(["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(fa), "-out", db_prefix])
    return db_prefix

def blast_hits_for(primer_seq: str, db_prefix: str):
    with tempfile.TemporaryDirectory() as td:
        qf = Path(td) / "q.fa"
        qf.write_text(">q\n"+primer_seq+"\n")
        out = _run([
            "blastn","-task","blastn-short","-db",db_prefix,"-query",str(qf),
            "-strand","both","-evalue","30000","-word_size","7",
            "-reward","1","-penalty","-1","-soft_masking","false","-dust","no",
            "-max_target_seqs","50000","-max_hsps","1",
            "-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand"
        ])
    hits=[]
    for line in out.strip().splitlines():
        qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,sstrand = line.split("\t")
        sstart=int(sstart); send=int(send); start,end=(sstart,send) if sstart<=send else (send,sstart)
        hits.append({
            "sseqid":sseqid,"pident":float(pident),"length":int(length),
            "mismatch":int(mismatch),"gapopen":int(gapopen),
            "qstart":int(qstart),"qend":int(qend),
            "sstart":sstart,"send":send,"start":start,"end":end,"sstrand":sstrand
        })
    return hits

def verify_hit(primer_seq: str, hit: dict, template_seq: str, min_perfect_3prime=15,
               nw_mismatch=-1, nw_gap=-1):
    # subject slice (1-based to 0-based)
    s_lo = min(hit["sstart"], hit["send"]) - 1
    s_hi = max(hit["sstart"], hit["send"])
    subj = template_seq[s_lo:s_hi]
    subj_oriented = subj if hit["sstrand"]=="plus" else rc(subj)
    _, aq, as_ = nw_global(primer_seq, subj_oriented, match=1, mismatch=nw_mismatch, gap=nw_gap)
    # 3' perfect window:
    n_ok=0; i=len(aq)-1
    while i>=0 and n_ok<min_perfect_3prime:
        if aq[i] != '-':
            if as_[i]=='-' or aq[i]!=as_[i]: return None
            n_ok+=1
        i-=1
    mism = sum(1 for a,b in zip(aq,as_) if a!='-' and b!='-' and a!=b)
    ident = 1.0 - mism/len(primer_seq)
    return {**hit, "nw_mismatch":mism, "nw_identity":ident}

# ---------- ranking pipeline ----------
def rank_primer_pairs_with_specificity(
    res: dict,
    template_seq: str,
    size_range=(1400,1500),
    min_perfect_3prime=15,
    min_identity=0.65,  # allow up to ~35% mismatches
    weights=None,
    top_k=3
):
    """
    Combine Primer3 quality + BLAST/NW specificity to score each pair.
    Returns a sorted list (best first).
    """
    if weights is None:
        weights = {
            "offtarget_pair": 1000.0,
            "single_primer_amp": 200.0,
            "compl_any": 1.0,      # per unit TH
            "compl_end": 2.0,      # per unit TH
            "left_self_end": 1.0,  # per unit TH
            "right_self_end": 1.0  # per unit TH
        }

    # Prep BLAST DB
    with tempfile.TemporaryDirectory() as td:
        db_prefix = build_db_from_template(template_seq, Path(td))
        # cache BLAST hits per unique primer
        cache = {}

        def get_hits_checked(primer_seq: str, expect_strand: str):
            if primer_seq not in cache:
                cache[primer_seq] = blast_hits_for(primer_seq, db_prefix)
            # keep only hits in expected binding orientation for pairing
            raw = [h for h in cache[primer_seq] if h["sstrand"] == expect_strand]
            ok = []
            for h in raw:
                vh = verify_hit(primer_seq, h, template_seq, min_perfect_3prime=min_perfect_3prime)
                if vh and vh["nw_identity"] >= min_identity:
                    ok.append(vh)
            return ok, cache[primer_seq]  # ok (filtered), raw (for single-primer check)

        results = []
        n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
        for i in range(n_pairs):
            Lseq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            Rseq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            Lpos = res[f"PRIMER_LEFT_{i}"]   # [start, length]
            Rpos = res[f"PRIMER_RIGHT_{i}"]  # [3' end idx, length]
            prod_size = res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE")
            p3_penalty = float(res.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0) or 0.0)

            # BLAST+ + NW verified hits for pairing
            L_ok, L_raw = get_hits_checked(Lseq, expect_strand="plus")
            R_ok, R_raw = get_hits_checked(Rseq, expect_strand="minus")

            # pair amplicons in window
            lo, hi = size_range
            pair_amps=[]
            for L in L_ok:
                for R in R_ok:
                    if L["sseqid"]!=R["sseqid"]: continue
                    start=min(L["start"], R["start"]); end=max(L["end"], R["end"])
                    size=end-start+1
                    if lo<=size<=hi and L["start"]<R["end"]:
                        pair_amps.append({
                            "size":size, "left":L, "right":R
                        })

            # define on-target as any amplicon whose size == Primer3 product_size (±3 bp tolerance)
            on_target = [a for a in pair_amps if prod_size and abs(a["size"]-prod_size)<=3]
            off_target_pairs = max(0, len(pair_amps) - len(on_target))

            # single-primer amplicons (same primer binds + and -)
            def single_amp_count(raw_hits):
                plus=[h for h in raw_hits if h["sstrand"]=="plus"]
                minus=[h for h in raw_hits if h["sstrand"]=="minus"]
                cnt=0
                for P in plus:
                    for M in minus:
                        if P["sseqid"]!=M["sseqid"]: continue
                        start=min(P["start"], M["start"]); end=max(P["end"], M["end"])
                        size=end-start+1
                        if lo<=size<=hi and P["start"]<M["end"]:
                            cnt+=1
                return cnt
            single_left = single_amp_count(L_raw)
            single_right= single_amp_count(R_raw)
            single_total = single_left + single_right

            # add soft penalties for thermodynamic risks (Primer3 metrics)
            compl_any = float(res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0.0) or 0.0)
            compl_end = float(res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0.0) or 0.0)
            left_self_end  = float(res.get(f"PRIMER_LEFT_{i}_SELF_END_TH", 0.0) or 0.0)
            right_self_end = float(res.get(f"PRIMER_RIGHT_{i}_SELF_END_TH", 0.0) or 0.0)

            score = (
                p3_penalty +
                weights["offtarget_pair"] * off_target_pairs +
                weights["single_primer_amp"] * single_total +
                weights["compl_any"] * compl_any +
                weights["compl_end"] * compl_end +
                weights["left_self_end"]  * left_self_end +
                weights["right_self_end"] * right_self_end
            )

            results.append({
                "rank": i,
                "left_seq": Lseq,
                "right_seq": Rseq,
                "product_size": prod_size,
                "primer3_penalty": p3_penalty,
                "offtarget_pairs": off_target_pairs,
                "single_primer_amps": single_total,
                "compl_any_th": compl_any,
                "compl_end_th": compl_end,
                "score": round(score, 3),
                "pair_amps_found": len(pair_amps),
                "on_target_found": len(on_target) > 0
            })

    # sort best first (lowest score wins)
    results.sort(key=lambda r: (r["score"], r["primer3_penalty"]))
    best = results[:top_k]

    # pretty print
    print("Best primer pairs (combined quality + specificity):")
    for j, r in enumerate(best, 1):
        print(f"{j}. P3rank={r['rank']}  score={r['score']}  size={r['product_size']}  "
              f"off-target pairs={r['offtarget_pairs']}  single-amps={r['single_primer_amps']}  "
              f"on-target={'yes' if r['on_target_found'] else 'no'}")
        print(f"   L: {r['left_seq']}")
        print(f"   R: {r['right_seq']}")
    return best


# Main Processing Logic
if __name__ == '__main__':
    if len(sys.argv) > 1:
        gene_ids = sys.argv[1]
    else:
        gene_ids = input("Enter the Gene IDs (comma-separated): ")
    last_exon_info = main(gene_ids)


# Fetch the DNA sequence for the last exon
if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    start = exon_end - 501
    end = exon_end + 351
    
    if exon_end < exon_start:
        start = exon_end - 351
        end = exon_end + 501
        exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    #print(f"New coordinates: Chromosome {chromosome}, Start {start}, End {end}")

    ensembl_server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?coord_system_version=GRCh38"
    
    headers = {"Content-Type": "application/json"}
    r = requests.get(ensembl_server + ext, headers=headers)


    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
    dna_sequence = r.json().get("seq")
    
    # If the start is greater than the end, print the reverse complement of the extended DNA sequence
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_DNA_extended = reverse_complement(dna_sequence)
        #print("Reverse complement of the extended DNA sequence:")
        gdna_sequence = reverse_DNA_extended
    else:
        #print(f"DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        gdna_sequence = dna_sequence

    #print(gdna_sequence)

    num_base_pairs = len(gdna_sequence)

    # Fetch the DNA sequence for the last exon itself
    exon_ext = f"/sequence/region/human/{chromosome}:{exon_start}..{exon_end}:1?coord_system_version=GRCh38"
    r = requests.get(ensembl_server + exon_ext, headers=headers)
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    exon_dna_sequence = r.json().get("seq")
    
    # If the start is greater than the end, print the reverse complement of the exon DNA sequence
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_exon = reverse_complement(exon_dna_sequence)
        exon_seq = Seq(reverse_complement_exon)
    else:
        exon_seq = Seq(exon_dna_sequence)

    num_exon_base_pairs = len(exon_dna_sequence)

    # Check if the length of the exon sequence is a multiple of three
    if num_exon_base_pairs % 3 != 0:
        # Calculate the number of bases to remove
        bases_to_remove = num_exon_base_pairs % 3
        
        # Remove bases from the start of the sequence
        exon_seq = exon_seq[bases_to_remove:]
    
    # Transcribe and translate the exon DNA sequence to amino acids
    amino_acid_seq = exon_seq.translate(to_stop=True)
    #print(f"Amino acid sequence for the last exon ({last_exon_info['exon_id']}):")
    #print(amino_acid_seq)

exon_ext = f"/sequence/region/human/{chromosome}:{exon_start}..{exon_end}:1?coord_system_version=GRCh38"
r = requests.get(ensembl_server + exon_ext, headers=headers)

if not r.ok:
    r.raise_for_status()
    sys.exit()

exon_dna_sequence = r.json().get("seq")

# If the start is greater than the end, print the reverse complement of the extended DNA sequence
if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_exon = reverse_complement(exon_dna_sequence)
        print("Reverse complement of the Exon sequence:")
        last_exon_seq = reverse_complement_exon
else:
        print(f"DNA sequence for the last exon({last_exon_info['exon_id']}):")
        last_exon_seq = exon_dna_sequence

print(last_exon_seq)

# Fetch extended exon sequence (800 bp upstream and downstream)
if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    upstream_start = exon_end - 800
    downstream_end = exon_end + 800
    
    if exon_end < exon_start:
        upstream_start = exon_end - 800
        downstream_end = exon_end + 800
        exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    #print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

    ensembl_server = "https://rest.ensembl.org"
    
    # Fetch the whole sequence including upstream and downstream
    whole_ext = f"/sequence/region/human/{chromosome}:{upstream_start}..{downstream_end}:1?coord_system_version=GRCh38"
    headers = {"Content-Type": "application/json"}
    r_whole = requests.get(ensembl_server + whole_ext, headers=headers)
    
    if not r_whole.ok:
        r_whole.raise_for_status()
        sys.exit()

    whole_dna_sequence_800 = r_whole.json().get("seq")
    
# If the start is greater than the end, print the reverse complement of the extended DNA sequence
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_extended = reverse_complement(whole_dna_sequence_800)
        #print("Reverse complement of the extended DNA sequence +/- 800 bp:")
        exon_seq_800 = reverse_complement_extended
    else:
        #print(f"DNA sequence for the region surrounding the last exon +/- 800 bp({last_exon_info['exon_id']}):")
        exon_seq_800 = whole_dna_sequence_800

#print(exon_seq_800)

# Fetch short upstream and downstream exon sequence (50 bp upstream, 32 bp downstream)
if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    upstream_start = exon_end - 50
    downstream_end = exon_end + 32
    
    if exon_end < exon_start:
        upstream_start = exon_end - 32
        downstream_end = exon_end + 50
        #exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    #print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

    ensembl_server = "https://rest.ensembl.org"
    
    # Fetch the whole sequence including upstream and downstream
    whole_ext = f"/sequence/region/human/{chromosome}:{upstream_start}..{downstream_end}:1?coord_system_version=GRCh38"
    headers = {"Content-Type": "application/json"}
    r_whole = requests.get(ensembl_server + whole_ext, headers=headers)
    
    if not r_whole.ok:
        r_whole.raise_for_status()
        sys.exit()

    whole_dna_sequence = r_whole.json().get("seq")
     
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_CRISPRi = reverse_complement(whole_dna_sequence)
        print(f"Reverse complement of Whole DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        CRISPRtgSearch = reverse_complement_CRISPRi
    else:
        print(f"Whole DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        CRISPRtgSearch = whole_dna_sequence  

print(CRISPRtgSearch)

#PCR Primer Design
if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    upstream_start = exon_end - 1000
    downstream_end = exon_end + 1000

    if exon_end < exon_start:
        upstream_start = exon_end - 1000
        downstream_end = exon_end + 1000
        #exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    #print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

    ensembl_server = "https://rest.ensembl.org"
    
    # Fetch the whole sequence including upstream and downstream
    whole_ext = f"/sequence/region/human/{chromosome}:{upstream_start}..{downstream_end}:1?coord_system_version=GRCh38"
    headers = {"Content-Type": "application/json"}
    r_whole = requests.get(ensembl_server + whole_ext, headers=headers)
    
    if not r_whole.ok:
        r_whole.raise_for_status()
        sys.exit()

    whole_dna_sequence = r_whole.json().get("seq")
     
    if last_exon_info['end'] < last_exon_info['start']:
        rev_PrimerTemplate = reverse_complement(whole_dna_sequence)
        print(f"Reverse complement of Primer Template ({last_exon_info['exon_id']}):")
        PrimerTemplate = rev_PrimerTemplate
    else:
        print(f"Primer Template ({last_exon_info['exon_id']}):")
        PrimerTemplate =  whole_dna_sequence

print(PrimerTemplate)

res = primer3.bindings.design_primers(
    {
        "SEQUENCE_ID": "MY_GENE",
        "SEQUENCE_TEMPLATE": PrimerTemplate,
    },
    {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 62.0,
        "PRIMER_MIN_TM": 60.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 40.0,
        "PRIMER_MAX_GC": 60.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[1400, 2000]],
    },
)

#print(res)

left = res["PRIMER_LEFT_0_SEQUENCE"]
right = res["PRIMER_RIGHT_0_SEQUENCE"]
size = res["PRIMER_PAIR_0_PRODUCT_SIZE"]

# Choose & print only the single best primer pair
best = rank_primer_pairs_with_specificity(
    res,
    template_seq=PrimerTemplate,    # same sequence used for design
    size_range=(1400, 1500),        # keep consistent with Primer3 settings
    min_perfect_3prime=15,
    min_identity=0.65,
    top_k=1
)

if not best:
    print("No suitable primer pair found.")
else:
    b = best[0]              # best scored pair
    i = b["rank"]            # Primer3's rank index

    # Pull more info for this pair from the Primer3 result
    Lseq  = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
    Rseq  = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    Lstart0, Llen = res[f"PRIMER_LEFT_{i}"]         # [start, length] (0-based)
    R3_0,   Rlen  = res[f"PRIMER_RIGHT_{i}"]        # [3' end index, length] (0-based)
    Rstart0 = R3_0 - Rlen + 1

    # 1-based inclusive spans for readability
    Lspan_1b = f"{Lstart0+1}-{Lstart0+Llen}"
    Rspan_1b = f"{Rstart0+1}-{R3_0+1}"

    Ltm   = res.get(f"PRIMER_LEFT_{i}_TM")
    Rtm   = res.get(f"PRIMER_RIGHT_{i}_TM")
    Lgc   = res.get(f"PRIMER_LEFT_{i}_GC_PERCENT")
    Rgc   = res.get(f"PRIMER_RIGHT_{i}_GC_PERCENT")

    print("=== BEST PRIMER PAIR ===")
    print(f"P3 rank      : {i}")
    print(f"Score        : {b['score']}  (lower is better)")
    print(f"Product size : {b['product_size']} bp  |  On-target: {'yes' if b['on_target_found'] else 'no'}")
    print(f"Off-target pairs : {b['offtarget_pairs']}   |   Single-primer amplicons : {b['single_primer_amps']}")
    print()
    print(f"Left  (5'→3'): {Lseq}")
    print(f"  Len/Tm/GC   : {Llen} / {Ltm:.2f} °C / {Lgc:.1f}%")
    print(f"  Position    : {Lspan_1b}  (1-based)")
    print()
    print(f"Right (5'→3'): {Rseq}")
    print(f"  Len/Tm/GC   : {Rlen} / {Rtm:.2f} °C / {Rgc:.1f}%")
    print(f"  Position    : {Rspan_1b}  (1-based)")

# Write sequence to file in the current directory
with open("sequence.fa", "w") as fa_file:
    fa_file.write(f">{gene_ids}_whole\n")
    fa_file.write(CRISPRtgSearch + "\n")

# Define the directory where crispor.py is located
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the directory of this script
crispor_dir = os.path.join(script_dir, "..", "crisporWebsite")  # Move to crisporWebsite directory
crispor_script = os.path.join(crispor_dir, "crispor.py")

# Check if crispor.py exists
if not os.path.exists(crispor_script):
    print(f"Error: crispor.py not found at {crispor_script}")
    sys.exit(1)
else:
    print(f"crispor.py found at {crispor_script}")

# Define the command and arguments for running crispor.py
crispor_command = [
    'python3', crispor_script, 'hg38', 'sequence.fa', 'CRISPRtg.output.scored.tsv'
]

# Run the crispor.py script with real-time output
try:
    print("Running crispor.py...")
    process = subprocess.Popen(
        crispor_command, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        text=True
    )

    # Stream output line by line
    for line in iter(process.stdout.readline, ''):
        print(line, end='')

    process.stdout.close()
    return_code = process.wait()
    
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, crispor_command)

    print("\ncrispor.py executed successfully")

except subprocess.CalledProcessError as e:
    print(f"crispor.py failed with return code {e.returncode}")
    exit(1)

# Define a function to calculate distance from exon for a target sequence
def calculate_distance(row):
    target_seq = row['targetSeq']
    orientation = row['orientation']

    if orientation == "forw":
        # Use the forward strand
        segment = target_seq
    else:
        # Use the reverse complement for reverse orientation
        segment = str(Seq(target_seq).reverse_complement())

    if last_exon_info['end'] < last_exon_info['start']:
        rvs_seq = reverse_complement_exon[-15:]
        last_exon_position = reverse_complement_CRISPRi.rfind(rvs_seq)
        last_letter_position = last_exon_position + len(rvs_seq) - 1
        # Find the position of the segment in the whole DNA sequence
        target_position = reverse_complement_CRISPRi.find(segment)
    else:
        fwd_seq = exon_dna_sequence[-15:]
        last_exon_position = whole_dna_sequence.rfind(fwd_seq)
        last_letter_position = last_exon_position + len(fwd_seq) - 1
        # Find the position of the segment in the whole DNA sequence
        target_position = whole_dna_sequence.find(segment)

    # Calculate the distance from the last letter of the exon sequence
    if orientation == "forw":
        distance = target_position - last_letter_position + 16
    else:
        distance = target_position - last_letter_position + 5

    return distance

try:
    # Load the TSV file
    df = pd.read_csv('CRISPRtg.output.scored.tsv', sep='\t')
    print("DataFrame loaded successfully:")
    print(df.head())


    # Split guideId into location and orientation
    df[['location', 'orientation']] = df['guideId'].str.extract(r'(\d+)(\w+)')


    # Convert location to numeric
    df['location'] = pd.to_numeric(df['location'])


    # Calculate the distance from exon
    df['Distance from Exon'] = df.apply(calculate_distance, axis=1)


    # Format the scores to three significant figures
    score_columns = [
        'mitSpecScore',
        'cfdSpecScore',
        'Doench \'16-Score',
        'Moreno-Mateos-Score',
        'Doench-RuleSet3-Score',
        'Out-of-Frame-Score',
        'Lindel-Score'
    ]


    for col in score_columns:
        df[col] = df[col].apply(lambda x: f"{x:.3g}")


    # Remove any rows where the absolute value of 'Distance from Exon' is greater than 23
    df = df[df['Distance from Exon'].abs() <= 100]

    # Keep only the top 20 entries without sorting
    filtered_df = df.head(20)

    # Define the columns to select
    selected_columns = [
        'targetSeq',
        'location',
        'orientation',
        'mitSpecScore',
        'cfdSpecScore',
        'Moreno-Mateos-Score',
        "Doench '16-Score",
        'Distance from Exon'
    ]

    # Check if the DataFrame is empty after filtering
    if filtered_df.empty:
        print("Error: No valid rows found after processing.")
    else:

        valid_targets = filtered_df[filtered_df["Distance from Exon"] <= 13]
        if valid_targets.empty:
            print("Error: No suitable targets with Distance from Exon <= 13.")
        else:
            # Select the first valid row
            selected_target = valid_targets.iloc[0][['targetSeq', 'orientation', 'Distance from Exon']]

            # Extract selected values
            selected_target_seq = selected_target['targetSeq']
            selected_orientation = selected_target['orientation']
            selected_distance_from_exon = selected_target['Distance from Exon']

            # Print the selected values
            print(f"Selected target sequence: {selected_target_seq}")
            print(f"Selected orientation: {selected_orientation}")
            print(f"Selected distance from exon: {selected_distance_from_exon}")
    
        # Check if the DataFrame is empty after filtering
    if filtered_df.empty:
        print("Error: No valid rows found after processing.")
        selected_scores = {}
    else:
        # Select the first row as the best guide
        selected_target = filtered_df.iloc[0]

        # Extract score columns safely
        score_columns = ['mitSpecScore', 'cfdSpecScore', 'offtargetCount', 'targetGenomeGeneLocus', "Doench '16-Score", 'Moreno-Mateos-Score', 'Doench-RuleSet3-Score', 'Out-of-Frame-Score', 'Lindel-Score', 'GrafEtAlStatus']
        selected_scores = {col: selected_target.get(col, 'NA') for col in score_columns}


except FileNotFoundError:
    print("CRISPRtg.output.scored.tsv file not found.")
except pd.errors.EmptyDataError:
    print("CRISPRtg.output.scored.tsv file is empty.")
except KeyError as e:
    print(f"One or more specified columns are not found in the DataFrame: {e}")
except Exception as e:
    print(f"An error occurred while processing the DataFrame: {e}")

# Function to find and extract sequences around the target
def find_and_extract_sequence(gdna, target, upstream, downstream):
    target_index = gdna.find(target)
    if target_index == -1:
        return None, None
    upstream_sequence = gdna[max(0, target_index + len(target) - upstream):target_index + len(target)]
    downstream_sequence = gdna[target_index + len(target) - 1:target_index + len(target) - 1 + downstream]
    return upstream_sequence, downstream_sequence

# Function to translate a DNA sequence into an amino acid sequence
def translate_dna_to_protein(dna_sequence):
    dna_seq = Seq(dna_sequence)
    return str(dna_seq.translate())

def find_sgrna_coordinates_in_exon(exon_sequence, sgrna_sequence):
    # Find the start index of the sgRNA sequence within the exon sequence
    start_index = exon_sequence.find(sgrna_sequence)
    if start_index == -1:
        return None, None
    end_index = start_index + len(sgrna_sequence)
    return start_index, end_index

def trim_to_multiple_of_three(sequence):
    """Trims the sequence from the start until its length is a multiple of three."""
    start_index = len(sequence) % 3
    return sequence[start_index:]

def codon_order(sequence):
    """Returns the codon order (123123123...) for the DNA sequence."""
    codon_indices = [(i % 3) + 1 for i in range(len(sequence))]
    return ''.join(map(str, codon_indices))

def check_sgrna_coordinates(sequence, sgrna_start, sgrna_end, codon_order_sequence):
    """Checks if the last two bases of the sgRNA (coordinates sgrna_end - 2 and sgrna_end - 1) are in codon positions 2 or 3."""
    if sgrna_start is None and sgrna_end is None:
        print("Error: sgRNA not found in exon sequence.")
        return None
    
    sgrna_extracted = sequence[sgrna_start:sgrna_end]
    #print(f"Extracted sgRNA Sequence from Exon: {sgrna_extracted}")

    if is_sgRNA_inverted == "Y":
        print(f"Checking sgRNA coordinates: start={sgrna_start}, end={sgrna_end}")
        position_1 = codon_order_sequence[sgrna_start]
        position_2 = codon_order_sequence[sgrna_start + 1]
        base_2 = sequence[sgrna_start + 1]
        base_3 = sequence[sgrna_start + 2]
        base_1 = sequence[sgrna_start]
        #print(f"Codon positions for first two bases: position_2={position_2}, position_1={position_1}")
        #print(f"Base pairs at these positions:  base_1={base_1}, base_2={base_2},base_3{base_3}")
        if position_1 == '3' or position_2 == '3':
            return "One or both of the last two bases of the sgRNA are in codon position 3."
        else:
            return "Neither of the last two bases of the sgRNA are in codon positions 3."
    else:   
        print(f"Checking sgRNA coordinates: start={sgrna_start}, end={sgrna_end}")
        position_2 = codon_order_sequence[sgrna_end - 2]
        position_3 = codon_order_sequence[sgrna_end - 1]
        base_2 = sequence[sgrna_end - 2]
        base_3 = sequence[sgrna_end - 1]
        base_1 = sequence[sgrna_end - 3]
        #print(f"Codon positions for last two bases: position_2={position_2}, position_3={position_3}")
        #print(f"Base pairs at these positions: base_2={base_2}, base_1={base_1},base_3{base_3}")
        if position_2 == '3' or position_3 == '3':
            return "One or both of the last two bases of the sgRNA are in codon position 3."
        else:
            return "Neither of the last two bases of the sgRNA are in codon positions 3."
    
    return "sgRNA coordinates are invalid or sgRNA not found in the sequence."

def get_alternative_codon(original_codon, amino_acid):
    """Returns an alternative codon for the given amino acid that is different from the original codon."""
    standard_table = CodonTable.unambiguous_dna_by_id[1]  # Standard codon table
    synonymous_codons = [codon for codon, aa in standard_table.forward_table.items() if aa == amino_acid]
    synonymous_codons = [codon for codon in synonymous_codons if codon != original_codon]
    return synonymous_codons[0] if synonymous_codons else None

def check_non_overlapping_codons(sequence, sgrna_start, sgrna_end, codon_order_sequence):
    """
    Checks non-overlapping codons in the sgRNA sequence and verifies their
    positions in the codon order (123123...).
    """
    if sgrna_start is not None and sgrna_end is not None:
        #print(f"Checking non-overlapping codons in sgRNA sequence: start={sgrna_start}, end={sgrna_end}")
        
        # Iterate through sgRNA sequence in non-overlapping codons
        for i in range(sgrna_start, sgrna_end, 3):
            codon = sequence[i:i + 3]
            codon_positions = codon_order_sequence[i:i + 3]
            
            # Ensure codon is complete (length of 3)
            if len(codon) == 3:
                print(f"Codon: {codon} at positions {codon_positions}")
                
                # Determine overlap with position 3
                if '3' in codon_positions:
                    print(f"Codon {codon} overlaps with position 3 in codon order.")
                else:
                    print(f"Codon {codon} does not overlap with position 3 in codon order.")
            else:
                print(f"Incomplete codon at the end: {codon}")
    else:
        print("sgRNA coordinates are invalid or sgRNA not found in the sequence.")

def find_partial_sgrna_coordinates_in_exon(exon_sequence, sgrna_sequence, min_match_length):
    for i in range(len(sgrna_sequence) - min_match_length + 1):
        partial_sequence = sgrna_sequence[i:i + min_match_length]
        start_index = exon_sequence.find(partial_sequence)
        if start_index != -1:
            return start_index, start_index + min_match_length
    return None, None

def modify_sgrna_codons(dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence):
    """
    Modify the sgRNA sequence codon-by-codon within the left_arm while preserving the amino acid sequence.
    
    Args:
        dna_sequence_trimmed (str): Trimmed DNA sequence.
        sgrna_start (int): Start position of sgRNA.
        sgrna_end (int): End position of sgRNA.
        codon_order_sequence (list): List of preferred codon order.
        left_arm (str): The left arm sequence.
        sgRNA_sequence (str): The sgRNA sequence.
    
    Returns:
        tuple: Modified left arm sequence and verification status.
    """
    # Locate the sgRNA sequence in the left_arm
    sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)

    # Check if sgRNA is found in the left_arm
    if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
        print(f"sgRNA sequence found in left_arm at positions: start={sgrna_in_left_arm_start}, end={sgrna_in_left_arm_end}")
    else:
        print("sgRNA sequence not found in the left_arm.")
        return left_arm, False
    
    # Extract the sgRNA sequence from the left_arm
    sgRNA_in_left_arm = left_arm[sgrna_in_left_arm_start:sgrna_in_left_arm_end]
    
    # Modify the sgRNA sequence codon-by-codon
    modified_left_arm = left_arm
    for i in range(0, len(sgRNA_in_left_arm), 3):
        codon_start = sgrna_in_left_arm_start + i
        codon_end = codon_start + 3
        original_codon = left_arm[codon_start:codon_end]
        
        if len(original_codon) == 3:  # Ensure the codon is complete
            # Translate the codon to an amino acid
            amino_acid = translate_dna_to_protein(original_codon)
            
            # Get an alternative codon
            alternative_codon = get_alternative_codon(original_codon, amino_acid)
            if alternative_codon:
                # Replace the codon in the left_arm sequence
                modified_left_arm = (
                    modified_left_arm[:codon_start] +
                    alternative_codon +
                    modified_left_arm[codon_end:]
                )
            else:
                print(f"No alternative codon available for {original_codon} ({amino_acid}).")
    
    # Trim sequences to multiples of three
    trimmed_left_arm = trim_to_multiple_of_three(left_arm)
    trimmed_edited_left_arm = trim_to_multiple_of_three(modified_left_arm)
    
    # Verify the amino acid sequence remains the same
    original_left_arm_aa = translate_dna_to_protein(trimmed_left_arm)
    edited_left_arm_aa = translate_dna_to_protein(trimmed_edited_left_arm)
    
    #print(f"Original left arm amino acid sequence: {original_left_arm_aa}")
    #print(f"Edited left arm amino acid sequence: {edited_left_arm_aa}")
    
    if original_left_arm_aa == edited_left_arm_aa:
        print("Verification successful: The amino acid sequence remains unchanged.")
        return modified_left_arm, True
    else:
        print("Verification failed: The amino acid sequence has changed.")
        return modified_left_arm, False

def modify_sgrna_PAM(left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence):
    """
    Modify the last three bases of the sgRNA sequence within the left arm while preserving the amino acid sequence.
    
    Args:
        left_arm (str): The left arm sequence containing sgRNA.
        last_three_start (int): Start position of the last three bases in sgRNA.
        last_three_bases (str): The last three bases of the sgRNA sequence.
        dna_sequence_trimmed (str): Trimmed DNA sequence.
        sgrna_start (int): Start position of sgRNA.
        sgrna_end (int): End position of sgRNA.
        codon_order_sequence (list): List of preferred codon order.
        sgRNA_sequence (str): The sgRNA sequence.
    
    Returns:
        tuple: Edited left arm sequence and verification status.
    """
    amino_acid = translate_dna_to_protein(last_three_bases)
    alternative_codon = get_alternative_codon(last_three_bases, amino_acid)
    
    if alternative_codon:
        edited_left_arm = (
            left_arm[:last_three_start] + alternative_codon + left_arm[last_three_start+3:]
        )
        
        trimmed_left_arm = trim_to_multiple_of_three(left_arm)
        trimmed_edited_left_arm = trim_to_multiple_of_three(edited_left_arm)
        
        original_left_arm_aa = translate_dna_to_protein(trimmed_left_arm)
        edited_left_arm_aa = translate_dna_to_protein(trimmed_edited_left_arm)
        
        if original_left_arm_aa == edited_left_arm_aa:
            print("Verification successful: The amino acid sequence remains unchanged.")
            return edited_left_arm, True
        else:
            print("Verification failed: The amino acid sequence has changed.")
            modified_left_arm, verification_status = modify_sgrna_codons(
                dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
            )
            return modified_left_arm, verification_status
    
    return left_arm, False

gene_name = gene_ids
sgRNA_sequence = selected_target_seq

if selected_orientation == "rev":
    is_sgRNA_inverted = "Y"
else:
    is_sgRNA_inverted = "N"

original_sgRNA = sgRNA_sequence

if len(last_exon_seq) >= 20:
    protein_coding_exon = last_exon_seq[-20:]
else:
    protein_coding_exon = last_exon_seq


if is_sgRNA_inverted == "Y":
    sgRNA_core = reverse_complement(original_sgRNA)[3:]  # Exclude first 3 bases (PAM)

else:
    sgRNA_core = original_sgRNA[:-3]  # Exclude last 3 bases (PAM)


left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)
#print(left_arm)

# Define the minimum match length as half of the sgRNA sequence length excluding PAM
min_match_length = len(sgRNA_sequence) // 2

def find_partial_sgrna_coordinates_in_exon(exon_sequence, sgrna_sequence, min_match_length):
    for i in range(len(sgrna_sequence) - min_match_length + 1):
        partial_sequence = sgrna_sequence[i:i + min_match_length]
        start_index = exon_sequence.find(partial_sequence)
        if start_index != -1:
            return start_index, start_index + min_match_length
    return None, None

# Modify the function to find at least half of the sgRNA sequence in left_arm
sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_partial_sgrna_coordinates_in_exon(left_arm, sgRNA_core, min_match_length)
#print(sgrna_in_left_arm_start, sgrna_in_left_arm_end)


if selected_distance_from_exon >= -3 or (selected_distance_from_exon <= -3 and sgrna_in_left_arm_end is None):
    # Pathway 1: Standard processing
    
    #print("sgRNA sequence is not found in the left_arm sequence.")

    # Add logic for handling this case
    #print(f"Last Exon Sequence{last_exon_seq}")
    #exon_amino_acid_sequence = translate_dna_to_protein(last_exon_seq)
    #print(f"Amino acid sequence of last exon: {'-'.join(exon_amino_acid_sequence)}")

    # Step 1: Print the last exon sequence
    #print(f"Area around Last exon sequence: {gdna_sequence}")

    # Step 2: Translate the last exon sequence to amino acids
    #amino_acid_sequence = translate_dna_to_protein(gdna_sequence)
    #print(f"Amino acid sequence of area around last exon: {'-'.join(amino_acid_sequence)}")

    # Original pathway for non-negative distance from exon
    pass

else:
    
    print("Handling negative distance from exon...")

    
    if is_sgRNA_inverted == "Y":
        #print(f"main sgRNA_Sequence {original_sgRNA}")
        #print(f"no PAM sgRNA_Sequence {sgRNA_core}")
        sgRNA_sequence = reverse_complement(original_sgRNA[-3:]) + sgRNA_core   # Reattach PAM at the end
        
    else:
        #print(f"main sgRNA_Sequence {original_sgRNA}")
        #print(f"no PAM sgRNA_Sequence {sgRNA_core}")
        sgRNA_sequence = sgRNA_core + original_sgRNA[-3:]

    print(f"sgRNA_Sequence {sgRNA_sequence}")
    # Trim DNA sequence to a multiple of three
    dna_sequence_trimmed = trim_to_multiple_of_three(last_exon_seq)
    print(f"Last Exon Sequence {dna_sequence_trimmed}")

    # Translate the trimmed sequence
    protein_sequence = translate_dna_to_protein(dna_sequence_trimmed)
    #print(f"Last Exon AA Sequence {protein_sequence}")

    # Generate codon order for the trimmed sequence
    codon_order_sequence = codon_order(dna_sequence_trimmed)

    # Find sgRNA in the trimmed sequence
    sgrna_start, sgrna_end = find_partial_sgrna_coordinates_in_exon(dna_sequence_trimmed, sgRNA_sequence, min_match_length)
    print(f"sgRNA Coord {sgrna_start} {sgrna_end}")


    # Check sgRNA coordinates in codon order
    sgrna_coordinates_check = check_sgrna_coordinates(dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence)

    #print(f"Trimmed DNA Sequence: {dna_sequence_trimmed}")
    #print(f"Protein Sequence: {protein_sequence}")
    #print(f"Codon Order: {codon_order_sequence}")
    
    if sgrna_start is not None:
        print("Testing for modifying PAM")
        print(f"sgRNA found at positions {sgrna_start} to {sgrna_end}")
        sgrna_sequence_found = dna_sequence_trimmed[sgrna_start:sgrna_end]
        left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)
        print(f"Base pairs in the sgRNA positions: {sgrna_sequence_found}")
        

        print(f"is_sgRNA_inverted: {is_sgRNA_inverted}")
        if is_sgRNA_inverted == "Y":
            position_1 = codon_order_sequence[sgrna_start]
            position_2 = codon_order_sequence[sgrna_start + 1]

            if position_2 == '3' or position_1 == '3':
            
                print("Pathway: One or both of the last two bases are in codon position 3.")
                left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)

                # Original pathway for non-negative distance from exon
                left_arm_start, left_arm_end = find_sgrna_coordinates_in_exon(gdna_sequence, left_arm)
                print(f"Left arm coordinates in exon: start={left_arm_start}, end={left_arm_end}")

                sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_partial_sgrna_coordinates_in_exon(left_arm, sgRNA_core, min_match_length)
                print(sgrna_in_left_arm_start, sgrna_in_left_arm_end)

                if position_2 == '3':                        
                        print("Attempting to Edit PAM")
                        print(f"Current left arm sequence: {left_arm}")
                        rev_left_arm = reverse_complement(left_arm)
                        print(f"rev comp of left arm {rev_left_arm}")
                        print(f"{sgrna_in_left_arm_start},{sgrna_in_left_arm_end}")

                            # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                        if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                            # Get the last three bases of the sgRNA sequence in the left arm
                            last_three_start = sgrna_in_left_arm_start-4
                            last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_start-1]
                            print(f"Last codon of sgRNA in left arm: {last_three_bases}")

                            if last_three_bases == "TGG":
                                modified_left_arm, verification_status = modify_sgrna_codons(
                                dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                                )
                            else:
                                modified_left_arm, verification_status = modify_sgrna_PAM(
                                left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                                sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                                )
                        else:
                            print("sgRNA sequence not found in the left arm.")
                        
                else:
                    print("Attempting to Edit PAM")
                    print(f"Current left arm sequence: {left_arm}")


                        # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                    if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                        # Get the last three bases of the sgRNA sequence in the left arm
                        last_three_start = sgrna_in_left_arm_start - 5
                        last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_start - 2]
                        print(f"Last three bases of sgRNA in left arm: {last_three_bases}")

                        if last_three_bases == "TGG":
                            modified_left_arm, verification_status = modify_sgrna_codons(
                                dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                                )
                        else:
                                modified_left_arm, verification_status = modify_sgrna_PAM(
                        left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                        sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                        )
                        
            else:
                #Start Here
                print("Pathway: Neither of the last two bases is in codon position 3.")
                check_non_overlapping_codons(dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence)
                    
                # Locate the sgRNA sequence in the left_arm
                sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_partial_sgrna_coordinates_in_exon(left_arm, sgRNA_core, min_match_length)
                print(sgrna_in_left_arm_start, sgrna_in_left_arm_end)

                modified_left_arm, verification_status = modify_sgrna_codons(
                dna_sequence_trimmed, sgrna_in_left_arm_start, sgrna_in_left_arm_end, codon_order_sequence, left_arm, sgRNA_sequence
                )            
        else:    
            # Pathway based on codon positions
            position_2 = codon_order_sequence[sgrna_end - 2]
            position_3 = codon_order_sequence[sgrna_end - 1]
            print(position_2,position_3)
            if position_2 == '3' or position_3 == '3':

                print("Pathway: One of the last two bases are in codon position 3.")
                left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)

                # Original pathway for non-negative distance from exon
                left_arm_start, left_arm_end = find_sgrna_coordinates_in_exon(gdna_sequence, left_arm)
                print(f"Left arm coordinates in exon: start={left_arm_start}, end={left_arm_end}")

                # Find sgRNA amino acids within the left arm sequence
                sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
                print(sgrna_in_left_arm_start,sgrna_in_left_arm_end)
                
                # Extract amino acid sequence within the left arm coordinates
                if left_arm_start is not None and left_arm_end is not None:
                    left_arm_amino_acids = translate_dna_to_protein(dna_sequence_trimmed[left_arm_start:left_arm_end])
                    print(f"Amino acids for left arm: {left_arm_amino_acids}")
                    if position_2 == '3':

                        # Find sgRNA amino acids within the left arm sequence
                        sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
                        print(sgrna_in_left_arm_start,sgrna_in_left_arm_end)

                        print(f"Current left arm sequence: {left_arm}")

                            # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                        if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                            # Get the last three bases of the sgRNA sequence in the left arm
                            last_three_start = sgrna_in_left_arm_end - 4
                            last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_end-1]
                            print(f"Last codon of sgRNA in left arm: {last_three_bases}")

                            if last_three_bases == "TGG":
                                modified_left_arm, verification_status = modify_sgrna_codons(
                                    dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                                )
                            else:
                                modified_left_arm, verification_status = modify_sgrna_PAM(
                                    left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                                    sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                                )
                        else:
                            print("sgRNA sequence not found in the left arm.")
                        
                    else:
                        # Find sgRNA amino acids within the left arm sequence
                        sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
                        print(sgrna_in_left_arm_start,sgrna_in_left_arm_end)

                        print(f"Current left arm sequence: {left_arm}")

                            # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                        if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                            # Get the last three bases of the sgRNA sequence in the left arm
                            last_three_start = sgrna_in_left_arm_end - 3
                            last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_end]
                            print(f"Last three bases of sgRNA in left arm: {last_three_bases}")

                            if last_three_bases == "TGG":
                                modified_left_arm, verification_status = modify_sgrna_codons(
                                    dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                                )
                            else:
                                modified_left_arm, verification_status = modify_sgrna_PAM(
                                    left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                                    sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                                )
                else:
                    print("sgRNA sequence not found in the left arm.")
            else:
                print("Pathway: Neither of the last two bases is in codon position 3.")

                modified_left_arm, verification_status = modify_sgrna_codons(
                    dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                )

    else:
        print("Pathway: Neither of the last two bases is in codon position 3.")

        modified_left_arm, verification_status = modify_sgrna_codons(
            dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
        )
 

# Use modified_left_arm if available, otherwise fallback to left_arm
left_arm_seq = modified_left_arm if 'modified_left_arm' in locals() else left_arm

print(f"sgRNA_sequence {sgRNA_sequence}" )

if is_sgRNA_inverted == "Y" and selected_orientation == "rev":
    rev_sgRNA_sequence = reverse_complement(sgRNA_sequence)
    print(f"Rev Compliment of sgRNA {rev_sgRNA_sequence}")

    # Try RC (7 bp) first
    _, right_arm = find_and_extract_sequence(gdna_sequence, rev_sgRNA_sequence[:7], 0, 321)
    print(f"Right_arm sequence (RC 7bp):\n{right_arm}")

    # Fallback: forward (18 bp) if no match
    if right_arm is None:
        _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
        print(f"Right_arm sequence (FWD 7bp fallback):\n{right_arm}")

elif is_sgRNA_inverted == "Y":
    # Keep your original behavior for the non-rev case
    _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)
    print(f"Right_arm sequence:\n{right_arm}")

else:
    # Not inverted: original behavior
    _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)
    print(f"Right_arm sequence:\n{right_arm}")



left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
if is_sgRNA_inverted == "Y":
    left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", left_arm_seq)
else:
    left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", left_arm_seq)
print(f">{gene_name}-Left-AI1\n{left_ai1_seq}")

right_ai1_template = "CCTGCGGTGTCTTTGCTTrycatgtGGTTCCATGGTGTAATGGTTAGCACTCTGGACTCTGAATCCAGCGATCCGAGTTCAAATCTCGGTGGAACCTxGTTTTAGAGCTAGAAATAGCAA"
if is_sgRNA_inverted == "Y":
    right_ai1_seq = (
        right_ai1_template.replace("r", right_arm)
        .replace("y", original_sgRNA)
        .replace("x", original_sgRNA[:20])
    )
else:
    right_ai1_seq = (
        right_ai1_template.replace("r", right_arm)
        .replace("y", original_sgRNA)
        .replace("x", original_sgRNA[:20])
    )


# Print sequences
print(f">{gene_name}-Left-AI1\n{left_ai1_seq}")
print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
print(f"Right_arm sequence:\n{right_arm}")
print(f"Left_arm sequence:\n{left_arm_seq}")

# Ensure the output directory exists
output_dir = "HDR_arms"
os.makedirs(output_dir, exist_ok=True)

# ---- helper: extract best pair either from our earlier ranking function or fall back to P3 rank 0 ----
def get_best_primer_from_res(res: dict):
    """
    Returns (index, summary_dict) for the 'best' pair.
    If a ranking function named rank_primer_pairs_with_specificity exists, uses that (top_k=1).
    Otherwise falls back to Primer3 rank 0.
    """
    try:
        # Use your BLAST+ specificity ranking if it's defined elsewhere in your script
        best = rank_primer_pairs_with_specificity(
            res,
            template_seq=PrimerTemplate,
            size_range=(1400, 1500),   # keep consistent with your Primer3 settings
            min_perfect_3prime=15,
            min_identity=0.65,
            top_k=1
        )
        if best:
            b = best[0]
            i = b["rank"]
            summary = {
                "score": b["score"],
                "product_size": b["product_size"],
                "on_target": b["on_target_found"],
                "offtarget_pairs": b["offtarget_pairs"],
                "single_primer_amplicons": b["single_primer_amps"],
            }
            return i, summary
    except NameError:
        # ranking function not present — fall through to P3 rank 0
        pass

    # Fallback: use Primer3 rank 0 only
    if int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0) == 0:
        return None, None
    i = 0
    summary = {
        "score": float(res.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0) or 0.0),
        "product_size": res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
        "on_target": True,      # unknown without BLAST; mark True for fallback
        "offtarget_pairs": None,
        "single_primer_amplicons": None,
    }
    return i, summary

def _fmt(x, nd=2):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return "NA" if x is None else str(x)

def _span_1b_left(res, i):
    try:
        start0, length = res[f"PRIMER_LEFT_{i}"]
        return f"{start0+1}-{start0+length}", length
    except Exception:
        return "NA", None

def _span_1b_right(res, i):
    try:
        end3_0, length = res[f"PRIMER_RIGHT_{i}"]
        start0 = end3_0 - length + 1
        return f"{start0+1}-{end3_0+1}", length
    except Exception:
        return "NA", None

# Check if the DataFrame is empty after filtering
if filtered_df.empty:
    print("Error: No valid rows found after processing.")
    selected_scores = {}
else:

    valid_targets = filtered_df[filtered_df["Distance from Exon"] <= 13].copy()
    if valid_targets.empty:
        print("Error: No suitable targets with Distance from Exon <= 13.")
        selected_scores = {}
    else:
        # Select the first valid row
        selected_target = valid_targets.iloc[0].copy()

        # Extract key values
        selected_target_seq = selected_target.get('targetSeq', 'NA')
        selected_orientation = selected_target.get('orientation', 'NA')
        selected_distance_from_exon = selected_target.get('Distance from Exon', 'NA')

        # Score columns to display
        score_columns = [
            'mitSpecScore',
            'cfdSpecScore',
            'offtargetCount',
            'targetGenomeGeneLocus',
            'Doench 2016-Score',
            'Moreno-Mateos-Score',
            'Doench-RuleSet3-Score',
            'Out-of-Frame-Score',
            'Lindel-Score',
            'GrafEtAlStatus'
        ]
        selected_scores = {col: selected_target.get(col, 'NA') for col in score_columns}

    # ---- pull best primer pair info (Primer3 + optional BLAST specificity) ----
    best_idx, best_summary = (None, None)
    try:
        best_idx, best_summary = get_best_primer_from_res(res)
    except Exception as e:
        print(f"[WARN] Primer best-pair extraction failed: {e}")

    # Build the single TXT output
    lines = [
        f"Selected target sequence: {selected_target_seq}",
        f"Selected orientation: {selected_orientation}",
        f"Selected distance from exon: {selected_distance_from_exon}",
        "",
        "Score summary:"
    ]
    lines += [f"- {col}: {selected_scores[col]}" for col in score_columns]

    # HDR arms (FASTA-style headers kept as you had them)
    lines += [
        "",
        f">{gene_name}-Left-AI1",
        left_ai1_seq,
        f">{gene_name}-Right-AI1",
        right_ai1_seq,
        "",
        "Left_arm sequence:",
        left_arm_seq,
        "",
        "Right_arm sequence:",
        right_arm
    ]

    # ===== ADD: Best primer pair directly in this same file (FASTA + summary) =====
    lines += ["", "===== Best primer pair ====="]
    if best_idx is None:
        lines += ["No primer pairs returned by Primer3."]
    else:
        Lseq = res.get(f"PRIMER_LEFT_{best_idx}_SEQUENCE", "NA")
        Rseq = res.get(f"PRIMER_RIGHT_{best_idx}_SEQUENCE", "NA")
        Lspan, Llen = _span_1b_left(res, best_idx)
        Rspan, Rlen = _span_1b_right(res, best_idx)
        Ltm  = _fmt(res.get(f"PRIMER_LEFT_{best_idx}_TM"))
        Rtm  = _fmt(res.get(f"PRIMER_RIGHT_{best_idx}_TM"))
        Lgc  = _fmt(res.get(f"PRIMER_LEFT_{best_idx}_GC_PERCENT"), 1)
        Rgc  = _fmt(res.get(f"PRIMER_RIGHT_{best_idx}_GC_PERCENT"), 1)
        prod = res.get(f"PRIMER_PAIR_{best_idx}_PRODUCT_SIZE", "NA")
        ppen = _fmt(res.get(f"PRIMER_PAIR_{best_idx}_PENALTY"))

        # FASTA-style primer sequences (so they can be copy-pasted for ordering)
        lines += [
            f">{gene_name}-BestPrimer-LEFT",
            Lseq,
            f">{gene_name}-BestPrimer-RIGHT",
            Rseq,
            "",
            "Best primer summary:",
            f"- Primer3 rank index: {best_idx}",
            f"- Product size: {prod} bp   |   Primer3 pair penalty: {ppen}",
            f"- Left  len/Tm/GC: {Llen} / {Ltm} °C / {Lgc}%   |   Pos (1-based): {Lspan}",
            f"- Right len/Tm/GC: {Rlen} / {Rtm} °C / {Rgc}%   |   Pos (1-based): {Rspan}",
        ]

        if isinstance(best_summary, dict):
            lines += [
                f"- Specificity combined score: {best_summary.get('score')}",
                f"- On-target found: {best_summary.get('on_target')}",
                f"- Off-target pairs: {best_summary.get('offtarget_pairs')}",
                f"- Single-primer amplicons: {best_summary.get('single_primer_amplicons')}",
            ]

    # ===== Keep full list of candidates in the same file =====
    lines += ["", "===== All Primer3 candidate pairs ====="]
    n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
    if n_pairs == 0:
        lines += ["(none)"]
    else:
        for i in range(n_pairs):
            Lseq = res.get(f"PRIMER_LEFT_{i}_SEQUENCE", "NA")
            Rseq = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "NA")
            Lspan, Llen = _span_1b_left(res, i)
            Rspan, Rlen = _span_1b_right(res, i)
            Ltm  = _fmt(res.get(f"PRIMER_LEFT_{i}_TM"))
            Rtm  = _fmt(res.get(f"PRIMER_RIGHT_{i}_TM"))
            Lgc  = _fmt(res.get(f"PRIMER_LEFT_{i}_GC_PERCENT"), 1)
            Rgc  = _fmt(res.get(f"PRIMER_RIGHT_{i}_GC_PERCENT"), 1)
            prod = res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", "NA")
            ppen = _fmt(res.get(f"PRIMER_PAIR_{i}_PENALTY"))

            lines += [
                f"[Pair {i}]  product={prod} bp  |  penalty={ppen}",
                f"  L: {Lseq}",
                f"     len/tm/gc: {Llen} / {Ltm} °C / {Lgc}%   pos: {Lspan}",
                f"  R: {Rseq}",
                f"     len/tm/gc: {Rlen} / {Rtm} °C / {Rgc}%   pos: {Rspan}",
                ""
            ]

    output_text = "\n".join(lines)

    # Define the output path (single TXT)
    txt_path  = os.path.join(output_dir, f"{gene_name}_selected_info.txt")

    # Write the single human-readable TXT (arms + best primer + full list)
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write(output_text)

    print(f"\nAll results written to: {txt_path}")
