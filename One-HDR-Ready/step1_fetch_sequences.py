#!/usr/bin/env python3
"""
Step 1 — Fetch sequences for One-HDR-Ready

Inputs
------
- Gene identifier string: either GeneCards gene symbol/name OR UniProt accession(s).

What this step does
-------------------
1) Get UniProt accession (directly or via mapping).
2) Query EBI Proteins API for exon coordinates and identify the *last exon*.
3) Use Ensembl REST API to fetch:
   - Last exon sequence
   - CRISPOR input sequence (50 bp upstream, 35 bp downstream)
   - Guide DNA window (501 bp upstream, 350 bp downstream)
   - Primer template (±1000 bp)
4) Saves outputs individually as FASTA and JSON files in ./results:
   - step1_exon.fa
   - step1_crispor_input.fa
   - step1_guide_window.fa
   - step1_primer_template.fa
   - step1_sequences.json
5) Returns a dictionary when imported and called by main.py.

All sequences are reverse-complemented if the gene is on the reverse strand so
that all outputs are in 5'→3' gene-sense orientation.
"""

from __future__ import annotations
import os, sys, json, time, re
from pathlib import Path
from typing import Dict, Any, List
import requests
from Bio.Seq import Seq

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
ENSEMBL_SERVER = "https://rest.ensembl.org"
COORD_SYS_VERSION = "GRCh38"

def reverse_complement(seq: str) -> str:
    t = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(t)[::-1]

def is_uniprot_accession(code: str) -> bool:
    if not code:
        return False
    base = code.strip().upper().split("-")[0]
    return bool(re.fullmatch(r"[A-Z0-9]{6}", base)) or bool(re.fullmatch(r"[A-Z0-9]{10}", base))

def uniprot_initiate_id_mapping(ids: str, from_db: str, to_db: str) -> str:
    url = "https://rest.uniprot.org/idmapping/run"
    r = requests.post(url, data={"ids": ids, "from": from_db, "to": to_db})
    r.raise_for_status()
    return r.json()["jobId"]

def uniprot_check_id_mapping(job_id: str) -> dict:
    url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()

def ebi_proteins_coordinates(uniprot_acc: str) -> dict:
    url = f"https://www.ebi.ac.uk/proteins/api/coordinates/{uniprot_acc}"
    r = requests.get(url, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

def ensembl_sequence_region(chrom: str, start: int, end: int, strand: int) -> str:
    ext = f"/sequence/region/human/{chrom}:{start}..{end}:{strand}?coord_system_version={COORD_SYS_VERSION}"
    r = requests.get(ENSEMBL_SERVER + ext, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json().get("seq", "")

def pick_last_exon(ebi_json: dict) -> Dict[str, Any]:
    valid_chr = {str(i) for i in range(1, 23)} | {"X", "Y"}
    infos: List[dict] = []
    for gene in ebi_json.get("gnCoordinate", []):
        ensembl_gene_id = gene.get("ensemblGeneId")
        loc = gene.get("genomicLocation", {})
        chrom = str(loc.get("chromosome"))
        if chrom not in valid_chr:
            continue
        reverse_strand = bool(loc.get("reverseStrand", False))
        for exon in loc.get("exon", []):
            start = exon.get("genomeLocation", {}).get("begin", {}).get("position")
            end = exon.get("genomeLocation", {}).get("end", {}).get("position")
            if start is None or end is None:
                continue
            infos.append({
                "ensembl_gene_id": ensembl_gene_id,
                "exon_id": exon.get("id"),
                "chromosome": chrom,
                "start": int(start),
                "end": int(end),
                "reverseStrand": reverse_strand,
            })
    if not infos:
        raise RuntimeError("No valid exon coordinates found on autosomes/X/Y.")
    reverse = infos[0]["reverseStrand"]
    last = min(infos, key=lambda x: x["start"]) if reverse else max(infos, key=lambda x: x["end"])
    return last

def compute_windows(exon: Dict[str, Any]) -> Dict[str, Any]:
    chrom = exon["chromosome"]
    start = int(exon["start"])
    end = int(exon["end"])
    reverse = bool(exon["reverseStrand"])
    stop_edge = start if reverse else end
    strand_arg = 1

    crispor_up, crispor_dn = stop_edge - 50, stop_edge + 35
    guide_up, guide_dn = stop_edge - 501, stop_edge + 350
    primer_up, primer_dn = stop_edge - 1000, stop_edge + 1000
    exon_lo, exon_hi = min(start, end), max(start, end)

    seq_exon = ensembl_sequence_region(chrom, exon_lo, exon_hi, strand_arg)
    seq_crispor = ensembl_sequence_region(chrom, crispor_up, crispor_dn, strand_arg)
    seq_guide = ensembl_sequence_region(chrom, guide_up, guide_dn, strand_arg)
    seq_primer = ensembl_sequence_region(chrom, primer_up, primer_dn, strand_arg)

    if reverse:
        seq_exon = reverse_complement(seq_exon)
        seq_crispor = reverse_complement(seq_crispor)
        seq_guide = reverse_complement(seq_guide)
        seq_primer = reverse_complement(seq_primer)

    return {
        "chromosome": chrom,
        "exon_start": start,
        "exon_end": end,
        "reverseStrand": reverse,
        "stop_edge": stop_edge,
        "seq_exon": seq_exon,
        "seq_crispor_input": seq_crispor,
        "seq_guide_window": seq_guide,
        "seq_primer_template": seq_primer,
    }

def run(gene_ids: str) -> Dict[str, Any]:
    gene_list = [g.strip() for g in gene_ids.split(",") if g.strip()]
    if not gene_list:
        raise ValueError("Empty gene identifier input.")

    if all(is_uniprot_accession(g) for g in gene_list):
        uniprot_list = gene_list
    else:
        job_id = uniprot_initiate_id_mapping(gene_ids, "GeneCards", "UniProtKB")
        for _ in range(60):
            res = uniprot_check_id_mapping(job_id)
            if "results" in res:
                break
            time.sleep(1.0)
        else:
            raise TimeoutError("UniProt ID mapping did not return results in time.")
        uniprot_list = [r["to"] for r in res.get("results", [])]
        if not uniprot_list:
            raise RuntimeError("No UniProt accessions returned from mapping.")

    acc = uniprot_list[0]
    ebi_json = ebi_proteins_coordinates(acc)
    last_exon = pick_last_exon(ebi_json)
    seqs = compute_windows(last_exon)
    seqs["gene_id"] = gene_ids            # what the user entered (e.g., "TUBB" or "P07437")
    seqs["uniprot_id"] = acc       


    out = {
        "input_gene": gene_ids,
        "uniprot_accessions": uniprot_list,
        "last_exon_info": last_exon,
        "sequences": seqs,
    }

    (RESULTS_DIR / "step1_last_exon.json").write_text(json.dumps(last_exon, indent=2))
    (RESULTS_DIR / "step1_sequences.json").write_text(json.dumps(seqs, indent=2))

    (RESULTS_DIR / "step1_exon.fa").write_text(f">{gene_ids}|EXON_last\n{seqs['seq_exon']}\n")
    (RESULTS_DIR / "step1_crispor_input.fa").write_text(f">{gene_ids}|CRISPOR_input_50u_35d\n{seqs['seq_crispor_input']}\n")
    (RESULTS_DIR / "step1_guide_window.fa").write_text(f">{gene_ids}|GUIDE_window_501u_350d\n{seqs['seq_guide_window']}\n")
    (RESULTS_DIR / "step1_primer_template.fa").write_text(f">{gene_ids}|PRIMER_template_pm1000\n{seqs['seq_primer_template']}\n")

    return out

if __name__ == "__main__":
    if len(sys.argv) > 1:
        arg = sys.argv[1]
    else:
        arg = input("Enter Gene symbol or UniProt accession(s) (comma-separated): ")
    try:
        data = run(arg)
        print("\n[Step 1] Completed. Outputs saved in ./results")
        for f in ["step1_last_exon.json", "step1_sequences.json", "step1_exon.fa", "step1_crispor_input.fa", "step1_guide_window.fa", "step1_primer_template.fa"]:
            print(f"- {f}")
        le = data["last_exon_info"]
        print(f"Last exon → chr{le['chromosome']}:{le['start']}-{le['end']}  reverseStrand={le['reverseStrand']}")
    except Exception as e:
        print(f"[Step 1] ERROR: {e}")
        sys.exit(1)