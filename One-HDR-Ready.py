import requests, sys, subprocess, csv,time, os
import pandas as pd
from Bio.Seq import Seq
from Bio.Data import CodonTable


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
    time.sleep(5)

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

        for gene in response_data.get("gnCoordinate", []):
            ensembl_gene_id = gene.get("ensemblGeneId")
            genomic_location = gene.get("genomicLocation", {})
            #print("Genomic Location:", genomic_location)  # Debugging

            if "exon" in genomic_location:
                #print("Exons Found:", genomic_location["exon"])  # Debugging

                for exon in genomic_location["exon"]:
                    #print("Processing Exon:", exon)  # Debugging
                    exon_id = exon.get("id")
                    chromosome = genomic_location.get("chromosome")
                    start = exon.get("genomeLocation", {}).get("begin", {}).get("position")
                    end = exon.get("genomeLocation", {}).get("end", {}).get("position")

                    # Debugging Start and End Positions
                    #print(f"Exon ID: {exon_id}, Start: {start}, End: {end}")

                    # Ensuring valid start and end positions
                    if start is not None and end is not None:
                        relevant_info.append({
                            "ensembl_gene_id": ensembl_gene_id,
                            "exon_id": exon_id,
                            "chromosome": chromosome,
                            "start": start,
                            "end": end
                        })

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
        #print("Reverse complement of the Exon sequence:")
        last_exon_seq = reverse_complement_exon
else:
        #print(f"DNA sequence for the last exon({last_exon_info['exon_id']}):")
        last_exon_seq = exon_dna_sequence

#print(last_exon_seq)

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
        #print(f"Reverse complement of Whole DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        CRISPRtgSearch = reverse_complement_CRISPRi
    else:
        #print(f"Whole DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        CRISPRtgSearch = whole_dna_sequence  

#print(CRISPRtgSearch)

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

# Write sequence to file in the current directory
with open("sequence.fa", "w") as fa_file:
    fa_file.write(f">{gene_ids}_whole\n")
    fa_file.write(CRISPRtgSearch + "\n")

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

    # Ensure numeric values for score columns
    score_columns = ['mitSpecScore', 'cfdSpecScore', 'Moreno-Mateos-Score', "Doench '16-Score"]
    for col in score_columns:
        filtered_df[col] = pd.to_numeric(filtered_df[col], errors='coerce')  # Convert to numeric, set errors to NaN

    # Drop rows where any score is missing (NaN)
    filtered_df = filtered_df.dropna(subset=score_columns)

    # Normalize scores for fair weighting
    filtered_df['mitSpecScore_norm'] = (filtered_df['mitSpecScore'] - filtered_df['mitSpecScore'].min()) / (filtered_df['mitSpecScore'].max() - filtered_df['mitSpecScore'].min())
    filtered_df['cfdSpecScore_norm'] = (filtered_df['cfdSpecScore'] - filtered_df['cfdSpecScore'].min()) / (filtered_df['cfdSpecScore'].max() - filtered_df['cfdSpecScore'].min())
    filtered_df['Moreno-Mateos-Score_norm'] = (filtered_df['Moreno-Mateos-Score'] - filtered_df['Moreno-Mateos-Score'].min()) / (filtered_df['Moreno-Mateos-Score'].max() - filtered_df['Moreno-Mateos-Score'].min())
    filtered_df["Doench '16-Score_norm"] = (filtered_df["Doench '16-Score"] - filtered_df["Doench '16-Score"].min()) / (filtered_df["Doench '16-Score"].max() - filtered_df["Doench '16-Score"].min())

    # Define weights: on-target scores have higher weight
    weight_mitSpecScore = 0.2
    weight_cfdSpecScore = 0.2
    weight_MorenoMateos = 0.3  # Higher weight for on-target score
    weight_Doench16 = 0.3       # Higher weight for on-target score

    # Calculate the weighted combined ranking score
    filtered_df['combined_score'] = (
        (filtered_df['mitSpecScore_norm'] * weight_mitSpecScore) +
        (filtered_df['cfdSpecScore_norm'] * weight_cfdSpecScore) +
        (filtered_df['Moreno-Mateos-Score_norm'] * weight_MorenoMateos) +
        (filtered_df["Doench '16-Score_norm"] * weight_Doench16)
    )

    # Check if the DataFrame is empty after filtering
    if filtered_df.empty:
        print("Error: No valid rows found after processing. Check if all values in score columns are numeric.")
    else:
        # Select the row with the highest combined score
        selected_target = filtered_df.loc[filtered_df['combined_score'].idxmax(), ['targetSeq', 'orientation', 'Distance from Exon']]

        # Extract selected values
        selected_target_seq = selected_target['targetSeq']
        selected_orientation = selected_target['orientation']
        selected_distance_from_exon = selected_target['Distance from Exon']

        # Print the selected values
        print(f"Selected target sequence: {selected_target_seq}")
        print(f"Selected orientation: {selected_orientation}")
        print(f"Selected distance from exon: {selected_distance_from_exon}")


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
print(left_arm)

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
print(sgrna_in_left_arm_start, sgrna_in_left_arm_end)


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

    if is_sgRNA_inverted == "Y":
        sgRNA_sequence = sgRNA_core + reverse_complement(original_sgRNA[:3])  # Reattach PAM at the end
    else:
        sgRNA_sequence = sgRNA_core + original_sgRNA[-3:]

    left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)

    if is_sgRNA_inverted == "Y":
        _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
    else:
        _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

    left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
    if is_sgRNA_inverted == "Y":
        left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", left_arm)
    else:
        left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", left_arm)
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
    print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
    print(f"Left_arm sequence:\n{left_arm}")
    pass

else:
    
    print("Handling negative distance from exon...")

    
    if is_sgRNA_inverted == "Y":
        print(f"main sgRNA_Sequence {original_sgRNA}")
        print(f"no PAM sgRNA_Sequence {sgRNA_core}")
        sgRNA_sequence = reverse_complement(original_sgRNA[-3:]) + sgRNA_core   # Reattach PAM at the end
        
    else:
        print(f"main sgRNA_Sequence {original_sgRNA}")
        print(f"no PAM sgRNA_Sequence {sgRNA_core}")
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
        #print(f"sgRNA found at positions {sgrna_start} to {sgrna_end}")
        sgrna_sequence_found = dna_sequence_trimmed[sgrna_start:sgrna_end]
        
        #print(f"Base pairs in the sgRNA positions: {sgrna_sequence_found}")
        if is_sgRNA_inverted == "Y":
            position_1 = codon_order_sequence[sgrna_start]
            position_2 = codon_order_sequence[sgrna_start + 1]

            if position_2 == '3' or position_1 == '3':
            
                #print("Pathway: One or both of the last two bases are in codon position 3.")
                left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)

                # Original pathway for non-negative distance from exon
                left_arm_start, left_arm_end = find_sgrna_coordinates_in_exon(gdna_sequence, left_arm)
                #print(f"Left arm coordinates in exon: start={left_arm_start}, end={left_arm_end}")

                sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_partial_sgrna_coordinates_in_exon(left_arm, sgRNA_core, min_match_length)
                print(sgrna_in_left_arm_start, sgrna_in_left_arm_end)

                if position_2 == '3':                        
                        print("Attempting to Edit PAM")
                        print(f"Current left arm sequence: {left_arm}")
                        #rev_left_arm = reverse_complement(left_arm)
                        #print(f"rev comp of left arm {rev_left_arm}")
                        #print(f"{sgrna_in_left_arm_start},{sgrna_in_left_arm_end}")

                            # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                        if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                            # Get the last three bases of the sgRNA sequence in the left arm
                            last_three_start = sgrna_in_left_arm_start-4
                            last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_start-1]
                            print(f"Last codon of sgRNA in left arm: {last_three_bases}")

                            modified_left_arm, verification_status = modify_sgrna_PAM(
                                left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                                sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                            )

                            if is_sgRNA_inverted == "Y":
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
                            else:
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

                            left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
                            if is_sgRNA_inverted == "Y":
                                left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", modified_left_arm)
                            else:
                                left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", modified_left_arm)
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
                            print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
                            print(f"Left_arm sequence:\n{modified_left_arm}")

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
                        #print(f"Last three bases of sgRNA in left arm: {last_three_bases}")

                        modified_left_arm, verification_status = modify_sgrna_PAM(
                        left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                        sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                        )

                        if is_sgRNA_inverted == "Y":
                            _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
                        else:
                            _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

                        left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
                        if is_sgRNA_inverted == "Y":
                            left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", modified_left_arm)
                        else:
                            left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", modified_left_arm)
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
                        print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
                        print(f"Left_arm sequence:\n{modified_left_arm}")
                

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


                if is_sgRNA_inverted == "Y":
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
                else:
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

                left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
                if is_sgRNA_inverted == "Y":
                    left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", modified_left_arm)
                else:
                    left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", modified_left_arm)
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
                print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
                print(f"Left_arm sequence:\n{modified_left_arm}")

            
        else:    
            # Pathway based on codon positions
            position_2 = codon_order_sequence[sgrna_end - 2]
            position_3 = codon_order_sequence[sgrna_end - 1]
            if position_2 == '3' or position_3 == '3':

                #print("Pathway: One or both of the last two bases are in codon position 3.")
                left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)

                # Original pathway for non-negative distance from exon
                left_arm_start, left_arm_end = find_sgrna_coordinates_in_exon(gdna_sequence, left_arm)
                #print(f"Left arm coordinates in exon: start={left_arm_start}, end={left_arm_end}")

                # Find sgRNA amino acids within the left arm sequence
                sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
                #print(sgrna_in_left_arm_start,sgrna_in_left_arm_end)
                
                # Extract amino acid sequence within the left arm coordinates
                if left_arm_start is not None and left_arm_end is not None:
                    #left_arm_amino_acids = translate_dna_to_protein(dna_sequence_trimmed[left_arm_start:left_arm_end])
                    #print(f"Amino acids for left arm: {left_arm_amino_acids}")
                    if position_2 == '3':

                        # Find sgRNA amino acids within the left arm sequence
                        sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
                        #print(sgrna_in_left_arm_start,sgrna_in_left_arm_end)

                        #print(f"Current left arm sequence: {left_arm}")

                            # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                        if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                            # Get the last three bases of the sgRNA sequence in the left arm
                            last_three_start = sgrna_in_left_arm_end - 4
                            last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_end-1]
                            #print(f"Last codon of sgRNA in left arm: {last_three_bases}")
                            modified_left_arm, verification_status = modify_sgrna_PAM(
                                left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                                sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                            )

                            if is_sgRNA_inverted == "Y":
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
                            else:
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

                            left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
                            if is_sgRNA_inverted == "Y":
                                left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", modified_left_arm)
                            else:
                                left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", modified_left_arm)
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
                            print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
                            print(f"Left_arm sequence:\n{modified_left_arm}")

                        else:
                            print("sgRNA sequence not found in the left arm.")
                        
                    else:
                        # Find sgRNA amino acids within the left arm sequence
                        sgrna_in_left_arm_start, sgrna_in_left_arm_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
                        #print(sgrna_in_left_arm_start,sgrna_in_left_arm_end)

                        #print(f"Current left arm sequence: {left_arm}")

                            # To extract the last three coordinates from the sgRNA sequence in the left arm sequence
                        if sgrna_in_left_arm_start is not None and sgrna_in_left_arm_end is not None:
                            # Get the last three bases of the sgRNA sequence in the left arm
                            last_three_start = sgrna_in_left_arm_end - 3
                            last_three_bases = left_arm[last_three_start:sgrna_in_left_arm_end]
                            #print(f"Last three bases of sgRNA in left arm: {last_three_bases}")

                            modified_left_arm, verification_status = modify_sgrna_PAM(
                            left_arm, last_three_start, last_three_bases, dna_sequence_trimmed, 
                            sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                            )

                            if is_sgRNA_inverted == "Y":
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
                            else:
                                _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

                            left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
                            if is_sgRNA_inverted == "Y":
                                left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", modified_left_arm)
                            else:
                                left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", modified_left_arm)
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
                            print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
                            print(f"Left_arm sequence:\n{modified_left_arm}")
                else:
                    print("sgRNA sequence not found in the left arm.")
        

    else:
        print("Pathway: Neither of the last two bases is in codon position 3.")

        modified_left_arm, verification_status = modify_sgrna_codons(
            dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
        )

        if is_sgRNA_inverted == "Y":
            _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
        else:
            _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

        left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
        if is_sgRNA_inverted == "Y":
            left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", modified_left_arm)
        else:
            left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", modified_left_arm)
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
        print(f">{gene_name}-Right-AI1\n{right_ai1_seq}")
        print(f"Left_arm sequence:\n{modified_left_arm}")
 

# Print the results

print(f"Right_arm sequence:\n{right_arm}")
