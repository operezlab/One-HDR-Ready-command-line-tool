import requests, sys, subprocess, csv,time
import pandas as pd
from Bio.Seq import Seq

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

# Get the Gene IDs from user input
gene_ids = input("Enter the Gene IDs (comma-separated): ")

# Initiate ID mapping from GeneCards to UniProtKB
job_id = initiate_id_mapping(gene_ids, "GeneCards", "UniProtKB")
print(f"Job ID: {job_id}")

# Poll the UniProt API for the results
while True:
    results = check_id_mapping_results(job_id)
    if "results" in results:
        break
    print("Waiting for results...")
    time.sleep(5)

# Extract UniProt accession codes
uniprot_accession_codes = [result["to"] for result in results["results"]]

print("UniProt Accession Codes:", uniprot_accession_codes)

# For simplicity, we'll use the first UniProt accession code
if uniprot_accession_codes:
    uniprot_accession_code = uniprot_accession_codes[0]
else:
    print("No UniProt accession codes found.")
    sys.exit()

# Function to get the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Define the URL for fetching exon coordinates
requestURL = f"https://www.ebi.ac.uk/proteins/api/coordinates/{uniprot_accession_code}"

r = requests.get(requestURL, headers={"Accept": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

response_data = r.json()

# Extract relevant information from each gene
relevant_info = []

for gene in response_data.get("gnCoordinate", []):
    ensembl_gene_id = gene.get("ensemblGeneId")
    genomic_location = gene.get("genomicLocation", {}) 

    if "exon" in genomic_location:
        for exon in genomic_location["exon"]:
            exon_id = exon.get("id")
            chromosome = genomic_location.get("chromosome")
            start = exon.get("genomeLocation", {}).get("begin", {}).get("position")
            end = exon.get("genomeLocation", {}).get("end", {}).get("position")
            relevant_info.append({
                "ensembl_gene_id": ensembl_gene_id,
                "exon_id": exon_id,
                "chromosome": chromosome,
                "start": start,
                "end": end
            })

# Filter relevant_info based on the extracted Ensembl gene ID
if relevant_info:
    ensembl_gene_id = relevant_info[0]['ensembl_gene_id']
    filtered_relevant_info = [info for info in relevant_info if info['ensembl_gene_id'] == ensembl_gene_id]

    # Print the filtered relevant information for debugging
    print("Filtered relevant information:")
    for info in filtered_relevant_info:
        print(info)

    # Label each exon with a sequential identifier
    for idx, info in enumerate(filtered_relevant_info, start=1):
        info['label'] = idx

    # Find the exon with the highest label
    last_exon_info = max(filtered_relevant_info, key=lambda x: x['label'])

    # Print the last exon information for debugging
    print("Last exon information:")
    print(last_exon_info)
else:
    print("No relevant exon information found.")

# Fetch the DNA sequence for the last exon
if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    start = exon_end - 500
    end = exon_end + 350
    
    if exon_end < exon_start:
        start = exon_end - 350
        end = exon_end + 500
        exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    print(f"New coordinates: Chromosome {chromosome}, Start {start}, End {end}")

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
        print("Reverse complement of the extended DNA sequence:")
        gdna_sequence = reverse_DNA_extended
    else:
        print(f"DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        gdna_sequence = dna_sequence

    print(gdna_sequence)

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
    print(f"Amino acid sequence for the last exon ({last_exon_info['exon_id']}):")
    print(amino_acid_seq)

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
    print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

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
        print("Reverse complement of the extended DNA sequence +/- 800 bp:")
        exon_seq_800 = reverse_complement_extended
    else:
        print(f"DNA sequence for the region surrounding the last exon +/- 800 bp({last_exon_info['exon_id']}):")
        exon_seq_800 = whole_dna_sequence_800

print(exon_seq_800)

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
    print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

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


with open("sequence.fa", "w") as fa_file:
            fa_file.write(f">{uniprot_accession_code}_whole\n")
            fa_file.write(CRISPRtgSearch + "\n")

# Define the command and arguments
discover_command = [
    'java', '-Xmx4g', '-jar', 'FlashFry-assembly-1.15.jar', 
    'discover', 
    '--database', 'GRCh38_cas9ngg_database', 
    '--fasta', 'sequence.fa', 
    '--output', 'CRISPRtg.output'
]

# Define the second command and arguments
score_command = [
    'java', '-Xmx4g', '-jar', 'FlashFry-assembly-1.15.jar', 
    'score', 
    '--input', 'CRISPRtg.output', 
    '--output', 'CRISPRtg.output.scored.tsv', 
    '--scoringMetrics', 'doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,moreno2015', 
    '--database', 'GRCh38_cas9ngg_database'
]

# Run the first command
try:
    print("Running discover command...")
    discover_result = subprocess.run(discover_command, capture_output=True, text=True, check=True)
    print("Discover command executed successfully")
    print("Discover Command Output:\n", discover_result.stdout)
    print("Discover Command Errors:\n", discover_result.stderr)
except subprocess.CalledProcessError as e:
    print(f"Discover command failed with return code {e.returncode}")
    print(f"Discover Command Output:\n{e.output}")
    print(f"Discover Command Errors:\n{e.stderr}")
    exit(1)

# Run the second command
try:
    print("Running score command...")
    score_result = subprocess.run(score_command, capture_output=True, text=True, check=True)
    print("Score command executed successfully")
    print("Score Command Output:\n", score_result.stdout)
    print("Score Command Errors:\n", score_result.stderr)
except subprocess.CalledProcessError as e:
    print(f"Score command failed with return code {e.returncode}")
    print(f"Score Command Output:\n{e.output}")
    print(f"Score Command Errors:\n{e.stderr}")
    exit(1)



# Define a function to calculate distance from exon for a target sequence
def calculate_distance(row):
    target_seq = row['target']
    orientation = row['orientation']
    print("target seq")
    print(target_seq)

    if orientation == "FWD":
        # Take 5 bp from the end of the target sequence
        segment = target_seq
    else:
        # Reverse complement the target sequence if orientation is RVS
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

    print("SegmentSearch")
    print(segment)
    print("This is where CRISPR is targeting")
    print(target_position)
    print("last_exon_position")
    print(last_exon_position)
    print("last letter position")
    print(last_letter_position)
    # Calculate the distance from the last letter of the exon sequence
    if orientation == "FWD":
        distance = target_position - last_letter_position + 16
    else:
        distance = target_position - last_letter_position + 5
        #if distance < 0:
            #distance = distance - 1
    print(distance)
    return distance

if last_exon_info['end'] < last_exon_info['start']:
        rvs_seq = reverse_complement_exon[-15:]
        print(rvs_seq)
        last_exon_position = reverse_complement_CRISPRi.rfind(rvs_seq)
        last_letter_position = last_exon_position + len(rvs_seq) - 1
else:
        fwd_seq = exon_dna_sequence[-15:]
        print(fwd_seq)
        last_exon_position = whole_dna_sequence.rfind(fwd_seq)
        last_letter_position = last_exon_position + len(fwd_seq) - 1


# Print the positions for debugging
print("Last exon sequence position in whole DNA sequence:", last_exon_position)
print("Last letter position of exon sequence:", last_letter_position)
print("CRISPR Search Input Seq", whole_dna_sequence)


try:
    df = pd.read_csv('CRISPRtg.output.scored.tsv', sep='\t')
    print("DataFrame loaded successfully:")
    print(df.head())
    
    df['Distance from Exon'] = df.apply(calculate_distance, axis=1)

    # Sort the DataFrame based on the specified criteria and keep only the top 5 entries
    df['Distance from Exon'] = df.apply(calculate_distance, axis=1)

    # Sort the DataFrame based on the specified criteria and keep only the top entries based on scores and distance
    filtered_df = df.sort_values(
        by=[
            'Hsu2013',  # Higher is better
            'DoenchCFD_maxOT',  # Lower is better
            'DoenchCFD_specificityscore',  # Higher is better
            'Moreno-Mateos2015OnTarget',  # Higher is better
            'Doench2014OnTarget',  # Higher is better
            'Distance from Exon'  # Closer to zero is better
        ],
        ascending=[False, True, False, False, False, True]  # Specify the sort order for each column
    ).head(20)
    
    filtered_df = filtered_df.iloc[filtered_df['Distance from Exon'].abs().argsort()[:20]]

    # Select the specified columns
    selected_columns = [
        'target', 
        'orientation', 
        'Doench2014OnTarget', 
        'DoenchCFD_maxOT', 
        'DoenchCFD_specificityscore', 
        'Hsu2013',
        'Moreno-Mateos2015OnTarget',
        'Distance from Exon'
        
    ]
    filtered_df = filtered_df[selected_columns]
    
    # Print the DataFrame with the new column
    print(filtered_df[selected_columns])

    # Save the filtered DataFrame to a CSV file
    filtered_df.to_csv(f'{uniprot_accession_code}_CRISPR_tgts.csv', index=False)
    print(f"Filtered DataFrame saved to filtered_crispr_data.csv {uniprot_accession_code}_CRISPR_tgts.csv")

    # Optionally, save the filtered DataFrame to an Excel file
    #filtered_df.to_excel('filtered_crispr_data.xlsx', index=False)
    #print(f"{uniprot_accession_code}_CRISPR_tgts.xlsx")

    new_selected_columns = ['target', 'orientation', 'Distance from Exon']

    while True:
        try:
            selected_row = int(input("Enter the row number: "))
            if selected_row in filtered_df.index:
                selected_target = filtered_df.loc[selected_row, new_selected_columns]
                print("\nYou have selected the following target:")
                print(selected_target)
                break
            else:
                print("Invalid row number. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a valid row number.")

    selected_target_seq = selected_target['target']
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
    print(f"An error occurred while loading the DataFrame: {e}")



def find_and_extract_sequence(gdna, target, upstream, downstream):
    target_index = gdna.find(target)
    if target_index == -1:
        return None, None
    upstream_sequence = gdna[max(0, target_index + len(target) - upstream):target_index + len(target)]
    downstream_sequence = gdna[target_index + len(target) - 1:target_index + len(target) - 1 + downstream]
    return upstream_sequence, downstream_sequence

gene_name = gene_ids
sgRNA_sequence = selected_target_seq

if selected_orientation == "RVS":
    is_sgRNA_inverted = "Y"
else:
    is_sgRNA_inverted = "N"

original_sgRNA = sgRNA_sequence

if is_sgRNA_inverted == "Y":
    sgRNA_sequence = reverse_complement(original_sgRNA)

if len(last_exon_seq) >= 20:
    protein_coding_exon = last_exon_seq[-20:]
else:
    protein_coding_exon = last_exon_seq

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
print(f"Right_arm sequence:\n{right_arm}")