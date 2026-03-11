import requests

print("[LOG] Initializing Pan-Dengue Target Finder...")

# UniProt IDs for the Polyproteins of Dengue Types 1 through 4
dengue_strains = {
    "DENV-1": "P33478",
    "DENV-2": "P29990",
    "DENV-3": "P27915",
    "DENV-4": "P09866"
}

sequences = {}

# 1. Fetching all sequences from the REST API
for strain, uid in dengue_strains.items():
    print(f"[LOG] Fetching data for {strain} (ID: {uid})...")
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Skip the first line (header) and join the sequence lines
        seq_lines = response.text.strip().split('\n')[1:]
        full_seq = "".join(seq_lines)
        sequences[strain] = full_seq
    else:
        print(f"[ERROR] Could not fetch {strain}")

# 2. Extracting Domain III (Approximated indices for demonstration)
# Note for developers: In production, you would use an alignment algorithm 
# like Clustal Omega here because the exact indices shift slightly due to mutations.
domain_3_targets = {}
for strain, seq in sequences.items():
    # Roughly extracting the Domain III region for comparison
    start_idx = 575
    end_idx = 675
    domain_3_targets[strain] = seq[start_idx:end_idx]

print("[LOG] --- ALIGNING SEQUENCES ---")

# 3. Finding the conserved (immutable) amino acids
# We compare letter by letter across all 4 strains
den1_target = domain_3_targets["DENV-1"]
den2_target = domain_3_targets["DENV-2"]
den3_target = domain_3_targets["DENV-3"]
den4_target = domain_3_targets["DENV-4"]

conserved_pattern = ""

for i in range(len(den1_target)):
    # If the amino acid is the exact same in all 4 types
    if den1_target[i] == den2_target[i] == den3_target[i] == den4_target[i]:
        conserved_pattern += den1_target[i]  # Add the letter
    else:
        conserved_pattern += "-"             # Add a dash for mutations

print("[LOG] Multiple Sequence Alignment complete.")
print("[LOG] Resulting Target Map (Letters are safe targets, dashes are dangerous mutations):")
print(f"[LOG] DENV-1: {den1_target}")
print(f"[LOG] SAFEMAP:{conserved_pattern}")