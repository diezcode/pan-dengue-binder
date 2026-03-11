import requests

print("[SYSTEM] Booting Pan-Dengue Target Finder...")

# The official UniProt IDs for the Dengue 1-4 Polyproteins
DENV_STRAINS = {
    "DENV-1": "P33478",
    "DENV-2": "P29990",
    "DENV-3": "P27915",
    "DENV-4": "P09866"
}

def fetch_sequence(strain_name, uniprot_id):
    """Fetches the raw amino acid string from the biological database."""
    print(f"[NETWORK] Fetching source code for {strain_name} ({uniprot_id})...")
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Strip the header and merge the multiline sequence into one massive string
        lines = response.text.strip().split('\n')[1:]
        return "".join(lines)
    else:
        raise Exception(f"Failed to download {strain_name}")

# 1. Download all sequences into memory
sequences = {}
for strain, uid in DENV_STRAINS.items():
    sequences[strain] = fetch_sequence(strain, uid)

print("\n[LOG] All sequences loaded successfully. Initiating alignment scan...")
print("[LOG] Scanning for immutable protein chains (minimum length: 10 amino acids)...\n")

def find_immutable_targets(seqs, min_length=10):
    """
    Scans DENV-1 and checks if that exact sequence exists in DENV 2, 3, and 4.
    This is an algorithmic search for biological conserved regions.
    """
    reference = seqs["DENV-1"]
    targets = []
    
    i = 0
    while i < len(reference) - min_length + 1:
        # Start with the minimum window size
        window = reference[i:i+min_length]
        
        # Check if this exact string is in all other Dengue strains
        if all(window in seqs[strain] for strain in ["DENV-2", "DENV-3", "DENV-4"]):
            
            # If it is, keep expanding the window to the right until it breaks
            current_target = window
            expand_idx = i + min_length
            
            while expand_idx < len(reference):
                test_window = reference[i:expand_idx + 1]
                if all(test_window in seqs[strain] for strain in ["DENV-2", "DENV-3", "DENV-4"]):
                    current_target = test_window
                    expand_idx += 1
                else:
                    break
            
            targets.append(current_target)
            # Jump ahead to avoid overlapping targets
            i = expand_idx 
        else:
            i += 1
            
    return targets

# 2. Run the algorithmic scan
safe_zones = find_immutable_targets(sequences, min_length=10)

print("[SUCCESS] Target Acquisition Complete.")
print(f"[RESULT] Found {len(safe_zones)} completely conserved target zones across all 4 strains.\n")

# 3. Output the targets for the AI generator
for idx, target in enumerate(safe_zones, 1):
    print(f"--- TARGET ZONE {idx} ---")
    print(f"Length: {len(target)} amino acids")
    print(f"Sequence: {target}")
    print("-" * 25)

print("\n[SYSTEM] Ready for Phase 2: Feeding 3D coordinates to RFdiffusion.")