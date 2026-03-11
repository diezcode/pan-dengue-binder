from Bio import SeqIO

# File from the previous step
fasta_file = "P33478_dengue_type1.fasta"

print("[LOG] Starting sequence analysis for Dengue Virus...")

# SeqIO parses the FASTA file. 
# Although FASTA files can contain multiple sequences, ours has just one.
for record in SeqIO.parse(fasta_file, "fasta"):
    print(f"[LOG] Successfully loaded record ID: {record.id}")
    print(f"[LOG] Total polyprotein string length: {len(record.seq)} amino acids")

    # The entire string is the polyprotein. 
    # The Envelope (E) protein roughly spans from amino acid 281 to 775.
    # Python uses 0-based indexing, so the 281st amino acid is index 280.
    start_index_e_protein = 280
    end_index_e_protein = 775
    
    # Slicing the string to isolate the E-Protein
    e_protein_seq = record.seq[start_index_e_protein:end_index_e_protein]
    print(f"[LOG] Extracted Envelope (E) protein. Length: {len(e_protein_seq)} amino acids")
    
    # --- FINDING THE DOCKING STATION ---
    # The E-Protein has three domains. Domain III (DIII) is the part that actually 
    # binds to the human cell receptors. It is located roughly in the last 
    # 100 amino acids of the E-Protein.
    
    domain3_start = 295  # Relative to the E-Protein string
    domain3_end = 395    # Relative to the E-Protein string
    
    # Slicing again to isolate the specific Receptor Binding Domain
    receptor_binding_domain = e_protein_seq[domain3_start:domain3_end]
    
    print("[LOG] --- TARGET ACQUIRED ---")
    print("[LOG] Receptor Binding Domain (Domain III) sequence extracted:")
    print(f"[LOG] {receptor_binding_domain}")
    print("[LOG] Note: This is the critical target zone for building a neutralizing binder.")