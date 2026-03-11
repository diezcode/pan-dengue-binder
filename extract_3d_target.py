import os
import requests
from Bio.PDB import PDBParser, PDBIO, Select

print("[SYSTEM] Booting 3D Coordinate Extractor...")

pdb_id = "1oan"
pdb_filename = f"{pdb_id}.pdb"

# 1. Download the 3D CAD file of the virus from the Protein Data Bank
if not os.path.exists(pdb_filename):
    print(f"[NETWORK] Downloading 3D structure for {pdb_id.upper()}...")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(pdb_filename, "w") as f:
            f.write(response.text)
        print("[SUCCESS] 3D structure downloaded.")
    else:
        raise Exception("Failed to download PDB file.")

# 2. Parse the 3D structure
print("[LOG] Parsing 3D coordinates into memory...")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("DENGUE_E_PROTEIN", pdb_filename)

# 3. Define the Target Zone
# In a full pipeline, you would map the string from the previous script to these numbers.
# For now, we are targeting Domain III (roughly residues 295 to 395).
target_start = 295
target_end = 395

print(f"[LOG] Slicing target zone: Residues {target_start} to {target_end}...")

# 4. Create a "Filter" to only save the atoms in our target zone
class TargetSelector(Select):
    def accept_residue(self, residue):
        # The ID of a residue in Biopython is a tuple; the sequence number is at index 1
        res_id = residue.get_id()[1] 
        if target_start <= res_id <= target_end:
            return True
        return False

# 5. Save the extracted 3D coordinates to a new file
io = PDBIO()
io.set_structure(structure)
output_filename = "dengue_target_motif.pdb"

# We pass our filter (TargetSelector) so it only writes the atoms we want
io.save(output_filename, TargetSelector())

print("\n[SUCCESS] Target extraction complete.")
print(f"[RESULT] Saved strictly the target 3D coordinates to: {output_filename}")
print("[SYSTEM] This file is now ready to be ingested by RFdiffusion.")