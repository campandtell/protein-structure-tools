from pymol import cmd
import os
import glob

def load_and_align():
    # Get all PDB files in current directory
    pdb_files = glob.glob("*.pdb")
    if not pdb_files:
        print("No PDB files found in current directory")
        return
        
    # Load all PDB files
    for pdb in pdb_files:
        name = os.path.splitext(pdb)[0]
        cmd.load(pdb, name)
    
    # Get the first structure as reference
    reference = os.path.splitext(pdb_files[0])[0]
    print(f"Using {reference} as reference structure")
    
    # Align all other structures to the reference
    for pdb in pdb_files[1:]:
        mobile = os.path.splitext(pdb)[0]
        cmd.align(mobile, reference)
        print(f"Aligned {mobile} to {reference}")
    
    # Center the view
    cmd.center()
    cmd.zoom()

# Run the function
load_and_align()
