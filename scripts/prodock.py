'''
Protein Docking Automation Script (config-based)
Author: Jordan Harrison
'''

#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
import stat
import subprocess
from rdkit import Chem

def remove_stereochemistry(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)

print("Starting ZDOCK docking automation (config mode)...")

base_dir = Path.cwd()
zdock_dir = base_dir / "ZDOCK"
scripts_dir = base_dir / "scripts"
shell_dir = base_dir / "shell"
protein_complexes_dir = base_dir / "complexes"
protein_complexes_dir.mkdir(exist_ok=True)

required_files = ["zdock", "create_lig", "create.pl", "mark_sur", "uniCHARMM", 'linkinvent.prior']

# === Read config.txt ===
config_path = base_dir / "config.txt"
if not config_path.exists():
    print("ERROR: config.txt not found!")
    exit(1)

with open(config_path, "r") as f:
    lines = [line.strip() for line in f if line.strip()]

if len(lines) < 4:
    print("ERROR: config.txt must have 4 lines: e3_ligasepdb, poipdb, e3ligand_smiles, poiligand_smiles")
    exit(1)

ligase_pdb, poi_pdb, lig1_smiles, lig2_smiles = lines

# Remove stereochemistry
lig1_smiles = remove_stereochemistry(lig1_smiles)
lig2_smiles = remove_stereochemistry(lig2_smiles)

complex_name = f"{Path(poi_pdb).stem}"
complex_dir = protein_complexes_dir / complex_name
complex_dir.mkdir(parents=True, exist_ok=True)

# Copy required ZDOCK files
for filename in required_files:
    src = zdock_dir / filename
    dest = complex_dir / filename
    shutil.copy(src, dest)
    dest.chmod(dest.stat().st_mode | stat.S_IEXEC)

# Copy scripts
for script_name in ["lig_dist.py", "prodock.py", "link_it.py", "dock.py", "analysis.py", "md.py"]:
    shutil.copy(scripts_dir / script_name, complex_dir)

for script_name in ["driver.sh", "link_it.sh"]:
    shutil.copy(shell_dir / script_name, complex_dir)

# Copy receptor and ligand PDBs
shutil.copy(poi_pdb, complex_dir / "receptor.pdb")
print(f"Copied {poi_pdb} to {complex_dir / 'receptor.pdb'}")
shutil.copy(ligase_pdb, complex_dir / "ligand.pdb")
print(f"Copied {ligase_pdb} to {complex_dir / 'ligand.pdb'}")

# Add dummy SEQRES
(complex_dir / "SEQRES").write_text("DUMMYSEQRES\n")

# Write smiles.smi
smiles_path = complex_dir / "smiles.smi"
with open(smiles_path, "w") as f:
    f.write(f"{lig1_smiles}|{lig2_smiles}\n")

# Execute SLURM script
slurm_script_path = complex_dir / "prodock.sh"
os.chmod(slurm_script_path, 0o755)

print("Submitting protein docking job...")
subprocess.run(["sbatch", "prodock.sh"])
print("Job submitted.")