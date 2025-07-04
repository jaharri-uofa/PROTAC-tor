#!/usr/bin/env python3
"""
ZDOCK docking automation script for a pair of protein complexes.
Author: Jordan Harrison

Usage:
    python3 run_docking.py path/to/complex1.pdb path/to/complex2.pdb
"""

import os
import shutil
from pathlib import Path
import stat
import subprocess
import argparse

# ----------------------------
# Parse command-line arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Run ZDOCK docking job on two input PDBs.")
parser.add_argument("complex1_pdb", help="Path to the first protein complex PDB file")
parser.add_argument("complex2_pdb", help="Path to the second protein complex PDB file")
args = parser.parse_args()

complex1_pdb = Path(args.complex1_pdb).resolve()
complex2_pdb = Path(args.complex2_pdb).resolve()

if not complex1_pdb.exists():
    raise FileNotFoundError(f"Complex 1 PDB not found: {complex1_pdb}")
if not complex2_pdb.exists():
    raise FileNotFoundError(f"Complex 2 PDB not found: {complex2_pdb}")

print(f"Setting up docking job for: {complex1_pdb.name} × {complex2_pdb.name}")

# ----------------------------
# Set up project directories
# ----------------------------
base_dir = Path.cwd()
protein_complexes_dir = base_dir / "Protein_Complexes"
complex_name = f"{complex1_pdb.stem}_{complex2_pdb.stem}"
complex_dir = protein_complexes_dir / complex_name
complex_dir.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Copy required binaries/scripts
# ----------------------------
required_files = ["zdock", "create_lig", "create.pl", "mark_sur", "uniCHARMM"]
for filename in required_files:
    src = base_dir / filename
    dest = complex_dir / filename
    shutil.copy(src, dest)
    dest.chmod(dest.stat().st_mode | stat.S_IEXEC)

# ----------------------------
# Copy PDB files
# ----------------------------
shutil.copy(complex1_pdb, complex_dir / "receptor.pdb")
shutil.copy(complex2_pdb, complex_dir / "ligand.pdb")

# ----------------------------
# Add dummy SEQRES file
# ----------------------------
(complex_dir / "SEQRES").write_text("DUMMYSEQRES\n") # Does this do anything?

# ----------------------------
# Write SLURM job script
# ----------------------------
slurm_script_path = complex_dir / "run_docking.sh"
with open(slurm_script_path, "w") as f:
    f.write(f"""#!/bin/bash
#SBATCH --job-name={complex_name}
#SBATCH --output={complex_dir}/zdock.out
#SBATCH --error={complex_dir}/zdock.err
#SBATCH --time=0-01:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

export LD_LIBRARY_PATH=/home/jordanha/zdock_libs/usr/lib64:$LD_LIBRARY_PATH

cp lig_dist.py "{complex_dir}"
cp link.py "{complex_dir}"
cd "{complex_dir}"

# Preprocess receptor and ligand
./mark_sur receptor.pdb {complex1_pdb.stem}.pdb
./mark_sur ligand.pdb {complex2_pdb.stem}.pdb

# Run ZDOCK
./zdock -R {complex1_pdb.stem}.pdb -L {complex2_pdb.stem}.pdb -o zdock_result.out
./create.pl zdock_result.out

# Calculate ligand distance
module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
python lig_dist.py

# Run LinkInvent
python run_linkinvent.py --config linkinvent_config.json
""")
os.chmod(slurm_script_path, 0o755)

# ----------------------------
# Submit SLURM job
# ----------------------------
print(f"Submitting: {slurm_script_path}")
subprocess.run(["sbatch", str(slurm_script_path)])

print("ZDOCK job submitted successfully.")
