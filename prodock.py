#!/usr/bin/env python3
'''
ZDOCK docking automation script.

Need to make this work in the form of calling the script with the directory containing PDB files as an argument.
Also consider what libraries are needed to run this script, and how to handle the output files.
'''

import os
import shutil
from pathlib import Path
import stat
import subprocess

print("Starting ZDOCK docking automation...")

# Add in logic to input pdb file names
base_dir = Path.cwd()
protein_complexes_dir = base_dir / "Protein_Complexes"
ad_kinases_dir = protein_complexes_dir / "AD_Kinases"
e3_dir = protein_complexes_dir / "E3_Gallus_gallus"

required_files = ["zdock", "create_lig", "create.pl", "mark_sur", "uniCHARMM"]

# Loop over AD Kinases and E3
for ad_pdb in ad_kinases_dir.glob("*.pdb"):
    for e3_pdb in e3_dir.glob("*.pdb"):
        # Create complex directory
        complex_name = f"{ad_pdb.stem}_{e3_pdb.stem}"
        complex_dir = protein_complexes_dir / complex_name
        complex_dir.mkdir(parents=True, exist_ok=True)

        # Copy required binaries/scripts into complex dir
        for filename in required_files:
            src = base_dir / filename
            dest = complex_dir / filename
            shutil.copy(src, dest)
            dest.chmod(dest.stat().st_mode | stat.S_IEXEC)

        # Copy ligand and receptor PDBs into directory with expected names
        # Can be used for ligand extraction
        shutil.copy(ad_pdb, complex_dir / "receptor.pdb")
        shutil.copy(e3_pdb, complex_dir / "ligand.pdb")

        # Add dummy SEQRES file
        (complex_dir / "SEQRES").write_text("DUMMYSEQRES\n")

        # Write SLURM script
        slurm_script_path = complex_dir / "run_docking.sh"
        with open(slurm_script_path, "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={complex_name}
#SBATCH --output={complex_dir}/zdock.out
#SBATCH --error={complex_dir}/zdock.err
#SBATCH --time=0-01:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

export LD_LIBRARY_PATH=/home/jordanha/zdock_libs/usr/lib64:$LD_LIBRARY_PATH

cp lig_dist.py "{complex_dir}"
cp link.py "{complex_dir}"
cd "{complex_dir}"

# Preprocess receptor and ligand
./mark_sur receptor.pdb kinase.pdb
./mark_sur ligand.pdb e3.pdb

# Run ZDOCK
./zdock -R kinase.pdb -L e3.pdb -o zdock_result.out
./create.pl zdock_result.out

# Calculate ligand distance
module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
python lig_dist.py

# LinkInvent
python run_linkinvent.py --config linkinvent_config.json
""")
        os.chmod(slurm_script_path, 0o755)

print("Submitting all docking jobs...")
for sh_file in Path("Protein_Complexes").rglob("run_docking.sh"):
    print(f"Submitting: {sh_file}")
    subprocess.run(["sbatch", str(sh_file)])

print("All docking jobs prepared. Check the Protein_Complexes directory for details.")