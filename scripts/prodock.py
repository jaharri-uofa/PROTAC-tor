#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
import stat
import subprocess

print("Starting ZDOCK docking automation...")

base_dir = Path.cwd()
zdock_dir = base_dir / "ZDOCK"
scripts_dir = base_dir / "scripts"
protein_complexes_dir = base_dir / "Protein_complexes"

# Get all .pdb files in the root of PROTAC-tor
pdb_files = list(base_dir.glob("*.pdb"))

# Make all pairwise combinations of pdb files
for i, pdb1 in enumerate(pdb_files):
    for j, pdb2 in enumerate(pdb_files):
        if i >= j:
            continue  # Avoid duplicates and self-pairing

        complex_name = f"{pdb1.stem}_{pdb2.stem}"
        complex_dir = protein_complexes_dir / complex_name
        complex_dir.mkdir(parents=True, exist_ok=True)

        # Copy required ZDOCK files from ZDOCK dir
        required_files = ["zdock", "create_lig", "create.pl", "mark_sur", "uniCHARMM"]
        for filename in required_files:
            src = zdock_dir / filename
            dest = complex_dir / filename
            shutil.copy(src, dest)
            dest.chmod(dest.stat().st_mode | stat.S_IEXEC)

        # Copy scripts
        for script_name in ["Lig_dist.py", "prodock.py", "link_it.py"]:
            shutil.copy(scripts_dir / script_name, complex_dir)

        # Copy receptor and ligand PDBs
        shutil.copy(pdb1, complex_dir / "receptor.pdb")
        shutil.copy(pdb2, complex_dir / "ligand.pdb")

        # Add dummy SEQRES
        (complex_dir / "SEQRES").write_text("DUMMYSEQRES\n")

        # Write SLURM script
        slurm_script_path = complex_dir / "run_docking.sh"
        with open(slurm_script_path, "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={complex_name}
#SBATCH --output=zdock.out
#SBATCH --error=zdock.err
#SBATCH --time=0-01:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

export LD_LIBRARY_PATH=/home/jordanha/zdock_libs/usr/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/REINVENT4:$PYTHONPATH

module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1


cd "{complex_dir}"

# Preprocess
./mark_sur receptor.pdb Receptor.pdb
./mark_sur ligand.pdb Ligand.pdb

# Run ZDOCK
./zdock -R Ligand.pdb -L Receptor.pdb -o zdock_result.out
./create.pl zdock_result.out

# Calculate distances
python Lig_dist.py

# LinkInvent
echo "Running LinkInvent..."
~/reinvent4/bin/python link_it.py --smiles_csv smiles.csv --dist_file input.txt --output_json linkinvent_config.json --slurm_script submit_linkinvent.sh
""")
        os.chmod(slurm_script_path, 0o755)

print("Submitting all docking jobs...")
for sh_file in protein_complexes_dir.rglob("run_docking.sh"):
    subprocess.run(["sbatch", str(sh_file)])

print("All jobs submitted.")
