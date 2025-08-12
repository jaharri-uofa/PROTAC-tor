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
for script_name in ["Lig_dist.py", "prodock.py", "link_it.py", "dock.py", "analysis.py", "md.py"]:
    shutil.copy(scripts_dir / script_name, complex_dir)

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

# Write SLURM script
slurm_script_path = complex_dir / "run_docking.sh"
with open(slurm_script_path, "w") as f:
    f.write(f"""#!/bin/bash
#SBATCH --job-name={complex_name}
#SBATCH --output=zdock.out
#SBATCH --error=zdock.err
#SBATCH --time=0-01:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.6
module load python-build-bundle/2025b

if [ ! -d "xxhash" ]; then
    pip install xxhash
fi

if [ ! -f "$HOME/dssp/build/mkdssp" ]; then
    cd $HOME
    git clone https://github.com/PDB-REDO/dssp.git
    cd dssp
    mkdir -p build && cd build
    cmake ..
    make -j4
fi

if [ ! -d "libcifpp_cache" ]; then
    mkdir -p $HOME/libcifpp_cache
    curl -o ~/libcifpp_cache/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
    curl -o ~/libcifpp_cache/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic
    curl -o ~/libcifpp_cache/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic
fi

export PATH=$HOME/dssp/build:$PATH
ln -s $HOME/dssp/build/mkdssp $HOME/dssp/build/dssp
export LD_LIBRARY_PATH=$HOME/PROTAC-tor/ZDOCK:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/.local/lib/python3.11/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/rdkit/2024.09.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/icu73/lib:$LD_LIBRARY_PATH
export LIBCIFPP_DATA_DIR=~/libcifpp_cache

cd "{complex_dir}"

echo "Preprocessing PDB files..."
./mark_sur receptor.pdb Receptor.pdb
./mark_sur ligand.pdb Ligand.pdb

echo "Running ZDOCK..."
./zdock -R Ligand.pdb -L Receptor.pdb -o zdock_result.out
./create.pl zdock_result.out
""")
os.chmod(slurm_script_path, 0o755)

print("Submitting docking job...")
subprocess.run(["sbatch", str(slurm_script_path)])
print("Job submitted.")