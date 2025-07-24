#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
import stat
import subprocess
import logging as log
import rdkit
from rdkit import Chem

def remove_stereochemistry(smiles):
    '''
    Remove stereochemistry from a SMILES string.
    :param smiles: Input SMILES string
    :return: SMILES string without stereochemistry
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)

print("Starting ZDOCK docking automation...")

base_dir = Path.cwd()
zdock_dir = base_dir / "ZDOCK"
scripts_dir = base_dir / "scripts"
protein_complexes_dir = base_dir / "complexes"
protein_complexes_dir.mkdir(exist_ok=True)


# Get all .pdb files in the root of PROTAC-tor
pdb_files = list(base_dir.glob("*.pdb"))

# Make all pairwise combinations of pdb files
for i, pdb1 in enumerate(pdb_files):
    for j, pdb2 in enumerate(pdb_files):
        if i >= j:
            continue  # Avoid duplicates and self-pairing

        # Create a directory for the protein complex
        complex_name = f"{pdb1.stem}_{pdb2.stem}"
        complex_dir = protein_complexes_dir / complex_name
        complex_dir.mkdir(parents=True, exist_ok=True)

        # Copy required ZDOCK files from ZDOCK dir
        required_files = ["zdock", "create_lig", "create.pl", "mark_sur", "uniCHARMM", 'linkinvent.prior'] # move this to top or main
        for filename in required_files:
            src = zdock_dir / filename
            dest = complex_dir / filename
            shutil.copy(src, dest)
            dest.chmod(dest.stat().st_mode | stat.S_IEXEC)

        # Copy scripts
        for script_name in ["Lig_dist.py", "prodock.py", "link_it.py", "dock.py"]:
            shutil.copy(scripts_dir / script_name, complex_dir)

        # Copy receptor and ligand PDBs
        receptor_pdb = input(f"Enter POI PDB file for {complex_name}: ").strip()
        shutil.copy(receptor_pdb, complex_dir / "receptor.pdb")
        print(f"Copied {receptor_pdb} to {complex_dir / 'receptor.pdb'}")
        ligase_pdb = input(f"Enter E3 ligase PDB file for {complex_name}: ").strip()
        shutil.copy(ligase_pdb, complex_dir / "ligand.pdb")
        print(f"Copied {ligase_pdb} to {complex_dir / 'ligand.pdb'}")

        # Add dummy SEQRES
        (complex_dir / "SEQRES").write_text("DUMMYSEQRES\n")

        # (Optional) input SMILES strings for ligands
        lig1_smiles = remove_stereochemistry(input("Enter SMILES for ligand 1 (or leave empty to skip): ").strip())
        lig2_smiles = remove_stereochemistry(input("Enter SMILES for ligand 2 (or leave empty to skip): ").strip())
        if lig1_smiles and lig2_smiles:
            with open("smiles.smi", "w") as f:
                f.write(f"{lig1_smiles}|{lig2_smiles}\n")
        if not lig1_smiles or not lig2_smiles:
            print("Skipping SMILES input for LinkInvent.")
            os.remove("smiles.smi") if os.path.exists("smiles.smi") else None

        # Add smiles.smi if provided
        if lig1_smiles and lig2_smiles:
            shutil.copy("smiles.smi", complex_dir / "smiles.smi")

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

# probably a better way of doing this
if [ ! -d "xxhash" ]; then
    pip install xxhash
fi

# Build DSSP if not built already
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
export LD_LIBRARY_PATH=/home/jordanha/zdock_libs/usr/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/.local/lib/python3.11/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/rdkit/2024.09.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/icu73/lib:$LD_LIBRARY_PATH
export LIBCIFPP_DATA_DIR=~/libcifpp_cache


cd "{complex_dir}"

# Preprocess
echo "Preprocessing PDB files..."
./mark_sur receptor.pdb Receptor.pdb
./mark_sur ligand.pdb Ligand.pdb

# Run ZDOCK
echo "Running ZDOCK..."
./zdock -R Ligand.pdb -L Receptor.pdb -o zdock_result.out
./create.pl zdock_result.out

# Calculate distances
echo "Calculating distances..."
python Lig_dist.py

# LinkInvent
echo "Running LinkInvent..."
# reload modules because the cluster is goofy
module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.6
module load python-build-bundle/2025b
python link_it.py --smiles_csv smiles.smi --dist_file input.txt --output_toml staged_linkinvent.toml

# PRosettaC
echo "Running PRosettaC..."
python dock.py

# Molecular Dynamics
echo "Running Molecular Dynamics..."

# MM/GBSA
echo "Running MM/GBSA..."

""")
        os.chmod(slurm_script_path, 0o755)

print("Submitting all docking jobs...")
for sh_file in protein_complexes_dir.rglob("run_docking.sh"):
    subprocess.run(["sbatch", str(sh_file)])

print("All jobs submitted.")
