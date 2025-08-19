#!/bin/bash
#SBATCH --job-name=prodock
#SBATCH --output=prodock.out
#SBATCH --error=prodock.err
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

# need to find out which of these modules are necessary for the scripts, and then pre-install them as narval doesnt have internet access
# use nibi as a reference point, and install any missing modules, and make sure you have certain things in path?
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

# Begin zdock protein preparation/docking
cd "complexes/*/"

echo "Preprocessing PDB files..."
./mark_sur receptor.pdb Receptor.pdb
./mark_sur ligand.pdb Ligand.pdb

echo "Running ZDOCK..."
./zdock -R Ligand.pdb -L Receptor.pdb -o zdock_result.out
./create.pl zdock_result.out

sbatch driver.sh