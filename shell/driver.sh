#!/bin/bash
# PROTACtor Driver Script
DIR=$(dirname $(realpath $0))

#SBATCH --job-name=${DIR##*/}
#SBATCH --output=protactor.out
#SBATCH --error=protactor.err
#SBATCH --mem=128M
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.2
module load python-build-bundle/2025b
module load gnina/1.3.1
module load openmpi/4.1.5


# === lig_dist.py ===
echo "Running lig_dist.py"
python lig_dist.py

# === link_it.py ===
echo "Running link_it.py"
python link_it.py --smiles_csv smiles.smi --dist_file input.txt

# === link_it.sh ===
echo "Submitting link_it.sh to SLURM..."
link_jobid=$(sbatch --parsable link_it.sh)
echo "Submitted link_it.sh as job $link_jobid"

# === dock.py ===
echo "Submitting dock.py to SLURM after link_it.sh completes..."
dock_jobid=$(sbatch --parsable --dependency=afterok:$link_jobid --job-name=dockpy --output=dockpy.out --error=dockpy.err --wrap="python dock.py")
echo "Submitted dock.py as job $dock_jobid (after link_it.sh)"

# === md.py ===
echo "Submitting md.py to SLURM after dock.py completes..."
# important that these are loaded after all other jobs are finished due to module dependencies
module load ambertools/25.0
module load amber-pmemd/24.3
md_jobid=$(sbatch --parsable --dependency=afterok:$dock_jobid --job-name=mdpy --output=mdpy.out --error=mdpy.err --wrap="python md.py")
echo "Submitted md.py as job $md_jobid (after dock.py)"

# === analysis.py ===
echo "Submitting analysis.py to SLURM after md.py completes..."
analysis_jobid=$(sbatch --parsable --dependency=afterok:$md_jobid --job-name=analyzpy --output=analyzpy.out --error=analyzpy.err --wrap="python analysis.py")
echo "Submitted analysis.py as job $analysis_jobid (after md.py)"