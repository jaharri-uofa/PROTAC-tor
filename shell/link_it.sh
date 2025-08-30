#!/bin/bash
#SBATCH --job-name=linkinvent_gpu
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module --force purge 
module load StdEnv/2023
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.6
module load python/3.11.5
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load python-build-bundle/2025b

# === link_it.py ===
source ~/reinvent4/bin/activate
export PATH=$HOME/.local/bin:$PATH

echo "Running REINVENT Link-INVENT sampling..."
reinvent -l staged.log sampling.toml

echo "Exit code: $?"

link_jobid=$SLURM_JOB_ID

# === dock.py ===
echo "Submitting dock.py to SLURM after link_it.sh completes..."
dock_jobid=$(sbatch --parsable --dependency=afterok:$link_jobid --job-name=dockpy --output=dockpy.out --error=dockpy.err --wrap="python dock.py")
echo "Submitted dock.py as job $dock_jobid (after link_it.sh)"