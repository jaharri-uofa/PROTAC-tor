#!/bin/bash
# PROTACtor Driver Script
#SBATCH --job-name=protactor
#SBATCH --output=protactor.out
#SBATCH --error=protactor.err
#SBATCH --mem=256M
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
module load amber/22.5-23.5


# === User settings ===
PYTHON=python3
SCRIPTS_DIR="$HOME/PROTAC-tor/scripts"
SLEEP_INTERVAL=60  # seconds between job checks
SELF_JOB_ID="$SLURM_JOB_ID"  # Capture our own SLURM job ID

# Run prodock.py
echo "=== Running: prodock.py ==="
$PYTHON "$SCRIPTS_DIR/prodock.py"

# Loop over all complexes/*/ directories and submit jobs
for target_dir in complexes/*/; do
    if [[ -d "$target_dir" ]]; then
        echo "Processing $target_dir"
        cd "$target_dir" || { echo "Failed to cd into $target_dir"; exit 1; }
        sbatch prodock.sh
        cd - >/dev/null
    fi
done


