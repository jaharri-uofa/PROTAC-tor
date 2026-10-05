#!/bin/bash
# PROTACtor Driver Script
#SBATCH --job-name=protactor
#SBATCH --output=protactor.out
#SBATCH --error=protactor.err
#SBATCH --mem=256M
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2023
module load python/3.11
module load scipy-stack/2023b
module load rdkit/2024.09.6
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.2
module load python-build-bundle/2025b
module load gnina/1.3.1
module load openmpi/4.1.5
module load amber-pmemd/24.3
module load ambertools/25.0

# === User settings ===
PYTHON=python3
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)}"
if [[ ! -d "$PROJECT_DIR" ]]; then
    echo "Project directory not found: $PROJECT_DIR" >&2
    exit 1
fi
cd -- "$PROJECT_DIR" || { echo "Failed to cd into project directory: $PROJECT_DIR" >&2; exit 1; }

SCRIPTS_DIR="$PROJECT_DIR/scripts"
SLEEP_INTERVAL=60  # seconds between job checks
SELF_JOB_ID="$SLURM_JOB_ID"  # Capture our own SLURM job ID

# Run prodock.py
echo "=== Running: prodock.py ==="
if [[ ! -f "$SCRIPTS_DIR/prodock.py" ]]; then
    echo "Required script not found: $SCRIPTS_DIR/prodock.py" >&2
    exit 1
fi
if ! "$PYTHON" "$SCRIPTS_DIR/prodock.py"; then
    echo "prodock.py failed; stopping before submitting docking jobs." >&2
    exit 1
fi

# Loop over all complexes/*/ directories and submit jobs
shopt -s nullglob
for target_dir in "$PROJECT_DIR"/complexes/*/; do
    if [[ ! -f "$target_dir/prodock.sh" ]]; then
        echo "Required job script not found: $target_dir/prodock.sh" >&2
        exit 1
    fi
    echo "Processing ${target_dir#"$PROJECT_DIR"/}"
    sbatch --chdir="$target_dir" "$target_dir/prodock.sh" || {
        echo "Failed to submit job script: $target_dir/prodock.sh" >&2
        exit 1
    }
done


