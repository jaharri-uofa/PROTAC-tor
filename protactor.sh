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

# Pipeline in order (Python scripts)
STEPS=(
    "prodock.py"
    "Lig_dist.py"
    "link_it.py"
    "dock.py"
    "md.py"
    "md_mmgbsa.py"
    "analysis.py"
)

# Helper: wait until no jobs (other than our own) remain
wait_for_jobs() {
    echo "Waiting for jobs to start..."
    # Wait until a new job appears (other than this driver)
    while true; do
        num_jobs=$(squeue -u "$USER" --noheader | awk -v self="$SELF_JOB_ID" '$1 != self' | wc -l)
        if [[ "$num_jobs" -gt 0 ]]; then
            break
        fi
        sleep 5
    done

    echo "Jobs detected. Waiting for them to finish..."
    # Now wait until they're all done
    while true; do
        num_jobs=$(squeue -u "$USER" --noheader | awk -v self="$SELF_JOB_ID" '$1 != self' | wc -l)
        if [[ "$num_jobs" -eq 0 ]]; then
            echo "All jobs finished."
            break
        fi
        echo "$num_jobs jobs still running..."
        sleep "$SLEEP_INTERVAL"
    done
}


# Main pipeline loop
for step in "${STEPS[@]}"; do
    echo "=== Running: $step ==="
    $PYTHON "$SCRIPTS_DIR/$step"
    wait_for_jobs
done

echo "Pipeline complete!"
