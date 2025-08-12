#!/bin/bash
#SBATCH --job-name=protactor_pipeline
#SBATCH --output=protactor_pipeline.log
#SBATCH --error=protactor_pipeline.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M

# Adjust if your environment needs module loading
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

# Path to the scripts directory
SCRIPTS_DIR="$HOME/PROTAC-tor/scripts"

# List of scripts in execution order
SCRIPTS=(
    "prodock.py"
    "Lig_dist.py"
    "link_it.py"
    "dock.py"
    "md.py"
    # Add in exact order you want them to run
)

# Function to submit each Python script as a SLURM job and wait for completion
submit_and_wait() {
    local script="$1"
    echo "Submitting job for: $script"

    # Create a temporary SLURM script
    JOB_SCRIPT=$(mktemp)
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=${script%.py}
#SBATCH --output=${script%.py}.out
#SBATCH --error=${script%.py}.err
#SBATCH --time=00:00:30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M

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
python "$SCRIPTS_DIR/$script"
EOF

    # Submit the job and capture the job ID
    job_id=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')
    echo "Submitted $script as job $job_id"

    # Wait for the job to complete
    while :
    do
        status=$(squeue --job "$job_id" --noheader | wc -l)
        if [ "$status" -eq 0 ]; then
            echo "Job $job_id completed."
            break
        fi
        sleep 30
    done

    rm "$JOB_SCRIPT"
}

# Main loop over all scripts
for script in "${SCRIPTS[@]}"; do
    submit_and_wait "$script"
done

echo "All scripts completed."
