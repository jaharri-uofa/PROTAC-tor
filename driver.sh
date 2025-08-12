#!/bin/bash

# === User settings ===
PYTHON=python3
SCRIPTS_DIR="$HOME/PROTAC-tor/scripts"
SLEEP_INTERVAL=60  # seconds between job checks

# Pipeline in order (Python scripts)
STEPS=(
    "prodock.py"
    "lig_dist.py"
    "link_it.py"
    "dock.py"
    "md.py"
    "md_mmgbsa.py"
)

# Helper: wait until no jobs with a certain name or user remain
wait_for_jobs() {
    echo "Waiting for jobs to finish..."
    while true; do
        # Get number of jobs for this user (change to match job name pattern if needed)
        num_jobs=$(squeue -u "$USER" --noheader | wc -l)
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
