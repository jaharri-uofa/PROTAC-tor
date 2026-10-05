#!/bin/bash
#SBATCH --job-name=linkinvent_gpu
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:h100:1
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
module load scipy-stack/2023b
module load rdkit/2024.09.6
module load python-build-bundle/2025b

# === link_it.py ===
REINVENT_VENV="$HOME/reinvent4"
if [[ ! -f "$REINVENT_VENV/bin/activate" ]]; then
    echo "REINVENT virtual environment not found: $REINVENT_VENV" >&2
    exit 1
fi
source "$REINVENT_VENV/bin/activate"
export PYTHONNOUSERSITE=1
REINVENT_BIN="$REINVENT_VENV/bin/reinvent"
if [[ ! -x "$REINVENT_BIN" ]]; then
    echo "REINVENT executable not found in the virtual environment: $REINVENT_BIN" >&2
    echo "Install REINVENT into $REINVENT_VENV before submitting this job." >&2
    exit 1
fi

echo "Running REINVENT Link-INVENT sampling..."
"$REINVENT_BIN" -l staged.log sampling.toml
reinvent_status=$?
echo "REINVENT exit code: $reinvent_status"
if [[ "$reinvent_status" -ne 0 ]]; then
    echo "REINVENT failed; not submitting dock.py." >&2
    exit "$reinvent_status"
fi

link_jobid=$SLURM_JOB_ID

# === dock.py ===
echo "Submitting dock.py to SLURM after link_it.sh completes..."
dock_jobid=$(sbatch --parsable --dependency=afterok:$link_jobid --mem=2G --job-name=dockpy --output=dockpy.out --error=dockpy.err --wrap="python dock.py")
echo "Submitted dock.py as job $dock_jobid (after link_it.sh)"