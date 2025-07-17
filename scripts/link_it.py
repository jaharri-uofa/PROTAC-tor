import toml
import argparse
import pandas as pd
import subprocess
import os

def generate_toml(smiles_csv, dist_file, output_toml):
    df = pd.read_csv(smiles_csv)

    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    config = {
        "run_type": "sampling",
        "json_out_config": "sampling_config.json",
        "device": "cuda:0",
        "parameters": {
            "model_file": "/home/jordanha/REINVENT4/priors/linkinvent.prior",
            "smiles_file": smiles_csv,
            "sample_strategy": "multinomial",
            "output_file": "sampled_linkers.csv",
            "num_smiles": 1000,
            "unique_molecules": True,
            "randomize_smiles": False
        }
    }

    with open(output_toml, 'w') as f:
        toml.dump(config, f)

    print(f"TOML configuration written to {output_toml}")

def write_slurm_script(output_toml, slurm_script="submit_linkinvent.sh"):
    slurm_contents = f"""#!/bin/bash
#SBATCH --job-name=linkinvent
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=0-00:30
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
module load gcc/13.3
module load cmake
module load cuda/12.6
python-build-bundle/2025b

source ~/reinvent4/bin/activate

~/reinvent4/bin/pip install tomli requests numpy --no-index -f /cvmfs/soft.computecanada.ca/custom/python/wheelhouse/gentoo2023/generic

echo "Python being used:"
~/reinvent4/bin/python -c "import sys; print(sys.executable)"
echo "Loaded modules:"
~/reinvent4/bin/pip list | grep -E

set -euo pipefail

echo "Running REINVENT Link-INVENT sampling..."
reinvent -l sampling.log {output_toml}

echo "Exit code: $?"
echo "Output folder contents:"
ls -lh
"""
    with open(slurm_script, 'w') as f:
        f.write(slurm_contents)

    print(f"SLURM script written to {slurm_script}")

def submit_job(slurm_script="submit_linkinvent.sh"):
    print("Submitting job with sbatch...")
    subprocess.run(["sbatch", slurm_script])
    print("Job submitted.")

def main():
    parser = argparse.ArgumentParser(description="Build and submit a Link-INVENT REINVENT sampling job using TOML + SLURM.")
    parser.add_argument("--smiles_csv", required=True, help="Input SMILES CSV (with fragment_1, fragment_2).")
    parser.add_argument("--dist_file", required=True, help="Distance file containing min,max.")
    parser.add_argument("--output_toml", default="sampling.toml", help="Output TOML config file.")
    parser.add_argument("--slurm_script", default="submit_linkinvent.sh", help="SLURM script filename.")
    args = parser.parse_args()

    assert os.path.exists(args.smiles_csv), f"Missing SMILES file: {args.smiles_csv}"
    assert os.path.exists(args.dist_file), f"Missing distance file: {args.dist_file}"

    generate_toml(args.smiles_csv, args.dist_file, args.output_toml)
    write_slurm_script(args.output_toml, args.slurm_script)
    submit_job(args.slurm_script)

main()
