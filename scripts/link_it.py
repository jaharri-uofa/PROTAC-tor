import json
import pandas as pd
import subprocess
import xxhash

def LinkInvent(smiles_csv='smiles.csv', dist_file='input.txt', output_json='linkinvent_config.json', slurm_script='submit_linkinvent.sh'):
    '''
    Generates a Link-INVENT configuration file and submits a SLURM job to run it.
    :param smiles_csv: Path to the CSV file containing SMILES strings
    :param dist_file: Path to the distance file
    :param output_json: Path to the output JSON file
    :param slurm_script: Path to the SLURM script file
    '''
    df = pd.read_csv(smiles_csv)
    fragment_1 = df.iloc[0, 0]
    fragment_2 = df.iloc[0, 1]
    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))
    config = {
    "version": 3,  # Important for REINVENT/LINKinvent compatibility
    "logging": {
        "log_level": "DEBUG",
        "log_file": "linkinvent_debug.log",
        "logging_path": "logs"
    },
    "model": {
        "path": "/home/jordanha/REINVENT4/priors/linkinvent.prior",
        "type": "LinkInvent",
        "model_parameters": {
            "batch_size": 64,
            "learning_rate": 0.0001,
            "num_epochs": 50,
            "hidden_size": 512,
            "num_layers": 3
        }
    },
    "run_type": "sample_linker",
    "input": {
        "source": smiles_csv,
        "columns": {
            "fragment_1": "fragment_1",
            "fragment_2": "fragment_2",
            "distance": "distance"  # Add if your CSV has distance column
        },
        "delimiter": ","
    },
    "output": {
        "save_to": "linkinvent_output",
        "save_every_n_epochs": 5,
        "num_samples": 1000
    },
    "scoring_function": {
        "name": "custom_sum",
        "parameters": [
            {
                "component_type": "LinkerLengthMatch",
                "name": "linker_length",
                "weight": 1,
                "specific_parameters": {
                    "min_length": int(min_dist),
                    "max_length": int(max_dist)
                }
            },
            # ... [keep your other scoring components] ...
        ]
    },
    "reinforcement_learning": {
        "enabled": False,  # Set to True if using RL
        "prior_weight": 0.5,
        "sigma": 128
        }
    }
    with open(output_json, 'w') as f:
        json.dump(config, f, indent=4)
    print(f"Link-INVENT config written to: {output_json}")
    with open(slurm_script, 'w') as f:
        f.write(f"""#!/bin/bash
#SBATCH --job-name=linkinvent
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=0-00:30
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca
                
~/reinvent4/bin/pip install tomli
~/reinvent4/bin/pip install requests
~/reinvent4/bin/pip install --no-index typing_extensions \
  -f /cvmfs/soft.computecanada.ca/custom/python/wheelhouse/gentoo2023/generic
                
echo "Python being used:"
~/reinvent4/bin/python -c "import sys; print(sys.executable)"

module load cuda/11.4
module load python/3.8
nvidia-smi
echo "CUDA version:"
nvcc --version

# Validate input files
echo "Checking input files:"
ls -lh ${SLURM_SUBMIT_DIR}/linkinvent_config.json
ls -lh ${SLURM_SUBMIT_DIR}/${smiles_csv}

# Run with error trapping
set -euo pipefail

echo "Running Link-INVENT..."
~/reinvent4/bin/python -m reinvent.runmodes.samplers.linkinvent \
    --config ${SLURM_SUBMIT_DIR}/linkinvent_config.json \
    2>&1 | tee linkinvent_run.log

# Post-run validation
echo "Exit code: $?"
echo "Output folder contents:"
ls -lh linkinvent_output

""")
    subprocess.run(["sbatch", slurm_script])
    print("Link-INVENT job submitted via SLURM.")

def main():
    print("executing LinkInvent script...")
    import argparse
    parser = argparse.ArgumentParser(description="Run Link-INVENT with specified parameters.")
    parser.add_argument('--smiles_csv', type=str, default='smiles.csv', help='Path to the CSV file containing SMILES strings.')
    parser.add_argument('--dist_file', type=str, default='input.txt', help='Path to the distance file.')
    parser.add_argument('--output_json', type=str, default='linkinvent_config.json', help='Path to the output JSON file.')
    parser.add_argument('--slurm_script', type=str, default='submit_linkinvent.sh', help='Path to the SLURM script file.')
    args = parser.parse_args()
    
    LinkInvent(args.smiles_csv, args.dist_file, args.output_json, args.slurm_script)

main()