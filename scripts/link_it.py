import json
import pandas as pd
import subprocess
import xxhash
import argparse
import os
import numpy as np

def LinkInvent(smiles_csv='smiles.csv', dist_file='input.txt', output_json='linkinvent_config.json', minimal_json='linkinvent_minimal.json', slurm_script='submit_linkinvent.sh'):
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
    "version": 4,  
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
            {"component_type": "LinkerLengthMatch","name": "linker_length","weight": 1,"specific_parameters": {"min_length": int(min_dist),"max_length": int(max_dist)}},
            {"component_type": "LinkerNumRings", "name": "max_one_ring", "weight": 1, "specific_parameters": {"min_num_rings": 0, "max_num_rings": 1}},
            {"component_type": "LinkerMW", "name": "mw_under_700", "weight": 1, "specific_parameters": {"min_mw": 0, "max_mw": 700}},
            {"component_type": "LinkerTPSA", "name": "tpsa_under_90", "weight": 1, "specific_parameters": {"min_tpsa": 0, "max_tpsa": 90}},
            {"component_type": "LinkerNumHBD", "name": "max_1_hbd", "weight": 1, "specific_parameters": {"min_hbd": 0, "max_hbd": 1}},
            {"component_type": "LinkerNumHBA", "name": "max_5_hba", "weight": 1, "specific_parameters": {"min_hba": 0, "max_hba": 5}},
            {"component_type": "LinkerLogP", "name": "logp_2_5", "weight": 1, "specific_parameters": {"min_logp": 2, "max_logp": 5}}
        ]
    },
    "reinforcement_learning": {
        "enabled": False,  # Set to True if using RL
        "prior_weight": 0.5,
        "sigma": 128
        }
    }
    
    print("Validating configuration...")
    print(json.dumps(config, indent=4))
    with open(output_json, 'w') as f:
        json.dump(config, f, indent=4)
    print(f"Link-INVENT config written to: {output_json}")
    with open(slurm_script, 'w') as f:
        f.write(f"""#!/bin/bash
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
                
module load StdEnv/2023 gcc/12.3 cuda/12.6 python/3.11.5 python-build-bundle/2024a scipy-stack/2024b rdkit/2024.03.4
source ~/reinvent4/bin/activate
                
~/reinvent4/bin/pip install tomli
~/reinvent4/bin/pip install requests
~/reinvent4/bin/pip install numpy
~/reinvent4/bin/pip install --no-index typing_extensions \
  -f /cvmfs/soft.computecanada.ca/custom/python/wheelhouse/gentoo2023/generic
                
echo "Python being used:"
~/reinvent4/bin/python -c "import sys; print(sys.executable)"
echo "Loaded modules:"
~/reinvent4/bin/pip list | grep -E

                
# Run with error trapping
set -euo pipefail

echo "Running Link-INVENT..."
~/reinvent4/bin/python -m reinvent.runmodes.samplers.linkinvent --config {output_json} 2>&1 | tee linkinvent_run.log
#  ~/reinvent4/bin/python -m reinvent.runmodes.samplers.linkinvent --config linkinvent_config.json 2>&1 | tee linkinvent_run.log

# Post-run validation
echo "Exit code: $?"
echo "Output folder contents:"
ls -lh linkinvent_output

""")
    subprocess.run(["sbatch", slurm_script])
    print("Link-INVENT job submitted via SLURM.")

def main():
    print("executing LinkInvent script...")

    parser = argparse.ArgumentParser(description="Run Link-INVENT with specified parameters.")
    parser.add_argument('--smiles_csv', type=str, default='smiles.csv', help='Path to the CSV file containing SMILES strings.')
    parser.add_argument('--dist_file', type=str, default='input.txt', help='Path to the distance file.')
    parser.add_argument('--output_json', type=str, default='linkinvent_config.json', help='Path to the output JSON file.')
    parser.add_argument('--minimal_json', type=str, default='linkinvent_minimal.json', help='Path to the minimal output JSON file.')
    parser.add_argument('--slurm_script', type=str, default='submit_linkinvent.sh', help='Path to the SLURM script file.')
    args = parser.parse_args()
    
    # Validate input files
    assert os.path.exists(args.smiles_csv), f"Input file {args.smiles_csv} not found"
    assert os.path.exists(args.dist_file), f"Distance file {args.dist_file} not found"

    LinkInvent(args.smiles_csv, args.dist_file, args.output_json, args.minimal_json, args.slurm_script)

main()