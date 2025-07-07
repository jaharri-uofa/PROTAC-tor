def LinkInvent(smiles_csv='smiles.csv', dist_file='input.txt', output_json='linkinvent_config.json', slurm_script='submit_linkinvent.sh'):
    '''
    Generates a Link-INVENT configuration file and submits a SLURM job to run it.
    :param smiles_csv: Path to the CSV file containing SMILES strings
    :param dist_file: Path to the distance file
    :param output_json: Path to the output JSON file
    :param slurm_script: Path to the SLURM script file
    '''
    import json
    df = pd.read_csv(smiles_csv)
    fragment_1 = df.iloc[0, 0]
    fragment_2 = df.iloc[0, 1]
    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))
    config = {
        "logging": {"log_level": "INFO"},
        "run_type": "sample_linker",
        "input": {
            "source": smiles_csv,
            "columns": {"fragment_1": "Ligand1", "fragment_2": "Ligand2"}
        },
        "output": {"save_to": "linkinvent_output"},
        "scoring_function": {
            "name": "custom_sum",
            "parameters": [
                {"component_type": "LinkerLengthMatch", "name": "linker_length", "weight": 1, "specific_parameters": {"min_length": int(min_dist), "max_length": int(max_dist)}},
                {"component_type": "LinkerNumRings", "name": "max_one_ring", "weight": 1, "specific_parameters": {"min_num_rings": 0, "max_num_rings": 1}},
                {"component_type": "LinkerMW", "name": "mw_under_700", "weight": 1, "specific_parameters": {"min_mw": 0, "max_mw": 700}},
                {"component_type": "LinkerTPSA", "name": "tpsa_under_90", "weight": 1, "specific_parameters": {"min_tpsa": 0, "max_tpsa": 90}},
                {"component_type": "LinkerNumHBD", "name": "max_1_hbd", "weight": 1, "specific_parameters": {"min_hbd": 0, "max_hbd": 1}},
                {"component_type": "LinkerNumHBA", "name": "max_5_hba", "weight": 1, "specific_parameters": {"min_hba": 0, "max_hba": 5}},
                {"component_type": "LinkerLogP", "name": "logp_2_5", "weight": 1, "specific_parameters": {"min_logp": 2, "max_logp": 5}}
            ]
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
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=0-04:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2020  gcc/11.3.0
module load cuda/11.8.0
                
~/reinvent4/bin/pip install tomli
~/reinvent4/bin/pip install requests


python -m reinvent.runmodes.samplers.linkinvent --config linkinvent_config.json

""")
    subprocess.run(["sbatch", slurm_script])
    print("Link-INVENT job submitted via SLURM.")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Run Link-INVENT with specified parameters.")
    parser.add_argument('--smiles_csv', type=str, default='smiles.csv', help='Path to the CSV file containing SMILES strings.')
    parser.add_argument('--dist_file', type=str, default='input.txt', help='Path to the distance file.')
    parser.add_argument('--output_json', type=str, default='linkinvent_config.json', help='Path to the output JSON file.')
    parser.add_argument('--slurm_script', type=str, default='submit_linkinvent.sh', help='Path to the SLURM script file.')
    args = parser.parse_args()
    
    LinkInvent(args.smiles_csv, args.dist_file, args.output_json, args.slurm_script)