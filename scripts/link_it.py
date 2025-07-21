import argparse
import pandas as pd
import subprocess
import os
import toml

def generate_toml(smiles_csv, dist_file, output_toml):
    df = pd.read_csv(smiles_csv)

    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    config = {
        "run_type": "staged_learning",
        "device": "cpu",
        "tb_logdir": "tb_logs",
        "json_out_config": "staged_linkinvent.json",
        "parameters": {
            "summary_csv_prefix": "linkinvent_stage",
            "use_checkpoint": False,
            "purge_memories": False,
            "prior_file": "/home/jordanha/PROTAC-tor/complexes/4ci2_len_JNK3_36/linkinvent.prior",
            "agent_file": "/home/jordanha/PROTAC-tor/complexes/4ci2_len_JNK3_36/linkinvent.prior",
            "smiles_file": smiles_csv,  # Format: warhead1_SMILES|warhead2_SMILES
            "batch_size": 64,
            "unique_sequences": True,
            "randomize_smiles": True,
            "tb_isim": False
        },
        "learning_strategy": {
            "type": "dap",
            "sigma": 128,
            "rate": 0.0001
        },
        "diversity_filter": {
            "type": "IdenticalMurckoScaffold",
            "bucket_size": 25,
            "minscore": 0.4,
            "minsimilarity": 0.4,
            "penalty_multiplier": 0.5
        },
        "stage": [{
            "chkpt_file": "stage1.chkpt",
            "termination": "simple",
            "max_score": 0.6,
            "min_steps": 50,
            "max_steps": 200,
            "scoring": {
                "type": "geometric_mean",
                "component": [
                    {"MolecularWeight": {
                        "endpoint": [{
                            "name": "Molecular weight",
                            "weight": 1,
                            "transform": {
                                "type": "double_sigmoid",
                                "high": 300.0,
                                "low": 50.0,
                                "coef_div": 500.0,
                                "coef_si": 20.0,
                                "coef_se": 20.0
                            }
                        }]
                    }},
                    {"TPSA": {
                        "endpoint": [{
                            "name": "TPSA",
                            "weight": 1,
                            "transform": {
                                "type": "double_sigmoid",
                                "high": 90.0,
                                "low": 0.0,
                                "coef_div": 140.0,
                                "coef_si": 20.0,
                                "coef_se": 20.0
                            }
                        }]
                    }},
                    {"HBondAcceptors": {
                        "endpoint": [{
                            "name": "Number of HB acceptors (Lipinski)",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 5,
                                "low": 0,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"HBondDonors": {
                        "endpoint": [{
                            "name": "Number of HB donors (Lipinski)",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 1,
                                "low": 0,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"NumRotBond": {
                        "endpoint": [{
                            "name": "Number of rotatable bonds",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 20,
                                "low": 5,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"NumRings": {
                        "endpoint": [{
                            "name": "Number of rings",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 1,
                                "low": 0,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"NumAromaticRings": {
                        "endpoint": [{
                            "name": "Number of aromatic rings",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 1,
                                "low": 0,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"NumAliphaticRings": {
                        "endpoint": [{
                            "name": "Number of aliphatic rings",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 1,
                                "low": 0,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"LargestRingSize": {
                        "endpoint": [{
                            "name": "Number of atoms in the largest ring",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 6,
                                "low": 5,
                                "k": 0.5
                            }
                        }]
                    }},
                    {"SAScore": {
                        "endpoint": [{
                            "name": "SA score",
                            "weight": 1
                        }]
                    }},
                    {"SlogP": {
                        "endpoint": [{
                            "name": "SlogP (RDKit)",
                            "weight": 1,
                            "transform": {
                                "type": "reverse_sigmoid",
                                "high": 5,
                                "low": 2,
                                "k": 0.5
                            }
                        }]
                    }}
                ]
            }
        }]
    }

    with open(output_toml, 'w') as f:
        toml.dump(config, f)

    print(f"Staged learning TOML configuration written to {output_toml}")

def write_slurm_script(output_toml, slurm_script="submit_linkinvent.sh"):
    slurm_contents = f"""#!/bin/bash
#SBATCH --job-name=linkinvent
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
##SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=3
#SBATCH --time=0-08:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module purge
module load StdEnv/2023 gcc/12.3 cuda/12.6 python/3.11.5 python-build-bundle/2025b scipy-stack/2025a rdkit/2024.09.6 ipykernel/2025a
source ~/reinvent4/bin/activate


set -euo pipefail

echo "Running REINVENT Link-INVENT sampling..."
reinvent -l staged.log {output_toml}

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
