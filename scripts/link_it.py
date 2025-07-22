import argparse
import pandas as pd
import subprocess
import os
import toml
import rdkit
from rdkit import Chem

def molecule_features(smiles):
    """
    Extracts molecular features from a SMILES string.
    :param smiles: Input SMILES string
    :return: Dictionary of molecular features
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    features = {
        "MolecularWeight": Chem.rdMolDescriptors.CalcExactMolWt(mol),
        "TPSA": Chem.rdMolDescriptors.CalcTPSA(mol),
        "HBondAcceptors": Chem.rdMolDescriptors.CalcNumHBA(mol),
        "HBondDonors": Chem.rdMolDescriptors.CalcNumHBD(mol),
        "NumRotBond": Chem.rdMolDescriptors.CalcNumRotatableBonds(mol),
        "NumRings": Chem.rdMolDescriptors.CalcNumRings(mol),
        "NumAromaticRings": Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
        "LargestRingSize": Chem.rdMolDescriptors.CalcLargestRingSize(mol),
        "SAScore": Chem.rdMolDescriptors.CalcSAScore(mol),
        "SlogP": Chem.Crippen.MolLogP(mol)
    }
    return features

def generate_toml(smiles_csv, dist_file, output_toml):
    df = pd.read_csv(smiles_csv)

    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    chem_data = molecule_features(smiles_csv)

    weight = chem_data.get("MolecularWeight", 0)
    TPSA = chem_data.get("TPSA", 0)
    HBondAcceptors = chem_data.get("HBondAcceptors", 0)
    HBondDonors = chem_data.get("HBondDonors", 0)
    NumRotBond = chem_data.get("NumRotBond", 0)
    NumRings = chem_data.get("NumRings", 0)
    NumAromaticRings = chem_data.get("NumAromaticRings", 0)
    LargestRingSize = chem_data.get("LargestRingSize", 0)
    SlogP = chem_data.get("SlogP", 0)

    # Define a base stage
    base_stage = {
        "termination": "simple",
        "chkpt_file": "",  # Will be set per stage
        "max_score": 0.6,
        "min_steps": 100,
        "max_steps": 500,
        "scoring": {
            "type": "geometric_mean",
            "component": [
                {"MolecularWeight": {
                    "endpoint": [{
                        "name": "Molecular weight",
                        "weight": 1,
                        "transform": {
                            "type": "double_sigmoid",
                            "high": weight + 300.0,
                            "low": weight + 50.0,
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
                            "high": TPSA + 90.0,
                            "low": TPSA + 0,
                            "coef_div": 140.0,
                            "coef_si": 20.0,
                            "coef_se": 20.0
                        }
                    }]
                }},
            {"LinkerLength": {
                "endpoint": [{
                    "name": "Linker Effective Length",
                    "weight": 1,
                    "transform": {
                        "type": "double_sigmoid",
                        "low": min_dist,
                        "high": max_dist,
                        "coef_div": 20.0,
                        "coef_si": 1.0,
                        "coef_se": 1.0
                    }
                }]
            }},
            {"HBondAcceptors": {
                "endpoint": [{
                    "name": "Number of HB acceptors (Lipinski)",
                    "weight": 1,
                    "transform": {
                        "type": "reverse_sigmoid",
                        "high": HBondAcceptors + 5,
                        "low": HBondAcceptors + 0,
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
                            "high": HBondDonors + 1,
                            "low": HBondDonors + 0,
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
                            "high": NumRotBond + 20,
                            "low": NumRotBond + 5,
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
                            "high": NumRings + 1,
                            "low": NumRings + 0,
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
                            "high": NumAromaticRings + 1,
                            "low": NumAromaticRings + 0,
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
                            "high": LargestRingSize + 6,
                            "low": LargestRingSize + 5,
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
                            "high": SlogP + 5,
                            "low": SlogP + 2,
                            "k": 0.5
                        }
                    }]
                }}
            ]
        }
    }

    # Generate 3 training stages to increase total sampling
    stages = []
    for i in range(3):
        stage = base_stage.copy()
        stage["chkpt_file"] = f"stage{i+1}.chkpt"
        # You can increase step sizes if you want slower transitions
        stage["min_steps"] = 100 + i * 100
        stage["max_steps"] = 500 + i * 200
        stages.append(stage)

    config = {
        "run_type": "staged_learning",
        "device": "cuda:0",
        "tb_logdir": "tb_logs",
        "json_out_config": "staged_linkinvent.json",
        "parameters": {
            "summary_csv_prefix": "linkinvent_stage",
            "use_checkpoint": False,
            "purge_memories": False,
            "prior_file": "/home/jordanha/PROTAC-tor/complexes/4ci2_len_JNK3_36/linkinvent.prior",
            "agent_file": "/home/jordanha/PROTAC-tor/complexes/4ci2_len_JNK3_36/linkinvent.prior",
            "smiles_file": smiles_csv,
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
        "stage": stages
    }

    with open(output_toml, 'w') as f:
        toml.dump(config, f)

    print(f"Staged learning TOML configuration written to {output_toml}")

def write_slurm_script(output_toml, slurm_script="submit_linkinvent.sh"):
    slurm_contents = f"""#!/bin/bash
#SBATCH --job-name=linkinvent_gpu
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-08:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2023 gcc/12.3 cuda/12.6 python/3.11.5 python-build-bundle/2025b scipy-stack/2025a rdkit/2024.09.6 ipykernel/2025a

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
