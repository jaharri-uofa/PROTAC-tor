'''
Builds and submits a linkinvent script for de novo linker design using REINVENT 4.0
Author: Jordan Harrison
'''

import argparse
import pandas as pd
import subprocess
import os
import toml
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
import re

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
        "SlogP": Chem.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    }
    return features

def extract_warhead_smiles(smiles):
    '''
    Extracts the warhead SMILES from the full PROTAC SMILES.
    Assumes the format is "warhead1|warhead2".
    :param smiles: smiles string in the format "warhead1|warhead2"
    :return: Tuple of warhead SMILES strings (warhead1, warhead2)
    '''
    parts = smiles.split('|')
    if len(parts) != 2:
        raise ValueError(f"Invalid PROTAC SMILES format: {smiles}")
    return parts[0].strip(), parts[1].strip()

def longest_path_length(smiles: str) -> int:
    '''
    Calculates the longest path length (number of bonds) in a molecule from its SMILES representation.
    :param smiles: SMILES string of the molecule
    :return: length of the longest path in bonds
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    dmat = rdmolops.GetDistanceMatrix(mol)
    return int(dmat.max())  # number of bonds in longest path

# need to find a way to make this more generalizable, just including distance is insufficient
# may also be worth it to have this as a seperate TOML file?
def generate_toml(smiles_csv, dist_file, output_toml):
    '''
    Generates a TOML configuration file for REINVENT 4.0 based on input SMILES and distance data.
    :param smiles_csv: Path to the input CSV file containing SMILES strings.
    :param dist_file: Path to the input distance file.
    :param output_toml: Path to the output TOML file.
    '''
    # df = pd.read_csv(smiles_csv)

    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    with open(smiles_csv, 'r') as f:
        smiles1, smiles2 = extract_warhead_smiles(f.readline().strip())

    chem_data = [molecule_features(smiles) for smiles in [smiles1, smiles2]]

    # hang onto this for future use
    weight = sum(data.get("MolecularWeight", 0) for data in chem_data)
    TPSA = sum(data.get("TPSA", 0) for data in chem_data)
    HBondAcceptors = sum(data.get("HBondAcceptors", 0) for data in chem_data)
    HBondDonors = sum(data.get("HBondDonors", 0) for data in chem_data)
    NumRotBond = sum(data.get("NumRotBond", 0) for data in chem_data)
    NumRings = sum(data.get("NumRings", 0) for data in chem_data)
    NumAromaticRings = sum(data.get("NumAromaticRings", 0) for data in chem_data)
    LargestRingSize = sum(data.get("LargestRingSize", 0) for data in chem_data)
    SlogP = sum(data.get("SlogP", 0) for data in chem_data)
    length = longest_path_length(smiles1) + longest_path_length(smiles2)

    # Define a base stage
    base_stage = {
        "termination": "simple",
        "chkpt_file": "",  # Will be set per stage
        "max_score": 0.6,
        "min_steps": 1000,
        "max_steps": 10000,
        "scoring": {
            "type": "geometric_mean",
            "component": [
                {"MolecularWeight": {
                    "endpoint": [{
                        "name": "Molecular weight",
                        "weight": 1,
                        "transform": {
                            "type": "reverse_sigmoid",
                            "high": weight + int(max_dist) * 2.0,
                            "low": weight + int(min_dist) * 1.5,
                            "k": 0.5
                        }
                    }]
                }},
                {"FragmentGraphLength": {
                    "endpoint": [{
                        "name": "Molecule length (number of bonds in longest path)",
                        "weight": 1,
                        "transform": {
                            "type": "sigmoid",
                            "high": int(max_dist) * 1.5,
                            "low": int(min_dist) * 1.5,
                            "k": 0.5
                        }
                    }]
                }},
                {"FragmentEffectiveLength": {
                        "endpoint": [{
                            "name": "Effective length (distance between anchor atoms)",
                            "weight": 1,
                            "transform": {
                                "type": "sigmoid",
                                "high": int(max_dist) + 6,
                                "low": int(min_dist) + 4,
                                "k": 0.5
                            }
                        }]
                    }
                },
                {
                    "FragmentLengthRatio": {
                        "endpoint": [{
                            "name": "Length ratio (effective / graph length)",
                            "weight": 1,
                            "transform": {
                                "type": "sigmoid",
                                "high": 1.0,
                                "low": 0.99,
                                "k": 0.5
                            }
                        }]
                    }
                },
                {"FragmentTPSA": {
                    "endpoint": [{
                        "name": "TPSA",
                        "weight": 1,
                        "transform": {
                            "type": "sigmoid",
                            "high": + 90.0,
                            "low": + 30.0,
                            "k": 0.5
                        }
                    }]
                }},
            {"FragmentHBondAcceptors": {
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
                {"FragmentHBondDonors": {
                    "endpoint": [{
                        "name": "Number of HB donors (Lipinski)",
                        "weight": 1,
                        "transform": {
                            "type": "sigmoid",
                            "high": 1,
                            "low": 0,
                            "k": 0.5
                        }
                    }]
                }},
                {"FragmentNumRotBond": {
                    "endpoint": [{
                        "name": "Number of rotatable bonds",
                        "weight": 1,
                        "transform": {
                            "type": "sigmoid",
                            "high": 10,
                            "low": 5,
                            "k": 0.5
                        }
                    }]
                }},
                {"FragmentNumRings": {
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
                {"FragmentNumAromaticRings": {
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
    }

    # Generate 3 training stages to increase total sampling
    stages = []
    for i in range(3):
        stage = base_stage.copy()
        stage["chkpt_file"] = f"stage{i+1}.chkpt"
        # You can increase step sizes if you want slower transitions
        stage["min_steps"] = 1000 + i * 1000
        stage["max_steps"] = 5000 + i * 2000
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
            "prior_file": "linkinvent.prior",
            "agent_file": "linkinvent.prior",
            "smiles_file": smiles_csv,
            "batch_size": 64,
            "unique_sequences": True,
            "randomize_smiles": True,
            "tb_isim": False
        },
        "learning_strategy": {
            "type": "dap",
            "sigma": 256,
            "rate": 0.0001
        },
        "diversity_filter": {
            "type": "ScaffoldSimilarity",
            "bucket_size": 100,
            "minscore": 0.4,
            "minsimilarity": 0.3,
            "penalty_multiplier": 1.0
        },
        "stage": stages
    }

    with open(output_toml, 'w') as f:
        toml.dump(config, f)

    print(f"Staged learning TOML configuration written to {output_toml}")

# need to remove this at some point
def write_slurm_script(output_toml, slurm_script="submit_linkinvent.sh"):
    slurm_contents = f"""#!/bin/bash
#SBATCH --job-name=linkinvent_gpu
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:30
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2023
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.6
module load python/3.11.5
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load python-build-bundle/2025b

echo "Running REINVENT Link-INVENT sampling..."
reinvent -l staged.log {output_toml}

echo "Exit code: $?"
echo "Job ran succesfully"
"""
    with open(slurm_script, 'w') as f:
        f.write(slurm_contents)

    print(f"SLURM script written to {slurm_script}")

# ditto
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

    generate_toml((args.smiles_csv), args.dist_file, args.output_toml)
    write_slurm_script(args.output_toml, args.slurm_script)

main()
