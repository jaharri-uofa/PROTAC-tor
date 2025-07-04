'''
This script calculates the distance between two ligands in two input PDB files,
and extracts their SMILES representations for input into LinkInvent.
Author: Jordan Harrison
'''

import numpy as np
import os
from rdkit import Chem
import subprocess
import pandas as pd
import re
import sys

def distance(lig1, lig2):
    return np.linalg.norm(lig1 - lig2)

def get_lig(pdb_file, lig_id):
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and lig_id in line[17:20]:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
    return np.array(coords)

def get_ligand_id(pdb_file):
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                return line[17:20].strip()
    return None

def min_max(dist):
    return [dist, dist]  # Only one measurement

def extract_ligand_smiles(pdb_file, lig_id, output_name="ligand_only.pdb"):
    lig_atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and lig_id in line[17:20]:
                lig_atoms.append(line)
    if not lig_atoms or len(lig_atoms) < 3:
        print(f"Ligand {lig_id} in {pdb_file} is too small or missing — skipping.")
        return None
    with open(output_name, 'w') as f:
        f.writelines(lig_atoms)
        f.write("END\n")
    convert_with_obabel(output_name, "ligand.mol")
    mol = Chem.MolFromMolFile("ligand.mol", sanitize=True)
    if mol is None:
        print(f"RDKit failed to read {output_name}")
        return None
    return Chem.MolToSmiles(mol)

def convert_with_obabel(input_pdb, output_mol):
    try:
        subprocess.run([
            "obabel", input_pdb,
            "-O", output_mol,
            "-h", "--gen3d", "--AddPolarH",
            "--partialcharge", "gasteiger"
        ], check=True)
        print(f"Converted {input_pdb} → {output_mol}")
    except subprocess.CalledProcessError as e:
        print(f"Conversion failed: {e}")

def LinkInvent(smiles_csv='smiles.csv', dist_file='input.txt', output_json='linkinvent_config.json', slurm_script='submit_linkinvent.sh'):
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
            "columns": {"fragment_1": "Kinase_Ligand", "fragment_2": "E3_Ligand"}
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

python run_linkinvent.py --config {output_json}
""")
    subprocess.run(["sbatch", slurm_script])
    print("Link-INVENT job submitted via SLURM.")

def main():
    """
    Loop through all PDB files in the current directory, find ligand pairs,
    compute closest distances, and run LinkInvent on the best complex.
    """
    print("Searching for ligand pairs in all PDB files...")

    # Holds file: distance
    distances = {}
    ligand_pairs = {}

    # Iterate through all .pdb files in the current directory
    for pdb_file in os.listdir('.'):
        if not pdb_file.endswith('.pdb'):
            continue

        print(f"Processing: {pdb_file}")

        # Get all unique ligand IDs in the PDB file
        ligand_ids = set()
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    lig_id = line[17:20].strip()
                    if lig_id and lig_id not in ligand_ids:
                        ligand_ids.add(lig_id)
        ligand_ids = list(ligand_ids)

        '''
        # Skip if we don't find exactly 2 ligands
        if len(ligand_ids) != 2:
            print(f"Skipping {pdb_file}: Expected 2 ligands, found {len(ligand_ids)}")
            continue
        '''

        top_file = list(teeny.keys())[0]
        lig1, lig2 = ligand_ids[top_file]
        coords1 = get_lig(pdb_file, lig1)
        coords2 = get_lig(pdb_file, lig2)

        if coords1.size == 0 or coords2.size == 0:
            continue

        # Calculate min distance between any atom in lig1 and lig2
        min_dist = min(distance(a1, a2) for a1 in coords1 for a2 in coords2)
        distances[pdb_file] = min_dist
        ligand_pairs[pdb_file] = (lig1, lig2)

    if not distances:
        print("No valid complexes found.")
        return

    # Sort by distance and select top one
    sorted_files = sorted(distances.items(), key=lambda x: x[1])
    top_file, top_dist = sorted_files[0]
    lig1, lig2 = ligand_pairs[top_file]

    print(f"Top complex: {top_file}")
    print(f"Ligands: {lig1}, {lig2}")
    print(f"Distance: {top_dist:.2f} Å")

    # Extract SMILES
    smiles1 = extract_ligand_smiles(top_file, lig1)
    smiles2 = extract_ligand_smiles(top_file, lig2)

    if smiles1 is None or smiles2 is None:
        print("Failed to extract SMILES. Aborting.")
        return

    # Save SMILES to CSV
    with open('smiles.csv', 'w') as f:
        f.write("Kinase_Ligand,E3_Ligand\n")
        f.write(f"{smiles1},{smiles2}\n")

    # Save input.txt with tight bounds
    min_val = max(1.0, top_dist - 1)
    max_val = top_dist + 1.5
    with open('input.txt', 'w') as f:
        f.write(f"{min_val},{max_val}\n")

    # Run LinkInvent
    LinkInvent()


main()