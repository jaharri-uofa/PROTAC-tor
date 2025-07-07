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

python run_linkinvent.py --config {output_json}
""")
    subprocess.run(["sbatch", slurm_script])
    print("Link-INVENT job submitted via SLURM.")

def extract_number(filename):
    match = re.search(r'complex\.(\d+)\.pdb', filename)
    return int(match.group(1)) if match else float('inf')

def get_main_ligand_id(pdb_file):
        """Return the ligand ID with the most atoms, skipping metal ions and solvent."""
        skip_residues = {
            'NA', 'CL', 'CA', 'MG', 'ZN', 'K', 'FE', 'CU', 'MN', 'HG',
            'HOH', 'WAT', 'SO4', 'PO4', 'HEM', 'DMS', 'ACE', 'NAG', 'GLC'
        }
        residue_atom_counts = {}
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    if resname not in skip_residues:
                        residue_atom_counts[resname] = residue_atom_counts.get(resname, 0) + 1
        return max(residue_atom_counts, key=residue_atom_counts.get) if residue_atom_counts else None

def main():
    """
    Loop through all PDB files in the current directory, find ligand pairs,
    compute closest distances, and run LinkInvent on the best complex.
    """
    print("Searching for ligand pairs in all PDB files...")
    
    current_dir = os.getcwd()
    receptor_pdb = os.path.join(current_dir, "receptor.pdb")
    ligand_pdb = os.path.join(current_dir, "ligand.pdb")
    lig1 = [get_main_ligand_id(receptor_pdb)]
    lig2 = [get_main_ligand_id(ligand_pdb)]
    print(f"Detected ligand 1 ID: {lig1}")
    print(f"Detected ligand 2 ID: {lig2}")

    teeny = {}
    ligand_ids = {}  # Track best matching ligand IDs per file
    cut = 20.0

    # This is dogshit
    for file in os.listdir('.'):
        if file.endswith('.pdb'):
            print(f"Processing {file}...")
            for lig in lig1:
                lig1_coords = get_lig(file, lig)
                if lig1_coords.size == 0:
                    continue
                print(f"Found ligand {lig} in {file}")
                for lig in lig2:
                    lig2_coords = get_lig(file, lig)
                    if lig2_coords.size == 0:
                        continue
                    for ligand1 in lig1_coords:
                        for ligand2 in lig2_coords:
                            dist = distance(ligand1, ligand2)
                            if dist < cut:
                                if file not in teeny or dist < teeny[file]:
                                    teeny[file] = dist
                                    ligand_ids[file] = (ligand1, ligand2)

    # Sort and truncate to top 3
    teeny = dict(sorted(teeny.items(), key=lambda item: item[1]))
    teeny = dict(list(teeny.items())[:10])
    lowest_files = sorted(teeny.keys(), key=extract_number)[:3]
    teeny = {k: teeny[k] for k in lowest_files}

    for pdb_file, dist in teeny.items():
        print(f"{pdb_file}: {dist:.2f} Å")

    print("\nTotal number of unique PDB files:", len(teeny))

    with open('lig_distances.txt', 'w') as f:
        for pdb_file, dist in teeny.items():
            f.write(f"{pdb_file}: {dist:.2f} Å\n")

    min_val, max_val = min_max('lig_distances.txt')
    print("Min distance:", min_val)
    print("Max distance:", max_val)

    # Get top complex and ligand IDs
    top_file = list(teeny.keys())[0]
    kinase_id, e3_id = ligand_ids[top_file]

    lig1_smiles = extract_ligand_smiles(top_file, kinase_id)
    lig2_smiles = extract_ligand_smiles(top_file, e3_id)

    print("Lig1 SMILES:", lig1_smiles)
    print("Lig2 SMILES:", lig2_smiles)

    with open('smiles.csv', 'w') as f:
        f.write("fragment_1,fragment_2\n")
        f.write(f"{lig1_smiles},{lig2_smiles}\n")

    with open('input.txt', 'w') as f:
        f.write(f"{min_val},{max_val}\n")

    LinkInvent()


main()