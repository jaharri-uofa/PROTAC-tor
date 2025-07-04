'''
This script calculates the distance between two ligands in a pdb file, and extracts their SMILES representations.
Author: Jordan Harrison

Also for the love of god make some comments
'''

import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from rdkit.Chem import rdmolfiles, rdmolops
import pandas as pd
import re

def extract_number(filename):
    match = re.search(r'complex\.(\d+)\.pdb', filename)
    return int(match.group(1)) if match else float('inf')

def distance(lig1, lig2):
    """
    Calculate the distance between two ligands.
    :param lig1: First ligand coordinates as a numpy array.
    :param lig2: Second ligand coordinates as a numpy array.
    :return: Distance between the two ligands.
    """
    return np.linalg.norm(lig1 - lig2)

def get_lig(pdb_file, lig_id):
    """
    Extract the coordinates of a ligand from a PDB file.
    :param pdb_file: Path to the PDB file.
    :param lig_id: Identifier for the ligand (e.g., 'LIG').
    :return: Numpy array of ligand coordinates.
    """
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and lig_id in line[17:20]:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
    return np.array(coords)

def get_ligands(pdb_file):
    """
    Extract the ligand IDs from a PDB file.
    :param pdb_file: Path to the PDB file.
    """
    with open(pdb_file, 'r') as f:
        while lig_id is None:
            lig_id = None
            # Read each line in the PDB file
            for line in f:
                if line.startswith('HETATM'):
                    lig_id = line[17:20].strip()
                    if lig_id:
                        lig_id = lig_id.strip()
                        break

    return lig_id
               

def min_max(text):
    """
    Extract the minimum and maximum values from a text file containing distances.
    :param text: Path to the text file.
    """
    min_val = float('inf')
    max_val = float('-inf')
    with open(text, 'r') as f:
        for line in f:
            if ":" not in line:
                continue
            try:
                # extract value after the colon, before "Å"
                value = float(line.split(":")[1].split()[0])
                min_val = min(min_val, value)
                max_val = max(max_val, value)
            except Exception as e:
                print(f"Skipping line: {line.strip()} — {e}")
    return [min_val, max_val]

def extract_ligand_smiles(pdb_file, lig_id, output_name="ligand_only.pdb"):
    '''
    This function gave me an aneurysm, it handles extraction of ligands from PDB files only if they are decent, and uncharged.
    Charged ligands make this function explode 
    '''
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
    mol = Chem.MolToSmiles(mol)
    print(mol)

    if mol is None:
        print(f"RDKit failed to read {output_name}")
        return None
    
    return mol

def convert_with_obabel(input_pdb, output_mol):
    ''' 
    Convert a PDB file to a MOL file using Open Babel.
    :param input_pdb: Path to the input PDB file.
    :param output_mol: Path to the output MOL file.
    '''
    try:
        subprocess.run([
            "obabel", input_pdb,
            "-O", output_mol,
            "-h",           # Add hydrogens
            "--gen3d",      # Generate 3D coordinates (if needed)
            "--AddPolarH",  # Add polar hydrogens
            "--partialcharge", "gasteiger"  # Add partial charges
        ], check=True)
        print(f"Converted {input_pdb} → {output_mol}")
    except subprocess.CalledProcessError as e:
        print(f"Conversion failed: {e}")

def LinkInvent(smiles_csv='smiles.csv', dist_file='input.txt', output_json='linkinvent_config.json', slurm_script='submit_linkinvent.sh'):
    import json
    import subprocess
    import pandas as pd

    # Load SMILES
    df = pd.read_csv(smiles_csv)
    fragment_1 = df.iloc[0, 0]
    fragment_2 = df.iloc[0, 1]

    # Load distance limits
    with open(dist_file, 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    # Build config
    config = {
        "logging": {"log_level": "INFO"},
        "run_type": "sample_linker",
        "input": {
            "source": smiles_csv,
            "columns": {
                "fragment_1": "Kinase_Ligand",
                "fragment_2": "E3_Ligand"
            }
        },
        "output": {"save_to": "linkinvent_output"},
        "scoring_function": {
            "name": "custom_sum",
            "parameters": [
                {"component_type": "LinkerLengthMatch", "name": "linker_length", "weight": 1,
                 "specific_parameters": {"min_length": int(min_dist), "max_length": int(max_dist)}},
                {"component_type": "LinkerNumRings", "name": "max_one_ring", "weight": 1,
                 "specific_parameters": {"min_num_rings": 0, "max_num_rings": 1}},
                {"component_type": "LinkerMW", "name": "mw_under_700", "weight": 1,
                 "specific_parameters": {"min_mw": 0, "max_mw": 700}},
                {"component_type": "LinkerTPSA", "name": "tpsa_under_90", "weight": 1,
                 "specific_parameters": {"min_tpsa": 0, "max_tpsa": 90}},
                {"component_type": "LinkerNumHBD", "name": "max_1_hbd", "weight": 1,
                 "specific_parameters": {"min_hbd": 0, "max_hbd": 1}},
                {"component_type": "LinkerNumHBA", "name": "max_5_hba", "weight": 1,
                 "specific_parameters": {"min_hba": 0, "max_hba": 5}},
                {"component_type": "LinkerLogP", "name": "logp_2_5", "weight": 1,
                 "specific_parameters": {"min_logp": 2, "max_logp": 5}}
            ]
        }
    }

    with open(output_json, 'w') as f:
        json.dump(config, f, indent=4)
    print(f"Link-INVENT config written to: {output_json}")

    # SLURM script content
    slurm_content = f"""#!/bin/bash
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

# Run LinkInvent
python run_linkinvent.py --config {output_json}
"""

    # Write SLURM script
    with open(slurm_script, 'w') as f:
        f.write(slurm_content)

    print(f"SLURM script written to: {slurm_script}")

    # Submit job
    subprocess.run(["sbatch", slurm_script])
    print("Link-INVENT job submitted via SLURM.")

def main():
    """
    Need to add some logic to extract ligand IDs from the PDB files. So you dont have the massive chain of if statements.
    probably make the code run better.
    """
    Kinase_lig = ['LIG', "*"]
    E3_ligs = ['LVY']
    teeny = {}
    ligand_ids = {}  # Track best matching ligand IDs per file
    cut = 20.0

    # This is dogshit
    for file in os.listdir('.'):
        if file.endswith('.pdb'):
            print(f"Processing {file}...")
            for kinase_lig in Kinase_lig:
                kinase_lig_coords = get_lig(file, kinase_lig)
                if kinase_lig_coords.size == 0:
                    continue
                print(f"Found kinase ligand {kinase_lig} in {file}")
                for E3_lig_id in E3_ligs:
                    E3_lig_coords = get_lig(file, E3_lig_id)
                    if E3_lig_coords.size == 0:
                        continue
                    for lig1 in kinase_lig_coords:
                        for lig2 in E3_lig_coords:
                            dist = distance(lig1, lig2)
                            if dist < cut:
                                if file not in teeny or dist < teeny[file]:
                                    teeny[file] = dist
                                    ligand_ids[file] = (kinase_lig, E3_lig_id)

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

    kinase_lig_smiles = extract_ligand_smiles(top_file, kinase_id)
    e3_lig_smiles = extract_ligand_smiles(top_file, e3_id)

    print("Kinase SMILES:", kinase_lig_smiles)
    print("E3 SMILES:", e3_lig_smiles)

    with open('smiles.csv', 'w') as f:
        f.write("Kinase_Ligand,E3_Ligand\n")
        f.write(f"{kinase_lig_smiles},{e3_lig_smiles}\n")

    with open('input.txt', 'w') as f:
        f.write(f"{min_val},{max_val}\n")

    LinkInvent()


main()