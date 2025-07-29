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
'''
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import make_dssp_dict
'''
from pathlib import Path

def distance(lig1, lig2):
    '''
    Calculates distance between two ligands
    :param lig1: Coordinates of the first ligand
    :param lig2: Coordinates of the second ligand
    :return: Distance between the two ligands
    '''
    return np.linalg.norm(lig1 - lig2)

def get_lig(pdb_file, lig_id):
    '''
    Extracts coordinates of a specific ligand from a PDB file.
    :param pdb_file: Path to the PDB file
    :param lig_id: ID of the ligand to extract
    :return: Numpy array of coordinates for the specified ligand
    '''
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
    '''
    Extracts the ligand ID from a PDB file.
    :param pdb_file: Path to the PDB file
    :return: Ligand ID if found, else None
    '''
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                return line[17:20].strip()
    return None

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
    Extracts the specified ligand from a PDB file and converts it to SMILES format.
    :param pdb_file: Path to the PDB file
    :param lig_id: ID of the ligand to extract
    :param output_name: Name of the output PDB file
    :return: SMILES representation of the ligand if successful, else None
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
    if mol is None:
        print(f"RDKit failed to read {output_name}")
        return None
    return Chem.MolToSmiles(mol)

def convert_with_obabel(input_pdb, output_mol):
    '''
    Converts a PDB file to a MOL file using Open Babel.
    :param input_pdb: Path to the input PDB file
    :param output_mol: Path to the output MOL file
    '''
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

def extract_number(filename):
    '''
    Extracts a number from the filename, specifically for complex PDB files.
    :param filename: Name of the PDB file
    :return: Extracted number or infinity if not found
    '''
    match = re.search(r'complex\.(\d+)\.pdb', filename)
    return int(match.group(1)) if match else float('inf')

def get_main_ligand_id(pdb_file):
    """
    Return the ligand ID with the most atoms, skipping metal ions and solvent.
    :param pdb_file: Path to the PDB file
    :return: Ligand ID with the most atoms or None if not found
    """
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

def find_surface_lysines(pdb_path, E3_ligand_path, asa_threshold=100):

# NOTE: There may be some issues with integrating this onto the cluster, modules may have to downloaded/installed per run
# Is there a better way to do this?

# Add distance calculation

    """
    Identifies surface lysines (residue K) from DSSP output using mkdssp.

    Args:
        pdb_path (str): Path to the cleaned receptor PDB file.
        poi_ligand_path (str): (Unused here, but kept for compatibility)
        asa_threshold (int): Accessible surface area (ASA) threshold to consider a residue 'surface exposed'.

    Returns:
        List of (chain_id, residue_number) tuples for surface lysines.
    """
    pdb_file = Path(pdb_path)
    dssp_file = pdb_file.with_suffix('.dssp')

    # Run mkdssp directly
    subprocess.run([
        "mkdssp", "--output-format", "dssp",
        str(pdb_file), str(dssp_file)
    ], check=True)

    # Parse the DSSP output manually
    surface_lysines = []
    parse_lines = False

    with open(dssp_file) as f:
        for line in f:
            if line.startswith("  #  RESIDUE"):
                parse_lines = True
                continue
            if not parse_lines or len(line) < 40:
                continue
            res_type = line[13]
            try:
                acc = int(line[35:38].strip())
                resnum = int(line[5:10].strip())
                coords = line[11:13].strip()
            except ValueError:
                continue

            if res_type == 'K' and acc >= asa_threshold:
                surface_lysines.append(resnum)

    return surface_lysines

def clean_pdb_for_dssp(input_pdb, output_pdb='cleaned.pdb'):

    with open(input_pdb, 'r') as f:
        lines = f.readlines()

    cleaned_lines = []
    for line in lines:
        # Skip ligand
        if ' LIG ' in line:
            continue
        # Replace HIE with HIS (important for DSSP)
        cleaned_lines.append(line.replace('HIE', 'HIS'))

    with open(output_pdb, 'w') as f:
        f.writelines(cleaned_lines)

    print(f"Cleaned PDB written to: {output_pdb}")

def lys_dist(lysines, pdb_path, lig_coords):
    """
    Calculate distances between surface lysines and a ligand.
    :param lysines: List of surface lysine residue numbers.
    :param pdb_path: Path to the PDB file.
    :param lig_coords: Coordinates of the ligand.
    :return: List of distances.
    """
    distances = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and int(line[22:26].strip()) in lysines:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                lys_coords = np.array([x, y, z])
                dist = distance(lys_coords, lig_coords)
                distances.append(dist)
    return distances

def remove_stereochemistry(smiles):
    '''
    Remove stereochemistry from a SMILES string.
    :param smiles: Input SMILES string
    :return: SMILES string without stereochemistry
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)

def main():
    """
    Loop through all PDB files in the current directory, find ligand pairs,
    compute closest distances, and run LinkInvent using the extracted SMILES.
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

    # Sort and truncate to top 10
    teeny = dict(sorted(teeny.items(), key=lambda item: item[1]))
    teeny = dict(list(teeny.items())[:10])
    lowest_files = sorted(teeny.keys(), key=extract_number)[:10]
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


    clean_pdb_for_dssp(receptor_pdb, 'cleaned_receptor.pdb')
    for file in list(teeny.keys()):
        lig2_coords = get_lig(file, lig2[0])
        if lig2_coords.size == 0:
            continue
        '''
        some_var = sorted(lys_dist(find_surface_lysines('cleaned_receptor.pdb', get_lig(top_file, lig2[0])), file, get_lig(file, lig2[0])))[:5]
        print(f"Surface lysines for {file}: {some_var}")
        with open('lysines.txt', 'w') as f:
            f.write({some_var})
        '''

    if 'smiles.smi' not in os.listdir('.'):
        print("Extracting SMILES for ligands...")
        lig1_smiles = remove_stereochemistry(extract_ligand_smiles(top_file, lig1[0]))
        lig2_smiles = remove_stereochemistry(extract_ligand_smiles(top_file, lig2[0]))

        print("Lig1 SMILES:", lig1_smiles)
        print("Lig2 SMILES:", lig2_smiles)

        lig1_smiles = lig1_smiles + '*'
        lig2_smiles = lig2_smiles + '*'

        with open('smiles.smi', 'w') as f:
            f.write(f"{lig1_smiles}|{lig2_smiles}\n")

    with open('input.txt', 'w') as f:
        f.write(f"{min_val},{max_val}\n")

    keep_files = set(lowest_files)
    for file in os.listdir('.'):
        if file.startswith('complex.') and file.endswith('.pdb') and file not in keep_files:
            try:
                os.remove(file)
                print(f"Deleted unused file: {file}")
            except Exception as e:
                print(f"Failed to delete {file}: {e}")

main()