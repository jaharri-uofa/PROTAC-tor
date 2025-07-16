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
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

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

def find_surface_lysines(POI_pdb_file, E3_lig_coords):

# NOTE: There may be some issues with integrating this onto the cluster, modules may have to downloaded/installed per run
# Is there a better way to do this?

    '''
    Takes a PDB file and returns a list of surface lysines that are acceasible for Ubiquitination by the E3 ligase.
    1.) Find all surface lysines in the PDB file using RSA (Shrake-Rupley algorithm).
    2.) Return a list of these lysines and their coordinates and distance to the E3 Ligase warhead
    :param POI_pdb_file: Path to the PDB file of the protein of interest
    :param E3_lig_coords: Coordinates of the E3 ligase warhead
    :return: List of tuples containing residue ID, coordinates, and distance to E3 ligase warhead
    '''
    rsa_threshold=0.3
    surface_lysines = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("POI", POI_pdb_file)

    # Call DSSP (make sure mkdssp is in PATH)
    dssp = DSSP(structure[0], POI_pdb_file)

    for key in dssp.keys():
        residue = dssp[key][0]  # Biopython residue object
        resname = residue.get_resname()
        rsa = dssp[key][3]  # Relative Solvent Accessibility

        if resname == "LYS" and rsa > rsa_threshold:
            if "CA" in residue:
                coords = residue["CA"].get_coord()
                dist = distance(coords, E3_lig_coords)
                resid = residue.get_id()[1]
                surface_lysines.append((resid, coords, dist))

    return surface_lysines

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
    
    print(f'Surface lysines found:', find_surface_lysines(receptor_pdb, get_lig(top_file, lig2[0])))

    if 'smiles.csv' not in os.listdir('.'):
        print("Extracting SMILES for ligands...")
        lig1_smiles = extract_ligand_smiles(top_file, lig1[0])
        lig2_smiles = extract_ligand_smiles(top_file, lig2[0])

        print("Lig1 SMILES:", lig1_smiles)
        print("Lig2 SMILES:", lig2_smiles)

        lig1_smiles = lig1_smiles + '*'
        lig2_smiles = lig2_smiles + '*'

        with open('smiles.csv', 'w') as f:
            f.write("fragment_1,fragment_2\n")
            f.write(f"{lig1_smiles},{lig2_smiles}\n")

    with open('input.txt', 'w') as f:
        f.write(f"{min_val},{max_val}\n")

main()