#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
import stat
import subprocess
import logging as log
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolfiles
import pandas as pd
import re
import numpy as np
import gzip
import csv

standard_residues = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'SEC', 'PYL', 'HIE', 'HIP', 'ASH',
    'GLH', 'CYM', 'CYX', 'LYN', 'ACE',
    'NME','HOH', 'WAT', 'H2O'  # common water residue names
}

skip_residues = {
    'NA', 'CL', 'CA', 'MG', 'ZN', 'K', 'FE', 'CU', 'MN', 'HG',
    'HOH', 'WAT', 'SO4', 'PO4', 'HEM', 'DMS', 'ACE', 'NAG', 'GLC'
}

def create_receptor_ligand_files(pdb_file):
    '''
    Takes a PDB file as input, identifies non-standard residues (e.g., ligands),
    and writes a new PDB file with those residues removed.

    Assumes ligand residues are not named like standard amino acids (e.g., 'LIG').
    '''

    # Standard residue names (3-letter codes for 20 AAs + common additions/variants)

    output_lines = []
    lig_lines = []
    print('starting to process PDB file')
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                residue = line[17:20].strip()
                if residue in standard_residues:
                    output_lines.append(line)
                elif residue not in standard_residues:
                    lig_lines.append(line)

    print('Writing receptor PDB file and ligand PDB file')
    output_file = f"receptor.pdb"
    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    lig_output_file = f"ligand.pdb"
    with open(lig_output_file, 'w') as f:
        f.writelines(lig_lines)

    print(f"Ligand-removed PDB written to: {output_file}")
    print(f"Ligand-only PDB written to: {lig_output_file}")
    return output_file, lig_output_file

def add_ligand(receptor_pdb_file, ligand_sdf_file, count):
    """
    Combines receptor PDB and ligand SDF into:
    - complex PDB (protein + ligand)
    - receptor PDB (protein only)
    - ligand PDB (ligand only)
    :param receptor_pdb_file: Path to the receptor PDB file
    :param ligand_sdf_file: Path to the ligand SDF file
    :param count: Counter for output file naming
    """
    # Read receptor PDB lines
    with open(receptor_pdb_file, 'r') as f:
        receptor_lines = [l for l in f if l.startswith(('ATOM', 'HETATM'))]

    # Convert ligand SDF to PDB using RDKit
    suppl = Chem.SDMolSupplier(ligand_sdf_file, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(f"Could not read ligand from {ligand_sdf_file}")

    ligand_pdb_block = Chem.MolToPDBBlock(mol)
    ligand_lines = [l for l in ligand_pdb_block.splitlines() if l.startswith(('ATOM', 'HETATM'))]

    # Write files
    base = receptor_pdb_file.replace('.pdb', '')
    combined = f"{base}_complex{count}.pdb"
    receptor = f"{base}_receptor{count}.pdb"
    ligand = f"{base}_ligand{count}.pdb"

    with open(combined, 'w') as f:
        for line in receptor_lines:
            f.write(line)
        for line in ligand_lines:
            f.write(line)
        f.write('END\n')
    with open(receptor, 'w') as f:
        for line in receptor_lines:
            f.write(line)
        f.write('END\n')
    with open(ligand, 'w') as f:
        for line in ligand_lines:
            f.write(line + '\n')
        f.write('END\n')

    print(f"combined pdb written to {combined}")
    print(f"receptor pdb written to {receptor}")
    print(f"ligand pdb written to {ligand}")
    return combined, receptor, ligand

def distance(lig1, lig2):
    '''
    Calculates distance between two ligands
    :param lig1: Coordinates of the first ligand
    :param lig2: Coordinates of the second ligand
    :return: Distance between the two ligands
    '''
    return np.linalg.norm(lig1 - lig2)

# why is this here?
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

def extract_sdf_gz(gz_path, out_path):
    '''
    unzips a file
    :param gz_path: Path to the input gzipped SDF file.
    :param out_path: Path to the output SDF file.
    '''
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        f_out.write(f_in.read())

def sdf_to_smiles_affinity(sdf_path):
    '''
    Converts SDF to SMILES and extracts affinity information.
    :param sdf_path: Path to the input SDF file.
    '''
    suppl = Chem.SDMolSupplier(sdf_path)
    results = []
    mol_count = 0

    for mol in suppl:
        if mol is None:
            continue

        mol_count += 1
        smiles = Chem.MolToSmiles(mol)
        props = mol.GetPropNames()
        print(f"  Molecule {mol_count}: properties = {list(props)}")
        affinity = None
        for key in ['minimizedAffinity']:
            if mol.HasProp(key):
                affinity = float(mol.GetProp(key))
                break
        if affinity is not None:
            results.append({'smiles': smiles, 'affinity': affinity, 'sdf_path': sdf_path})
    print(f"  Total molecules read: {mol_count}")
    return results

def get_main_ligand_id(pdb_file):
    """
    Return the ligand ID with the most atoms, skipping metal ions and solvent.
    :param pdb_file: Path to the PDB file
    :return: Ligand ID with the most atoms or None if not found
    """
    residue_atom_counts = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname not in skip_residues:
                    residue_atom_counts[resname] = residue_atom_counts.get(resname, 0) + 1
    return max(residue_atom_counts, key=residue_atom_counts.get) if residue_atom_counts else None

def get_control_ligand_id(pdb_file):
    '''
    Gets the ligands from the control file and returns them
    Note: This function is explicitly for proteins with multiple ligands and no HETATM tags
    '''
    previous_line = None
    ligands = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname not in standard_residues and resname not in skip_residues and line != previous_line:
                    ligands.append(resname)
            previous_line = line
    return ligands

def split_control_pdb(pdb):
    '''
    Takes the control pdb and splits it into the receptor file, and individual ligand files
    Note: The complex.pdb file doesnt have HETATM headers
    '''
    with open(pdb, 'r') as f:
        lines = f.readlines()

    # Write the receptor file
    with open('receptor.pdb', 'w') as f:
        for line in lines:
            resname = line[17:20].strip()
            if line.startswith("ATOM") and resname in standard_residues:
                f.write(line)

    # Write individual ligand files
    ligands = get_control_ligand_id(pdb)
    for lig in ligands:
        with open(f'ligand_{lig}.pdb', 'w') as f:
            for line in lines:
                if line.startswith("ATOM"):
                    resname = line[17:20].strip()
                    if resname == lig:
                        f.write(line)

def main():
    print("Starting main process...")

    # collect docking directories
    docking_dirs = [d for d in os.listdir() if os.path.isdir(d) and d.startswith('dock')]
    print(f"Found docking directories: {docking_dirs}")
    all_complexes = []
    print(os.listdir())
    for i, dock_dir in enumerate(docking_dirs):
        print(f"\nProcessing docking directory: {dock_dir}") 
        gz_file = os.path.join(dock_dir, 'docked.sdf.gz')
        sdf_file = os.path.join(dock_dir, 'docked.sdf')

        # unpack files
        if not os.path.exists(gz_file):
            print(f"  ERROR: {gz_file} does not exist. Skipping {dock_dir}.")
            continue
        try:
            extract_sdf_gz(gz_file, sdf_file)
            print(f"  Extracted {gz_file} to {sdf_file}")
        except Exception as e:
            print(f"  ERROR extracting {gz_file}: {e}")
            continue
        try:
            complexes = sdf_to_smiles_affinity(sdf_file)
            print(f"  Found {len(complexes)} complexes in {sdf_file}")
        except Exception as e:
            print(f"  ERROR parsing {sdf_file}: {e}")
            continue
        for c in complexes:
            c['dock_dir'] = dock_dir
        all_complexes.extend(complexes)

    if not all_complexes:
        print("No complexes found in any docking directory. Exiting.")
        return

    try:
        all_complexes.sort(key=lambda x: x['affinity'])
        top5 = all_complexes[:5]
        top_dirs = set(c['dock_dir'] for c in top5)
        print(f"\nTop 5 complexes (by affinity):")
        for i, c in enumerate(top5, 1):
            print(f"  {i}. Dir: {c['dock_dir']}, Affinity: {c['affinity']}, SMILES: {c['smiles']}")
    except Exception as e:
        print(f"ERROR sorting/selecting top complexes: {e}")
        return

    # Remove duplicate entries with the same SMILES, keeping the one with the best (lowest) affinity
    best_by_smiles = {}
    for entry in all_complexes:
        smiles = entry['smiles']
        affinity = entry['affinity']
        if smiles not in best_by_smiles or affinity < best_by_smiles[smiles]['affinity']:
            best_by_smiles[smiles] = entry
    unique = list(best_by_smiles.values())
    unique.sort(key=lambda x: x['affinity'])  # sort by best affinity

    print(f"\nFiltered to {len(unique)} unique complexes (best affinity per SMILES).")

    # Save unique top complexes as CSV (optionally limit to top N)
    try:
        pd.DataFrame(unique[:5]).to_csv('top5_complexes.csv', index=False)
        print("\nTop unique complexes saved to top5_complexes.csv")
    except Exception as e:
        print(f"ERROR saving CSV: {e}")

    with open('top5_complexes.csv', 'r') as f:
        reader = csv.DictReader(f)
        pdb_paths = []
        for i, row in enumerate(reader):
            smiles = row['smiles']
            affinity = row['affinity']
            sdf_path = row['sdf_path']
            pdb_path = row['dock_dir']
            pdb_path = str(pdb_path) + '/' + str(pdb_path).split('docking_')[-1].split('g_')[0] + 'g.pdb'
            print(f"SMILES: {smiles}, Affinity: {affinity}, SDF: {sdf_path}, PDB: {pdb_path}")
            pdb_paths.append(pdb_path)
            # need to find the corresponding ligand in protac.sdf
            with open(sdf_path, 'r') as sdf_file:
                molecules = sdf_file.read().split("$$$$\n")

            found = False
            for mol in molecules:
                if f"<minimizedAffinity>\n{affinity}" in mol:
                    out_name = f"ligand.sdf"
                    with open(out_name, 'w') as out:
                        out.write(mol.strip() + "\n$$$$\n")
                    print(f"Extracted ligand for {smiles} to {out_name}")
                    found = True
                    break

            with open(pdb_path, 'r') as pdb_file:
                pdb_content = pdb_file.read()
                pdb_out_name = f'ternary.pdb'
                with open(pdb_out_name, 'w') as pdb_out:
                    pdb_out.write(pdb_content)

            if not found:
                print(f"Affinity {affinity} not found in {sdf_path}")

            outdir = f"ternary_complex{i+1}"
            os.makedirs(outdir, exist_ok=True)

            combined, receptor, ligand = add_ligand('ternary.pdb', 'ligand.sdf', i+1)

            ligand_resname = get_main_ligand_id(combined)
            if ligand_resname is None:
                print("ERROR: No ligand found in ternary.pdb")
                sys.exit(1)
            with open('ligand_resname.txt', 'w') as f:
                f.write(ligand_resname)

            shutil.move(combined, os.path.join(outdir, os.path.basename(combined)))
            shutil.move(receptor, os.path.join(outdir, os.path.basename(receptor)))
            shutil.move(ligand, os.path.join(outdir, os.path.basename(ligand)))
            shutil.move('ligand_resname.txt', os.path.join(outdir, 'ligand_resname.txt'))
    
            shutil.move('ligand.sdf', os.path.join(outdir, 'ligand.sdf'))
            shutil.move('ternary.pdb', os.path.join(outdir, 'ternary.pdb'))

            subprocess.run(
                ['python', '../md_mmgbsa.py', os.path.basename(combined), os.path.basename(receptor), os.path.basename(ligand)],
                cwd=outdir
                )
    # this is broken :)
    
    freq = {}
    print
    for pdb in pdb_paths:
        # Extract numeric ID from filename using regex
        match = re.search(r'(\d+)', os.path.basename(pdb))
        if match:
            num = match.group(1)
            freq[num] = freq.get(num, 0) + 1

    control_dir = f"ternary_complex_control"
    os.makedirs(control_dir, exist_ok=True)

    most_common = max(freq, key=freq.get)
    print(f"Most common PDB ID: {most_common}")

    # Look for corresponding PDB file in current directory
    candidate = f"complex.{most_common}.pdb"
    if os.path.exists(candidate):
        print(f"Found corresponding PDB: {candidate}")
    else:
        print(f"Could not find {candidate} in directory.")

    # test case, delete this later
    candidate = 'complex.165.pdb'

    ligand_resname = get_control_ligand_id(candidate)
    print(ligand_resname)
    if ligand_resname is None:
        print("ERROR: No (s) found in ternary.pdb")
        sys.exit(1)
    with open('ligand_resname.txt', 'w') as f:
        for lig in ligand_resname:
            f.write(lig + "\n")

    split_control_pdb(candidate)
    ligands = get_control_ligand_id(candidate)
    print(ligands)

    shutil.copy(candidate, os.path.join(control_dir, 'complex.pdb'))
    shutil.move('receptor.pdb', os.path.join(control_dir, 'receptor.pdb'))
    shutil.move(f'ligand_{ligands[0]}.pdb', os.path.join(control_dir, 'ligand1.pdb'))
    shutil.move(f'ligand_{ligands[1]}.pdb', os.path.join(control_dir, 'ligand2.pdb'))
    shutil.move('ligand_resname.txt', os.path.join(control_dir, 'ligand_resname.txt'))

    cwd=os.path.join(os.getcwd(), control_dir)

    subprocess.run(
        ['python', '../control_md.py', os.path.basename('complex.pdb'), os.path.basename('receptor.pdb'), os.path.basename('ligand1.pdb'), os.path.basename('ligand2.pdb')],
        cwd=cwd
    )

main()