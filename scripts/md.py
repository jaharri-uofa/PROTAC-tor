#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
import stat
import subprocess
import logging as log
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

def add_ligand(pdb_file, sdf):
    '''
    Within the docking complex this will take the sdf file and pdb and combine the two into a single pdb
    file for MD, as well as saving off the ligand.pdb and receptor.pdb
    :param pdb_file: pdb file of docked proteins
    :param sdf: an sdf file of the compound
    :return: a pdb with the protac docked onto it
    '''
    lig_pdb = sdf.replace('.sdf', '_lig.pdb')
    os.system(f'obabel {sdf} -O {lig_pdb}')

    with open(pdb_file, 'r') as f:
        protein_lines = f.readlines()
    with open(lig_pdb, 'r') as f:
        ligand_lines = f.readlines()

    combined = pdb_file.replace('.pdb', '_complex.pdb')
    receptor = pdb_file.replace('.pdb', '_receptor.pdb')
    ligand = pdb_file.replace('.pdb', '_ligand.pdb')
    with open(combined, 'w') as f:
        for line in protein_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                f.write(line)
        for line in ligand_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                f.write(line)
        f.write('END\n')
    with open(receptor, 'w') as f:
        for line in protein_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                f.write(line)
        f.write('END\n')
    with open(combined, 'w') as f:
        for line in ligand_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                f.write(line)
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
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        f_out.write(f_in.read())

def sdf_to_smiles_affinity(sdf_path):
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

def main():
    print("Starting main process...")
    docking_dirs = [d for d in os.listdir() if os.path.isdir(d) and d.startswith('dock')]
    print(f"Found docking directories: {docking_dirs}")
    all_complexes = []
    print(os.listdir())
    for dock_dir in docking_dirs:
        print(f"\nProcessing docking directory: {dock_dir}") 
        gz_file = os.path.join(dock_dir, 'docked.sdf.gz')
        sdf_file = os.path.join(dock_dir, 'docked.sdf')
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

    # Sort by affinity (lower is better, adjust if higher is better)
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
        for row in reader:
            smiles = row['smiles']
            affinity = row['affinity']
            sdf_path = row['sdf_path']
            pdb_path = row['dock_dir']
            pdb_path = str(pdb_path) + str(pdb_path).split('docking_')[-1]
            print(f"SMILES: {smiles}, Affinity: {affinity}, SDF: {sdf_path}, PDB: {pdb_path}")
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
            
            add_ligand('ternary.pdb', 'ligand.sdf')
            subprocess.run(['python', 'md_mmgbsa.py', f'ternary_complex.pdb', f'ternary_receptor.pdb', f'ternary_ligand.pdb'])
                        

main()