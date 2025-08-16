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

def add_ligand(receptor_pdb_file, ligand_sdf_file, count):
    """
    Combines receptor PDB and ligand SDF into:
    - complex PDB (protein + ligand)
    - receptor PDB (protein only)
    - ligand PDB (ligand only)
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
            f.write(line + '\n')
        for line in ligand_lines:
            f.write(line + '\n')
        f.write('END\n')
    with open(receptor, 'w') as f:
        for line in receptor_lines:
            f.write(line + '\n')
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

def main():
    print("Starting main process...")
    docking_dirs = [d for d in os.listdir() if os.path.isdir(d) and d.startswith('dock')]
    print(f"Found docking directories: {docking_dirs}")
    all_complexes = []
    print(os.listdir())
    for i, dock_dir in enumerate(docking_dirs):
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
        for i, row in enumerate(reader):
            smiles = row['smiles']
            affinity = row['affinity']
            sdf_path = row['sdf_path']
            pdb_path = row['dock_dir']
            pdb_path = str(pdb_path) + '/' + str(pdb_path).split('docking_')[-1].split('g_')[0] + 'g.pdb'
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

            # Create a unique directory for each result
            outdir = f"ternary_complex{i+1}"
            os.makedirs(outdir, exist_ok=True)

            combined, receptor, ligand = add_ligand('ternary.pdb', 'ligand.sdf', i+1)

            # Save ligand_resname for use in md_mmgbsa.py (e.g., write to a file)
            ligand_resname = get_main_ligand_id(combined)
            if ligand_resname is None:
                print("ERROR: No ligand found in ternary.pdb")
                sys.exit(1)
            # Write ligand_resname.txt directly into outdir
            with open(('ligand_resname.txt'), 'w') as f:
                f.write(ligand_resname)

            # Move the generated files into the output directory
            shutil.move(combined, os.path.join(outdir, os.path.basename(combined)))
            shutil.move(receptor, os.path.join(outdir, os.path.basename(receptor)))
            shutil.move(ligand, os.path.join(outdir, os.path.basename(ligand)))
            shutil.move('ligand_resname.txt', os.path.join(outdir, 'ligand_resname.txt'))
    
            # Also move the ligand.sdf and ternary.pdb for reference
            shutil.move('ligand.sdf', os.path.join(outdir, 'ligand.sdf'))
            shutil.move('ternary.pdb', os.path.join(outdir, 'ternary.pdb'))

            # Run md_mmgbsa.py in the output directory by calling from parent dir
            subprocess.run(
                ['python', '../md_mmgbsa.py', os.path.basename(combined), os.path.basename(receptor), os.path.basename(ligand)],
                cwd=outdir
)
main()