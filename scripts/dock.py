'''
This script takes two pdb files, the POI and an E3 ligase, and the PROTAC SMILES string to run PRosettaC,
a PROTAC docking software. It creates a directory for the complex based off the PROTAC number and saves the 
best docking results. 
'''

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

def get_ligand_sdf(smiles_list, name):
    '''
    Convert a list of SMILES strings to a single SDF file containing all molecules.
    :param smiles_list: List of SMILES strings
    :param name: Output SDF file base name (no extension)
    :return: SDF filename
    '''
    writer = Chem.SDWriter(f'{name}.sdf')
    for smile in smiles_list:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            print(f"Warning: Invalid SMILES string skipped: {smile}")
            continue
        mol = Chem.AddHs(mol)  # Add hydrogens for better 3D
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(mol)
        writer.write(mol)
    writer.close()
    return f'{name}.sdf'

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

def pdb_to_sdf(input_pdb, output_sdf):
    '''
    Takes a pdb as input and converts it to an .sdf file
    :param pdb: pdb file
    :return: sdf file of the pdb
    '''
    try:
        subprocess.run([
            "obabel", input_pdb,
            "-O", output_sdf,
            "-h", "--gen3d"
        ], check=True)
        print(f"Converted {input_pdb} → {output_sdf}")
    except subprocess.CalledProcessError as e:
        print(f"Conversion failed: {e}")
    


def get_PROTAC(csv_file, output_path='top_smiles.txt', top_n=10):
    '''
    Get the PROTAC information from a CSV file.
    :param csv_file: Path to the CSV file
    '''
    df = pd.read_csv(csv_file)
    sorted_df = df.sort_values(by='Score', ascending=False)

    top_smiles = sorted_df.head(top_n)['SMILES'].tolist()
    with open(output_path, 'w') as f:
        for smiles in top_smiles:
            f.write(smiles + '\n')
        print(f"SMILES strings saved to {output_path}")
    return top_smiles

def minimize_and_measure(sdf_file, anchor_atom_indices, output_sdf):
    """
    Parameters:
        sdf_file (str): Path to the input SDF file (assumes one molecule).
        anchor_atom_indices (tuple/list of 2 ints): Atom indices of the two anchors (0-based).
        
    Returns:
        mol_minimized (rdkit.Chem.Mol): Molecule with minimized conformation.
        distance (float): Distance between the two anchor atoms in Angstroms after minimization.
    """
    # Load molecule
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol = suppl[0]
    if mol is None:
        raise ValueError("Could not load molecule from SDF.")
    
    # Make sure molecule has 3D coords or generate them
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    
    # Minimize the molecule (UFF or MMFF)
    if AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol)
    else:
        AllChem.MMFFOptimizeMolecule(mol)
    
    # Get minimized conformer
    conf = mol.GetConformer()
    
    # Measure distance between anchor atoms
    idx1, idx2 = anchor_atom_indices
    pos1 = conf.GetAtomPosition(idx1)
    pos2 = conf.GetAtomPosition(idx2)
    distance = pos1.Distance(pos2)

    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()
    
    return mol, distance

def get_anchor_atoms(smiles):
    """
    Finds anchor atom indices in each ligand, and modifies the SMILES to add '[' before '*' and ';]' after '*'.
    Returns anchor atom indices as a tuple.
    """
    ligands = smiles.split('|')
    anchor_atoms = []
    modified_ligands = []

    for smi in ligands:
        # Add [ before '*' and ;] after '*'
        if '*' in smi:
            smi_mod = smi.replace('*', '[*;]')
        else:
            smi_mod = smi
        modified_ligands.append(smi_mod)

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi}")

        anchor_idx = None
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == '*' or atom.GetAtomMapNum() > 0:
                anchor_idx = atom.GetIdx()
                break  # Take first anchor atom found

        if anchor_idx is None:
            raise ValueError(f"No anchor atom (atom map or '*') found in ligand: {smi}")
        anchor_atoms.append(anchor_idx)

    # Optionally, you can return the modified ligands as well
    # return anchor_atoms[0], anchor_atoms[1], '|'.join(modified_ligands)
    return anchor_atoms[0], anchor_atoms[1]

def remove_ligand(pdb_file):
    '''
    Takes a PDB file as input, identifies non-standard residues (e.g., ligands),
    and writes a new PDB file with those residues removed.

    Assumes ligand residues are not named like standard amino acids (e.g., 'LIG').
    '''

    # Standard residue names (3-letter codes for 20 AAs + common additions)
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
        'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO',
        'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HOH', 'WAT', 'H2O'  # common water residue names
    }

    output_lines = []
    lig_lines = []

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                residue = line[17:20].strip()
                if residue in standard_residues:
                    output_lines.append(line)
                elif residue not in standard_residues:
                    lig_lines.append(line)

    output_file = f"{pdb_file.replace('.pdb', '')}_nolig.pdb"
    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    lig_output_file = f"{pdb_file.replace('.pdb', '')}_lig.pdb"
    with open(lig_output_file, 'w') as f:
        f.writelines(lig_lines)

    print(f"Ligand-removed PDB written to: {output_file}")
    print(f"Ligand-only PDB written to: {lig_output_file}")
    return output_file, lig_output_file

def main():
    with open("smiles.smi", "r") as f:
        smiles = f.read().strip()
    protac_smiles_list = get_PROTAC('linkinvent_stage_1.csv', output_path='top_smiles.txt', top_n=10)
    sdf_file = get_ligand_sdf(protac_smiles_list, 'protac')

    anchor_indices = get_anchor_atoms(smiles)
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)

    count = 0
    for mol in suppl:
        if mol is None:
            continue
        # Save each molecule to a temporary SDF for minimization
        temp_sdf = f"min_protac_{count}.sdf"
        writer = Chem.SDWriter(temp_sdf)
        writer.write(mol)
        writer.close()

        mol_min, distance = minimize_and_measure(temp_sdf, anchor_indices, temp_sdf)

        with open('input.txt', 'r') as f:
            min_dist, max_dist = map(float, f.readline().strip().split(','))

        print(f'Max Complex Distance: {max_dist}, Min Complex Distance: {min_dist}, PROTAC minimized distance: {distance}')

        if distance <= max_dist and distance >= min_dist:
            print('valid protac conformation acheived, proceeding to docking')
        elif distance > max_dist:
            print('protac conformation is too large! adjust linkinvent model')
        elif distance < min_dist:
            print('protac conformation is too small! adjust linkinvent model')

        proteins = []
            
        with open('lig_distances.txt', 'r') as f:
            for line in f:
                if ':' in line:
                    filename = line.split(':')[0].strip()
                    proteins.append(filename)

        for p in proteins[:1]:
            ternary, ligands = remove_ligand(p)
            ligand = pdb_to_sdf(ligands, 'ligand.sdf')
            job_dir = f"docking_{os.path.splitext(os.path.basename(ternary))[0]}_{count}"
            os.makedirs(job_dir, exist_ok=True)

            # Write config file in job_dir
            config = f'''receptor = {ternary}
ligand = protac.sdf
autobox_ligand = {ternary}
out = docked.sdf.gz
log = log
cnn_scoring = none
num_modes = 25
exhaustiveness = 32
pose_sort_order = Energy
                        '''
            try:
                with open(os.path.join(job_dir, 'config'), 'w') as config_file:
                    config_file.write(config)
            except IOError:
                print(f"Error: Could not write config file in {job_dir}.")
                exit(1)

            # Write job script in job_dir
            job_script = f'''#!/bin/bash
#SBATCH --job-name={ternary}
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=256M
##SBATCH --gres=gpu:1
#SBATCH --time=1:00:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module load StdEnv/2023
module load python/3.11
module load scipy-stack/2025a
module load rdkit/2024.09.6
module load openbabel/3.1.1
module load gcc/12.3
module load cmake
module load cuda/12.2
module load python-build-bundle/2025b
module load gnina/1.3.1

gnina --config config
        '''
            try:
                with open(os.path.join(job_dir, 'job.sh'), 'w') as job_script_file:
                    job_script_file.write(job_script)
            except IOError:
                print(f"Error: Could not write job script in {job_dir}.")
                exit(1)

            # Copy protac.sdf to job_dir if needed
            if not os.path.exists(os.path.join(job_dir, 'protac.sdf')):
                try:
                    shutil.copy('protac.sdf', os.path.join(job_dir, 'protac.sdf'))
                except Exception as e:
                    print(f"Error copying protac.sdf to {job_dir}: {e}")
                    exit(1)

            # Copy stripped protein complex to the directory
            try:
                shutil.copy(ternary, os.path.join(job_dir, os.path.basename(ternary)))
            except Exception as e:
                print(f"Error copying {ternary} to {job_dir}: {e}")
                exit(1)

            # Submit job from job_dir
            try:
                os.system(f'cd {job_dir} && sbatch job.sh')
            except Exception as e:
                print(f"Error: Could not submit job in {job_dir}: {e}")
                exit(1)
        
        count = count + 1
    
    

main()