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
import pandas as pd

def create_param(ligand_pdb, receptor_pdb, warhead1, warhead2, anchor1, anchor2, protac):
    '''
    Create a parameter file for the PROTAC docking.
    :param ligand_pdb: Path to the ligand PDB file.
    :param receptor_pdb: Path to the receptor PDB file.
    :param warhead1: SMILES string for the first warhead.
    :param warhead2: SMILES string for the second warhead.
    :param anchor_atoms: Anchor atoms for the linker attachment.
    :param smiles: Full SMILES string for the PROTAC.
    '''

    with open('parameters.txt', 'w') as f:
        f.write(f'''
        For main.py / extended.py:
        Structures: {ligand_pdb} {receptor_pdb}
        Chains: A B  
        Heads: headA.sdf headB.sdf 
        Anchor atoms: {anchor1} {anchor2}
        Protac: {protac}
        Full: True
        ''')

    return 'parameters.txt'

def get_ligand_sdf(smiles, name):
    '''
    Convert a SMILES string to an SDF string.
    :param smiles: Input SMILES string
    :return: SDF string
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    with open(f'{name}.sdf', 'w') as f:
        f.write(Chem.MolToMolBlock(mol))
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
    ligands = smiles.split('|')
    anchor_atoms = []

    for smi in ligands:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smi}")
        
        anchor_idx = None
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                anchor_idx = atom.GetIdx()
                break  # Take first anchor atom found
        
        if anchor_idx is None:
            raise ValueError(f"No anchor atom (atom map) found in ligand: {smi}")
        anchor_atoms.append(anchor_idx)

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

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                residue = line[17:20].strip()
                if residue in standard_residues:
                    output_lines.append(line)

    output_file = f"{pdb_file.replace('.pdb', '')}_nolig.pdb"
    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    print(f"Ligand-removed PDB written to: {output_file}")
    return output_file

def main():
    with open ("smiles.smi", "r") as f:
        smiles = f.read().strip()
    protac = get_ligand_sdf((get_PROTAC('linkinvent_stage_1.csv', output_path='top_smiles.txt', top_n=10)[0]), 'protac')
    mol, distance = minimize_and_measure(protac, get_anchor_atoms(smiles), 'min_protac')

    with open('input.txt', 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

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
    
    for p in proteins:
        ternary = remove_ligand(p)

        config = f'''receptor = {ternary}
                    ligand = protac.sdf
                    autobox_ligand = {ternary}
                    out = docked.sdf.gz
                    log = log
                    cnn_scoring = rescore
                    num_modes = 25
                    exhaustiveness = 32
                    pose_sort_order = Energy
                    '''
        try:
            with open('config', 'w') as config_file:
                config_file.write(config)
        except IOError:
            print("Error: Could not write config file.")
            exit(1)

        job_script = f'''#!/bin/bash
    #SBTACH --job-name={ternary}
    #SBATCH --cpus-per-task=3
    #SBATCH --mem-per-cpu=4000M
    #SBATCH --gres=gpu:1
    #SBATCH --time=1:00:00
    #SBATCH --account=def-aminpour
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=jaharri1@ualberta.ca

    module load gnina/1.3.1
    module load StdEnv/2023
    module load python/3.11
    module load scipy-stack/2025a
    module load rdkit/2024.09.6
    module load openbabel/3.1.1
    module load gcc/12.3
    module load cmake
    module load cuda/12.2
    module load python-build-bundle/2025b

    gnina --config config
    '''
        
        try:
            with open('job.sh', 'w') as job_script_file:
                job_script_file.write(job_script)
        except IOError:
            print("Error: Could not write job script.")
            exit(1)

        # Submit job
        try:
            os.system('sbatch job.sh')
        except Exception as e:
            print(f"Error: Could not submit job: {e}")
            exit(1)

main()