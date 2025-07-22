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

def create_param(ligand_pdb, receptor_pdb, warhead1, warhead2, anchor_atoms, smiles):
    '''
    Create a parameter file for the PROTAC docking.
    :param ligand_pdb: Path to the ligand PDB file.
    :param receptor_pdb: Path to the receptor PDB file.
    :param warhead1: SMILES string for the first warhead.
    :param warhead2: SMILES string for the second warhead.
    :param anchor_atoms: Anchor atoms for the linker attachment.
    :param smiles: Full SMILES string for the PROTAC.
    '''

    with open('params.txt', 'w') as f:
        f.write(f'''
        /////////////////////////////
        For main.py / extended.py:
        Structures: {ligand_pdb} {receptor_pdb}
        Chains: A B  # need to figure out how to get chains? Maybe?
        Heads: {warhead1} {warhead2} 
        Anchor atoms: {anchor_atoms}
        Protac: {smiles}
        Full: True
        ////////////////////////////
        ''')

    return 'parameters.txt'

def get_ligand_sdf(smiles):
    '''
    Convert a SMILES string to an SDF string.
    :param smiles: Input SMILES string
    :return: SDF string
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return Chem.MolToMolBlock(mol)

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

def get_anchor_atoms(smiles):
    '''
    Extract the anchor atoms from a SMILES string.
    Anchor atoms is defined as the atoms that are part of the warhead.
    :param smiles: Input SMILES string
    :return: number of atoms in the warhead
    '''
    warhead1, warhead2 = extract_warhead_smiles(smiles)
    return Chem.MolFromSmiles(warhead1).GetNumAtoms(), Chem.MolFromSmiles(warhead2).GetNumAtoms()

def main():
    ligand = "ligand.pdb" # Path to the E3 ligase PDB file
    receptor = "receptor.pdb" # Path to the POI PDB file
    with open ("smiles.smi", "r") as f:
        smiles = f.read().strip()
    warhead1, warhead2 = extract_warhead_smiles(smiles)
    anchor1, anchor2 = get_anchor_atoms(smiles)
    warhead1 = get_ligand_sdf(warhead1)
    warhead2 = get_ligand_sdf(warhead2)
    param_file = create_param(ligand, receptor, warhead1, warhead2, (anchor1, anchor2), smiles)
    print(f"Parameter file created: {param_file}")

main()