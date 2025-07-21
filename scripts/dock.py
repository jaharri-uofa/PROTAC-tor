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

    return 'params.txt'

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


def main():
    pass

main()