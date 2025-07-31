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

def add_ligand(pdb_file, sdf):
    '''
    Within the docking complex this will take the sdf file and pdb and combine the two into a single pdb
    file for MD
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
    with open(combined, 'w') as f:
        for line in protein_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                f.write(line)
        for line in ligand_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                f.write(line)
        f.write('END\n')

    print(f"combined pdb written to {combined}")
    return combined


def main():
    proteins = []
    base = {}
            
    with open('lig_distances.txt', 'r') as f:
        for line in f:
            if ':' in line:
                filename = line.split(':')[0].strip()
                proteins.append(filename)




main()