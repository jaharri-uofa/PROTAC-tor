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

def delta_G(pdb_file):
    '''
    Calculates the free energy of a complex
    :param pdb_file: a pdb file with a ligand(s) and two docked proteins
    :return: the free energy of the complex
    '''

    sdf_tmp = pdb_file.replace('.pdb', '_tmp.sdf')
    os.system(f'obabel {pdb_file} -O {sdf_tmp} --gen3d')

    mol = Chem.SDMolSupplier(sdf_tmp, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"RDKit failed to load molecule from {sdf_tmp}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)

    if not AllChem.UFFHasAllMoleculeParams(mol):
        raise ValueError("Molecule has unsupported atoms for UFF.")

    ff = AllChem.UFFGetMoleculeForceField(mol)
    energy = ff.CalcEnergy()
    print(f"Estimated total system energy: {energy:.2f} kcal/mol")
    return energy


def main():
    proteins = []
    base = {}
            
    with open('lig_distances.txt', 'r') as f:
        for line in f:
            if ':' in line:
                filename = line.split(':')[0].strip()
                proteins.append(filename)
    
    for p in proteins:
        base[p] = delta_G(p)

    for dir in os.listdir():
        if os.path.isdir(dir) and dir.startswith("docking_"):
            docked_sdf = os.path.join(dir, "docked.sdf.gz")
            ternary_pdb = None
            for file in os.listdir(dir):
                if file.endswith("_nolig.pdb"):
                    ternary_pdb = os.path.join(dir, file)
                    break
            if ternary_pdb and os.path.exists(docked_sdf):
                # Unzip docked.sdf.gz if needed
                sdf_out = os.path.join(dir, "docked.sdf")
                if not os.path.exists(sdf_out):
                    subprocess.run(["gunzip", "-c", docked_sdf], stdout=open(sdf_out, "wb"))
                # Combine protein and ligand
                complex_pdb = add_ligand(ternary_pdb, sdf_out)
                # Calculate free energy
                try:
                    energy = delta_G(complex_pdb)
                    print(f"Free energy for {complex_pdb}: {energy:.2f} kcal/mol")
                    # Find the matching base complex (by name)
                    base_name = os.path.basename(ternary_pdb).replace("_nolig.pdb", ".pdb")
                    base_energy = base.get(base_name)
                    if base_energy is not None and energy < base_energy:
                        print(f"Deleting {dir} (energy {energy:.2f} < base {base_energy:.2f})")
                        shutil.rmtree(dir)
                except Exception as e:
                    print(f"Failed to calculate free energy for {complex_pdb}: {e}")



main()