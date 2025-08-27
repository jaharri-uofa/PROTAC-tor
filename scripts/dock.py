'''
This script takes two pdb files, the POI and an E3 ligase, and the PROTAC SMILES string to run PRosettaC,
a PROTAC docking software. It creates a directory for the complex based off the PROTAC number and saves the 
best docking results. 
Author: Jordan Harrison
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

def get_ligand_sdf(smiles_list, prefix):
    '''
    Generates a 3D sdf file from a SMILES string
    :param smiles_list: list of smiles strings
    :param prefix: prefix for the output sdf file
    :return: path to the generated sdf file
    '''
    sdf_path = f"{prefix}.sdf"
    writer = Chem.SDWriter(sdf_path)

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"[WARN] Invalid SMILES skipped: {smi}")
            continue

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Embed 3D coords
        status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if status != 0:
            print(f"[WARN] Embedding failed for: {smi}")
            continue

        # Optimize geometry
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except ValueError as e:
            print(f"[WARN] Optimization failed for: {smi} ({e})")
            continue

        # Set original SMILES as property
        mol.SetProp("_Name", smi)

        # Write to SDF
        writer.write(mol)

    writer.close()
    return sdf_path

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

# this isnt used anywhere
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
    
def get_PROTAC(csv_file, output_path='top_smiles.txt', top_n=50):
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

# needs to be fixed, only returns zero. anchor atom inputs are also hard to get, may not even be useful?
def minimize_and_measure(sdf_file, anchor_atom_indices, output_sdf):
    """
    :param sdf_file: Path to the input SDF file (assumes one molecule).
    :param anchor_atom_indices: Atom indices of the two anchors (0-based).

    :return: Tuple of (mol_minimized, distance)
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

# Needs to be fixed
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

# if a pdb gets chopped to bits this is probably the function, need to add in more residues/variants
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
        'SEC', 'PYL', 'HIE', 'HIP', 'ASH',
        'GLH', 'CYM', 'CYX', 'LYN', 'ACE',
        'NME','HOH', 'WAT', 'H2O'  # common water residue names
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
    # Get top 10 PROTAC SMILES from CSV
    protac_smiles_list = get_PROTAC('linkinvent_stage_1.csv',
                                    output_path='top_smiles.txt',
                                    top_n=50)

    # Read warheads from smiles.smi
    with open("smiles.smi", "r") as f:
        raw_smiles = f.read().strip()

    # Split into warheads
    lig1, lig2 = raw_smiles.split('|')

    # Clean '*' for RDKit 3D building
    lig1_clean = lig1.replace('*', '')
    lig2_clean = lig2.replace('*', '')

    # Combine PROTACs + warheads
    all_smiles = protac_smiles_list + [lig1_clean, lig2_clean]

    # Write all 52 ligands into protac.sdf
    sdf_file = get_ligand_sdf(all_smiles, 'protac')

    print(f"[INFO] protac.sdf generated with {len(all_smiles)} ligands")
    print(f"      (50 PROTACs from CSV + 2 warheads from smiles.smi)")

    # Read distance constraints
    with open('input.txt', 'r') as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    print(f'Max Complex Distance: {max_dist}, Min Complex Distance: {min_dist}')

    # Read protein filenames from lig_distances.txt
    proteins = []
    with open('lig_distances.txt', 'r') as f:
        for line in f:
            if ':' in line:
                filename = line.split(':')[0].strip()
                proteins.append(filename)

    # For each protein, create a docking directory and submit a job
    for count, p in enumerate(proteins):
        ternary, ligands = remove_ligand(p)
        job_dir = f"docking_{os.path.splitext(os.path.basename(ternary))[0]}_{count}"
        os.makedirs(job_dir, exist_ok=True)

        # Write config file in job_dir
        config = f'''receptor = {ternary}
ligand = protac.sdf
autobox_ligand = {ternary}
out = docked.sdf.gz
log = log
cnn_scoring = none
num_modes = 5
exhaustiveness = 64
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
#SBATCH --gpus=nvidia_h100_80gb_hbm3_1g.10gb:1
#SBATCH --mem-per-cpu=128M
#SBATCH --time=2:30:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

module purge
module load StdEnv/2023
module load cuda/12.2
module load gcc/12.3
module load gnina/1.3.1

gnina --config config

dock_jobid=$SLURM_JOB_ID

# === md.py ===
echo "Submitting md.py to SLURM after dock.py completes..."
# important that these are loaded after all other jobs are finished due to module dependencies
module load ambertools/25.0
module load amber-pmemd/24.3
md_jobid=$(sbatch --parsable --dependency=afterok:$dock_jobid --job-name=mdpy --output=mdpy.out --error=mdpy.err --wrap="python ../md.py")
echo "Submitted md.py as job $md_jobid (after dock.py)"
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

        try:
            os.system(f'cd {job_dir} && sbatch job.sh')
        except Exception as e:
            print(f"Error: Could not submit job in {job_dir}: {e}")
            exit(1)

main()