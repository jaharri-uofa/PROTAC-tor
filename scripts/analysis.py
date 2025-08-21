'''
Analyis of PROTACtor run and output
Author: Jordan Harrison
'''

import os
import pandas as pd

#Protein protein docking analysis
def pp_compatibility():
    """
    assess compatibility of protein complexes by looking at the rate of pocket to pocket binds compared to other outputs
    """
    score = []
    for file in os.listdir('.'):
        if file.startswith('complex.') and file.endswith('.pdb') and not file.endswith('_nolig.pdb') and not file.endswith('_lig.pdb'):
            num = file.split('.')[1]
            score.append(float(num))
    avg = (sum(score) / len(score)) / 1000
    return avg


def min_max(text):
    """
    Extract the minimum and maximum values from a text file containing distances.
    :param text: Path to the text file.
    """
    min_val = float('inf')
    max_val = float('-inf')
    with open(text, 'r') as f:
        for line in f:
            if ":" not in line:
                continue
            try:
                # extract value after the colon, before "Å"
                value = float(line.split(":")[1].split()[0])
                min_val = min(min_val, value)
                max_val = max(max_val, value)
            except Exception as e:
                print(f"Skipping line: {line.strip()} — {e}")
    return [min_val, max_val]


def get_lysines():
    with open('lysines.txt', 'r') as f:
        out = f.readlines()
        return out

def get_warhead_smiles():
    with open('smiles.smi', 'r') as f:
        out = f.readlines()
        return out
    
def get_total_linkers():
    df = pd.read_csv('linkinvent_stage_1.csv')
    return df.shape[0]

def get_max_linker_score():
    df = pd.read_csv('linkinvent_stage_1.csv')
    return df['Score'].max()

def get_warheads_binding_affinity(csv_file="top5_complexes.csv"):
    df = pd.read_csv(csv_file, header=None)

    # Extract numeric part(s) from the affinity column
    def extract_numbers(x):
        # Split on "-" and keep parts that are numeric
        parts = x.replace("affinity", "").split("-")
        nums = []
        for p in parts:
            try:
                nums.append(float(p))
            except ValueError:
                pass
        return nums

    # Apply extraction and flatten list
    all_affinities = []
    for val in df[1]:
        all_affinities.extend(extract_numbers(str(val)))

    affinities = pd.Series(all_affinities)

    mean_affinity = affinities.mean()
    std_affinity = affinities.std()

    return mean_affinity, std_affinity


def get_top_protac_results(csv_file="top5_protacs.csv"):
    df = pd.read_csv(csv_file, header=None)
    return df

def get_highest_protac_binding_affinity(csv_file="top5_protacs.csv"):
    df = pd.read_csv(csv_file, header=None)
    return df[1].max()

def get_trajectory_rmsd():
    pass

def get_failed_trajectories():
    complex_dirs = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('ternary')]
    failed = []
    for dir in complex_dirs:
        md_dir = os.path.join(dir, "md")
        if not os.path.isdir(md_dir):
            failed.append(dir)
            continue
        if not os.path.isfile(os.path.join(md_dir, 'md.dcd')) or os.path.isfile(os.path.join(md_dir, '06_prod.nc')):
            failed.append(dir)
    return failed

def get_lysine_accessibility_score():
    pass

def get_mmgbsa_scores():
    pass

def main():
    with open('output.txt', 'w') as f:
        f.write(f'Analysis of PROTACtor Output\n'
    f'#Protein-Protein Docking\n'
    f'  compatibility score: {pp_compatibility()}\n'
    f'  Min and Max distance values: {min_max("lig_distances.txt")}\n'
    f'#Linker Generation\n'
    f'  Warhead SMILES (linkinventformat): {get_warhead_smiles()}\n'
    f'  Total number of linkers sampled: {get_total_linkers()}\n'
    f'  Max score of best linker: {get_max_linker_score()}\n'
    f'#PROTAC Docking\n'
    f'  Warheads binding affinity: {get_warheads_binding_affinity()}\n'
    f'  Top PROTAC results: {get_top_protac_results()}\n'
    f'  Highest PROTAC binding affinity: {get_highest_protac_binding_affinity()}\n'
    f'#MD and MM/GBSA\n'
    f'  Trajectory RMSD: {get_trajectory_rmsd()}\n'
    f'  Failed trajectories: {get_failed_trajectories()}\n'
    f'  Lysine accessibility score: {get_lysine_accessibility_score()}\n'
    f'  MM/GBSA scores: {get_mmgbsa_scores()}\n')



main()