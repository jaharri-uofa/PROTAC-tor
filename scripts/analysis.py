#!/usr/bin/env python3
"""
Analysis of PROTACtor run and output
Author: Jordan Harrison
"""

import pandas as pd
from pathlib import Path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os


def pp_compatibility() -> float:
    """
    Assess compatibility of protein complexes by looking at the rate of
    pocket-to-pocket binds compared to other outputs.
    """
    scores = []
    for file in Path(".").glob("complex.*.pdb"):
        if not (file.name.endswith("_nolig.pdb") or file.name.endswith("_lig.pdb")):
            try:
                num = float(file.name.split(".")[1])
                scores.append(num)
            except ValueError:
                continue

    if not scores:
        return float("nan")

    return (sum(scores) / len(scores)) / 1990

def min_max(path: str | Path) -> list[float]:
    """
    Extract the minimum and maximum values from a text file containing distances.
    Lines must contain a colon (":") and values are assumed to be before "Å".
    """
    path = Path(path)
    min_val, max_val = float("inf"), float("-inf")

    with path.open() as f:
        for line in f:
            if ":" not in line:
                continue
            try:
                value = float(line.split(":")[1].split()[0])
                min_val = min(min_val, value)
                max_val = max(max_val, value)
            except Exception:
                continue

    return [min_val, max_val]


def get_lysines() -> list[str]:
    return Path("lysines.txt").read_text().splitlines()


def get_warhead_smiles() -> list[str]:
    return str(Path("smiles.smi").read_text()).split('|')


def get_total_linkers() -> int:
    return pd.read_csv("linkinvent_stage_1.csv").shape[0]


def get_max_linker_score() -> float:
    return pd.read_csv("linkinvent_stage_1.csv")["Score"].max()


def _extract_affinities(values) -> list[float]:
    """Helper: parse affinity strings like 'affinity-13.9-13.7' into floats."""
    affinities = []
    for val in values:
        for part in str(val).replace("affinity", "").split("-"):
            try:
                affinities.append(float(part))
            except ValueError:
                continue
    return affinities


def get_warheads_binding_affinity(csv_file: str | Path = "top5_complexes.csv") -> tuple[float, float]:
    df = pd.read_csv(csv_file, header=None)
    affinities = pd.Series(_extract_affinities(df[1]))
    return affinities.mean(), affinities.std()

import pandas as pd

def display_top_results(csv_path = "top5_complexes.csv", top_n=5):
    df = pd.read_csv(csv_path)

    # Clean column names (remove leading/trailing spaces and lowercase)
    df.columns = df.columns.str.strip().str.lower()

    # Ensure the required columns exist
    if "smiles" not in df.columns or "affinity" not in df.columns:
        raise ValueError(f"CSV must contain 'smiles' and 'affinity' columns. Found: {list(df.columns)}")

    # Convert affinity to float and round
    df["affinity"] = pd.to_numeric(df["affinity"], errors="coerce").round(2)

    # Sort by affinity (assuming lower is better binding)
    df = df.sort_values(by="affinity").reset_index(drop=True)

    # Add rank
    df.insert(0, "Rank", range(1, len(df) + 1))

    # Keep only smiles and affinity
    df = df[["Rank", "smiles", "affinity"]]

    # Get top N results
    return df.head(top_n).to_string(index=False)

def get_highest_protac_binding_affinity(csv_file: str | Path = "top5_complexes.csv") -> float:
    df = pd.read_csv(csv_file, header=None)
    return df[1].max()


def get_trajectory_rmsd(dirs):
    # cd into the md directory and read off the rmsd.dat file and create a figure
    # need to run cpptraj rmsd
    # take the output file which is a tab seperated file with two columns, frame number, and rmsd
    #get into the directory 'md'
    for dir in dirs:
        os.chdir('md')
        subprocess.run(["cpptraj", "-i", "rmsd.in"], check=True)
        frame = []
        rmsd = []
        with open("rmsd1.dat", 'r') as f:
            for line in f:
                _ = line.strip().split()
                frame.append(int(_[0]))
                rmsd.append(float(_[1]))

        plt.plot(np.array(frame), np.array(rmsd), label = 'RMSD', color = 'blue') # To add multiple functions to the same graph repeat this line with different arguments
        plt.xlabel("Frame #")
        plt.ylabel("RMSD (Å)")
        plt.title("RMSD over Time")
        plt.grid(False) # Set to false or delete line if you do not want a grid
        plt.legend()
        plt.show()

def distance(lig1, lig2):
    '''
    Calculates distance between two points (typically ligands)
    :param lig1: Coordinates of the first ligand
    :param lig2: Coordinates of the second ligand
    :return: Distance between the two ligands
    '''
    return np.linalg.norm(lig1 - lig2)

def get_lysine_accessibility_score():
    with open('lysines.txt', 'r') as f:
        lysines = [line.strip() for line in f.readlines()]

    pdb = []
    ligase = []
    poi = []
    new_pdb = False
    centroids = []
    lys_dist = {lys: [] for lys in lysines}

    # needs to be in md directorys
    with open('ensemble.pdb', 'r') as f:
        for line in f:
            pdb.append(line)
            if line.startswith('TER'):
                with open('ensemble_frag.pdb', 'w') as out:
                    out.writelines(pdb)
                with open('ensemble_frag.pdb', 'r') as frag:
                    lines = frag.readlines()
                    atom_num = 0
                    for _ in lines:
                        if _.startswith('ATOM'):
                            atom = int(_[6:11].strip())
                            if atom_num == atom - 1 and not new_pdb:   # this is the first chain
                                ligase.append(_)
                                atom_num = atom
                            elif atom_num != atom -1:  # This means we have hit a new chain
                                poi.append(_)
                                atom_num = atom
                                new_pdb = True
                            elif atom_num == atom -1 and new_pdb:  # this is the second chain
                                poi.append(_)
                                atom_num = atom
                            else:
                                print('Something is wrong with the pdb, try again :)')
                    centroid_ligase = centroid(ligase)
                    centroid_poi = centroid(poi)
                    centroids.append((centroid_ligase, centroid_poi))
                    for _ in poi:
                        res_id = f"{_[17:20].strip()}_{int(_[22:26].strip())}"
                        atom_type = _[12:16].strip()
                        if res_id in lysines and atom_type == 'NZ':
                            dist = distance(centroid_ligase, np.array([float(_[30:38].strip()), float(_[38:46].strip()), float(_[46:54].strip())]))
                            lys_dist[res_id].append(dist)
                    ligase = []
                    poi = []
                pdb = []
                new_pdb = False
    return {k: (np.mean(v), np.std(v)) if v else (float('nan'), float('nan')) for k, v in lys_dist.items()}

def centroid(pdb: list):
    """
    Calculate the centroid of a set of atomic coordinates from a PDB file.
    :param pdb: List of lines from a PDB file.
    :return: Numpy array representing the centroid coordinates (x, y, z).
    """
    coords = []
    for line in pdb:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except ValueError:
                continue

    if not coords:
        return np.array([0.0, 0.0, 0.0])

    coords_array = np.array(coords)
    return np.mean(coords_array, axis = 0)

def get_mmgbsa_scores():
    # cd into the md directorys and read off the out1.dat file
    # this is what I would do if I had an out1.dat file
    return None


def main():
    dirs = [d for d in Path(".").iterdir() if d.is_dir() and d.name.startswith("ternary")]
    
    print(dirs)
    output_lines = [
        "Analysis of PROTACtor Output",
        "#Protein-Protein Docking",
        f"  Compatibility score: {pp_compatibility()}",
        f"  Min and Max distance values: {min_max('lig_distances.txt')}",
        "#Linker Generation",
        f"  Warhead SMILES (linkinvent format): {get_warhead_smiles()[0]} | {get_warhead_smiles()[1]}",
        f"  Total number of linkers sampled: {get_total_linkers()}",
        f"  Max score of best linker: {get_max_linker_score()}",
        "#PROTAC Docking",
        f"  Warheads binding affinity: {get_warheads_binding_affinity()[0]} ± {get_warheads_binding_affinity()[1]}",
        f"  Top PROTAC results:\n{display_top_results('top5_complexes.csv')}",
        f"  Highest PROTAC binding affinity: {get_highest_protac_binding_affinity()}",
    ]

    Path("output.txt").write_text("\n".join(map(str, output_lines)))


if __name__ == "__main__":
    main()
