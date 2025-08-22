#!/usr/bin/env python3
"""
Analysis of PROTACtor run and output
Author: Jordan Harrison
"""

import pandas as pd
from pathlib import Path


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

    return (sum(scores) / len(scores)) / 1000


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
    return Path("smiles.smi").read_text().splitlines()


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


def get_top_protac_results(csv_file: str | Path = "top5_complexes.csv") -> pd.DataFrame:
    return pd.read_csv(csv_file, header=None)


def get_highest_protac_binding_affinity(csv_file: str | Path = "top5_complexes.csv") -> float:
    df = pd.read_csv(csv_file, header=None)
    return df[1].max()


def get_trajectory_rmsd():
    return None


def get_failed_trajectories() -> list[str]:
    """Check all `ternary*` directories for missing MD outputs."""
    failed = []
    for d in Path(".").glob("ternary*"):
        if not d.is_dir():
            continue
        md_dir = d / "md"
        if not md_dir.exists():
            failed.append(d.name)
            continue
        if not (md_dir / "06_prod.nc").exists():
            failed.append(d.name)
    return failed


def get_lysine_accessibility_score():
    return None


def get_mmgbsa_scores():
    return None


def main():
    output_lines = [
        "Analysis of PROTACtor Output",
        "#Protein-Protein Docking",
        f"  Compatibility score: {pp_compatibility()}",
        f"  Min and Max distance values: {min_max('lig_distances.txt')}",
        "#Linker Generation",
        f"  Warhead SMILES (linkinvent format): {get_warhead_smiles()}",
        f"  Total number of linkers sampled: {get_total_linkers()}",
        f"  Max score of best linker: {get_max_linker_score()}",
        "#PROTAC Docking",
        f"  Warheads binding affinity: {get_warheads_binding_affinity()}",
        f"  Top PROTAC results:\n{get_top_protac_results()}",
        f"  Highest PROTAC binding affinity: {get_highest_protac_binding_affinity()}",
        "#MD and MM/GBSA",
        f"  Trajectory RMSD: {get_trajectory_rmsd()}",
        f"  Failed trajectories: {get_failed_trajectories()}",
        f"  Lysine accessibility score: {get_lysine_accessibility_score()}",
        f"  MM/GBSA scores: {get_mmgbsa_scores()}",
    ]

    Path("output.txt").write_text("\n".join(map(str, output_lines)))


if __name__ == "__main__":
    main()
