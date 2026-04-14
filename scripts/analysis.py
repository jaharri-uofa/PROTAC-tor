#!/usr/bin/env python3
"""
protactor_analysis.py — Unified PROTACtor Analysis
====================================================
Combines static PROTACtor output analysis with MD trajectory analysis.

Stages run in order:
  1. Protein–Protein Docking  — compatibility score, ligand distances
  2. Linker Generation        — warhead SMILES, linker count & top score
  3. PROTAC Docking           — binding affinities, top-ranked complexes
  4. MD Analysis              — RMSD/RMSF via cpptraj, MM/GBSA, ΔΔG

All output is written to output.txt and plots/data are saved to results/.

Author: Jordan Harrison
"""

import os
import re
import sys
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # headless — safe on SLURM nodes
import matplotlib.pyplot as plt


# ── Constants ──────────────────────────────────────────────────────────────────

PROD_TRAJ       = '07-prod.nc'
PRMTOP          = 'complex.prmtop'
CPPTRAJ_TIMEOUT = 900   # seconds


# ── Formatting helpers ─────────────────────────────────────────────────────────

def fmt_angstrom(mean, std=None):
    if mean is None or (isinstance(mean, float) and np.isnan(mean)):
        return 'N/A'
    if std is not None:
        return f'{mean:.2f} ± {std:.2f} Å'
    return f'{mean:.2f} Å'

def fmt_kcal(mean, std=None):
    if mean is None:
        return 'N/A'
    if std is not None:
        return f'{mean:.2f} ± {std:.2f} kcal/mol'
    return f'{mean:.2f} kcal/mol'

def section(title: str) -> str:
    return f'\n## {title}'

def divider() -> str:
    return '=' * 65


# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — Static PROTACtor Output Analysis
# ══════════════════════════════════════════════════════════════════════════════

def pp_compatibility() -> float:
    """
    Assess compatibility of protein complexes by looking at the rate of
    pocket-to-pocket binds compared to other outputs.
    """
    scores = []
    for file in Path('.').glob('complex.*.pdb'):
        if not (file.name.endswith('_nolig.pdb') or file.name.endswith('_lig.pdb')):
            try:
                num = float(file.name.split('.')[1])
                scores.append(num)
            except ValueError:
                continue

    if not scores:
        return float('nan')

    return (sum(scores) / len(scores)) / 1990


def min_max(path: str | Path) -> list[float]:
    """
    Extract the minimum and maximum values from a text file containing distances.
    Lines must contain a colon (":") and values are assumed to be before "Å".
    """
    path = Path(path)
    min_val, max_val = float('inf'), float('-inf')

    with path.open() as f:
        for line in f:
            if ':' not in line:
                continue
            try:
                value = float(line.split(':')[1].split()[0])
                min_val = min(min_val, value)
                max_val = max(max_val, value)
            except Exception:
                continue

    return [min_val, max_val]


def get_lysines() -> list[str]:
    return Path('lysines.txt').read_text().splitlines()


def get_warhead_smiles() -> list[str]:
    return str(Path('smiles.smi').read_text()).split('|')


def get_total_linkers() -> int:
    return pd.read_csv('linkinvent_stage_1.csv').shape[0]


def get_max_linker_score() -> float:
    return pd.read_csv('linkinvent_stage_1.csv')['Score'].max()


def _extract_affinities(values) -> list[float]:
    """Helper: parse affinity strings like 'affinity-13.9-13.7' into floats."""
    affinities = []
    for val in values:
        for part in str(val).replace('affinity', '').split('-'):
            try:
                affinities.append(float(part))
            except ValueError:
                continue
    return affinities


def get_warheads_binding_affinity(csv_file: str | Path = 'top5_complexes.csv') -> tuple[float, float]:
    df = pd.read_csv(csv_file, header=None)
    affinities = pd.Series(_extract_affinities(df[1]))
    return affinities.mean(), affinities.std()


def display_top_results(csv_path='top5_complexes.csv', top_n=5) -> str:
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()

    if 'smiles' not in df.columns or 'affinity' not in df.columns:
        raise ValueError(
            f"CSV must contain 'smiles' and 'affinity' columns. Found: {list(df.columns)}"
        )

    df['affinity'] = pd.to_numeric(df['affinity'], errors='coerce').round(2)
    df = df.sort_values(by='affinity').reset_index(drop=True)
    df.insert(0, 'Rank', range(1, len(df) + 1))
    df = df[['Rank', 'smiles', 'affinity']]
    return df.head(top_n).to_string(index=False)


def get_highest_protac_binding_affinity(csv_file: str | Path = 'top5_complexes.csv') -> float:
    df = pd.read_csv(csv_file, header=None)
    return df[1].max()


def get_lysine_accessibility_score() -> dict:
    with open('lysines.txt', 'r') as f:
        lysines = [line.strip() for line in f.readlines()]

    pdb = []
    ligase = []
    poi = []
    new_pdb = False
    centroids = []
    lys_dist = {lys: [] for lys in lysines}

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
                            if atom_num == atom - 1 and not new_pdb:
                                ligase.append(_)
                                atom_num = atom
                            elif atom_num != atom - 1:
                                poi.append(_)
                                atom_num = atom
                                new_pdb = True
                            elif atom_num == atom - 1 and new_pdb:
                                poi.append(_)
                                atom_num = atom
                            else:
                                print('Something is wrong with the pdb, try again :)')
                    centroid_ligase = centroid(ligase)
                    centroid_poi    = centroid(poi)
                    centroids.append((centroid_ligase, centroid_poi))
                    for _ in poi:
                        res_id    = f"{_[17:20].strip()}_{int(_[22:26].strip())}"
                        atom_type = _[12:16].strip()
                        if res_id in lysines and atom_type == 'NZ':
                            dist = distance(
                                centroid_ligase,
                                np.array([
                                    float(_[30:38].strip()),
                                    float(_[38:46].strip()),
                                    float(_[46:54].strip()),
                                ])
                            )
                            lys_dist[res_id].append(dist)
                    ligase  = []
                    poi     = []
                pdb     = []
                new_pdb = False

    return {
        k: (np.mean(v), np.std(v)) if v else (float('nan'), float('nan'))
        for k, v in lys_dist.items()
    }


def distance(lig1, lig2) -> float:
    """Calculate Euclidean distance between two coordinate arrays."""
    return np.linalg.norm(lig1 - lig2)


def centroid(pdb: list) -> np.ndarray:
    """Calculate the centroid of atomic coordinates from PDB lines."""
    coords = []
    for line in pdb:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except ValueError:
                continue

    if not coords:
        return np.array([0.0, 0.0, 0.0])

    return np.mean(np.array(coords), axis=0)


# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — MD Analysis
# ══════════════════════════════════════════════════════════════════════════════

# ── Directory discovery ────────────────────────────────────────────────────────

def find_completed_md_dirs(base: Path) -> dict:
    """
    Walk base looking for ternary_complex* and ternary_complex_control dirs
    whose md/ subdirectory contains a completed production trajectory.

    Returns:
        {
          'protac':  [Path, ...],   # ternary_complex1/md, ...
          'control': Path | None,   # ternary_complex_control/md or None
        }
    """
    protac_dirs = []
    control_dir = None

    for d in sorted(base.iterdir()):
        if not d.is_dir():
            continue
        md = d / 'md'
        if not md.is_dir():
            continue
        if not (md / PROD_TRAJ).exists():
            continue

        if d.name == 'ternary_complex_control':
            control_dir = md
        elif d.name.startswith('ternary_complex'):
            protac_dirs.append(md)

    return {'protac': protac_dirs, 'control': control_dir}


# ── Ligand residue name ────────────────────────────────────────────────────────

def read_ligand_resname(md_dir: Path) -> str | None:
    """Return the first ligand residue name from ligand_resname.txt."""
    for p in [
        md_dir.parent / 'ligand_resname.txt',
        md_dir.parent / 'prep' / 'ligand_resname.txt',
        md_dir.parent.parent / 'ligand_resname.txt',
    ]:
        if p.exists():
            lines = p.read_text().strip().splitlines()
            return lines[0].strip() if lines else None
    return None


# ── cpptraj RMSD ──────────────────────────────────────────────────────────────

def write_rmsd_cpptraj(md_dir: Path, lig_resname: str | None) -> Path:
    """
    Write a cpptraj script that computes:
      - Backbone (Cα/C/N) RMSD vs frame 1    → rmsd1.dat
      - Ligand heavy-atom RMSD (no-fit)       → rmsd_ligand.dat  (if resname known)
      - Per-residue backbone RMSF             → rmsf_byres.dat
    """
    lines = [
        f'parm {PRMTOP}',
        f'trajin {PROD_TRAJ}',
        'autoimage',
        'strip :WAT,Na+,Cl-',
        'rms first @CA,C,N out rmsd1.dat mass',
    ]
    if lig_resname:
        lines.append(
            f'rms first :{lig_resname}&!@H= out rmsd_ligand.dat mass nofit'
        )
    lines += [
        'atomicfluct @CA,C,N byres out rmsf_byres.dat',
        'run',
        'quit',
    ]
    inp = md_dir / 'rmsd_analysis.in'
    inp.write_text('\n'.join(lines) + '\n')
    return inp


def run_cpptraj(md_dir: Path, inp: Path) -> bool:
    """Execute cpptraj inside md_dir. Returns True on success."""
    try:
        r = subprocess.run(
            ['cpptraj', '-i', inp.name],
            cwd=str(md_dir),
            capture_output=True, text=True,
            timeout=CPPTRAJ_TIMEOUT,
        )
        if r.returncode != 0:
            print(f'  [cpptraj WARN] {md_dir.parent.name}: {r.stderr[:300]}')
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f'  [cpptraj TIMEOUT] {md_dir.parent.name}')
        return False
    except FileNotFoundError:
        print('  [cpptraj ERROR] cpptraj not found in PATH.')
        return False


def parse_two_col_dat(filepath: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse a whitespace-delimited two-column cpptraj output file."""
    xs, ys = [], []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
                except ValueError:
                    continue
    return np.array(xs, dtype=float), np.array(ys, dtype=float)


# ── Plotting ───────────────────────────────────────────────────────────────────

def save_rmsd_plot(x: np.ndarray, y: np.ndarray,
                   title: str, ylabel: str, out: Path):
    """Save a time-series plot with an overlaid rolling mean."""
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(x, y, color='steelblue', linewidth=0.7, alpha=0.7, label='RMSD')

    w = max(1, len(y) // 50)
    if len(y) >= w * 2:
        roll   = np.convolve(y, np.ones(w) / w, mode='valid')
        roll_x = x[:len(roll)]
        ax.plot(roll_x, roll, color='tomato', linewidth=1.5,
                label=f'Rolling mean (w={w})')

    ax.axhline(np.mean(y), color='grey', linewidth=0.8,
               linestyle='--', label=f'Mean={np.mean(y):.2f} Å')
    ax.set_xlabel('Frame')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    fig.savefig(str(out), dpi=150)
    plt.close(fig)


def save_rmsf_plot(residues: np.ndarray, rmsf: np.ndarray,
                   title: str, out: Path):
    """Save a per-residue RMSF bar/line plot."""
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.fill_between(residues, rmsf, alpha=0.5, color='teal')
    ax.plot(residues, rmsf, color='teal', linewidth=0.8)
    ax.set_xlabel('Residue number')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title(title)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    fig.savefig(str(out), dpi=150)
    plt.close(fig)


# ── MM/GBSA parsing ────────────────────────────────────────────────────────────

_DELTA_TOTAL = re.compile(
    r'DELTA\s+TOTAL\s+([-+]?\d+\.\d+)\s+(\d+\.\d+)'
)


def parse_mmgbsa_dat(dat_file: Path) -> tuple[float, float] | tuple[None, None]:
    """
    Parse an MMPBSA.py output file.
    Returns (DELTA_TOTAL_mean, DELTA_TOTAL_std) or (None, None).
    """
    if not dat_file.exists():
        return None, None
    with open(dat_file) as f:
        for line in f:
            m = _DELTA_TOTAL.search(line)
            if m:
                return float(m.group(1)), float(m.group(2))
    return None, None


def collect_mmgbsa(md_dir: Path) -> dict[str, tuple[float, float]]:
    """
    Find all out_ligand*.dat files in md_dir and parse them.
    Returns {label: (mean, std)}.
    """
    results = {}
    for candidate in ['out_ligand.dat', 'out_ligand1.dat', 'out_ligand2.dat']:
        mean, std = parse_mmgbsa_dat(md_dir / candidate)
        if mean is not None:
            results[candidate.replace('.dat', '')] = (mean, std)
    return results


# ── ΔΔG calculation ────────────────────────────────────────────────────────────

def compute_ddg(protac_mmgbsa: dict, ctrl_mmgbsa: dict) -> dict:
    """
    ΔΔG = ΔG(PROTAC ternary) − ΔG(control warhead-only complex).
    """
    out = {}

    if 'out_ligand' in protac_mmgbsa:
        p_mean, p_std = protac_mmgbsa['out_ligand']
    elif 'out_ligand1' in protac_mmgbsa:
        p_mean, p_std = protac_mmgbsa['out_ligand1']
    else:
        return out

    if 'out_ligand' in ctrl_mmgbsa:
        c_mean, c_std = ctrl_mmgbsa['out_ligand']
    else:
        return out

    ddg     = p_mean - c_mean
    ddg_std = np.sqrt(p_std**2 + c_std**2)

    out['PROTAC_deltaG']  = (p_mean, p_std)
    out['Control_deltaG'] = (c_mean, c_std)
    out['ddG_linker']     = (ddg, ddg_std)
    return out


# ── Per-directory MD analysis ──────────────────────────────────────────────────

def analyse_md_dir(md_dir: Path, label: str, results_dir: Path) -> dict:
    """
    Full analysis for one md/ directory.
    Runs cpptraj, generates plots, parses MM/GBSA.
    Returns a summary dict.
    """
    print(f'\n  Analysing: {md_dir.parent.name}/{md_dir.name}')
    summary = {'label': label}

    lig_res = read_ligand_resname(md_dir)

    inp = write_rmsd_cpptraj(md_dir, lig_res)
    run_cpptraj(md_dir, inp)

    rmsd_dat = md_dir / 'rmsd1.dat'
    if rmsd_dat.exists():
        frames, rmsd = parse_two_col_dat(rmsd_dat)
        summary['bb_rmsd_mean'] = float(np.mean(rmsd))
        summary['bb_rmsd_std']  = float(np.std(rmsd))
        save_rmsd_plot(
            frames, rmsd,
            title=f'Backbone RMSD — {label}',
            ylabel='RMSD (Å)',
            out=results_dir / f'{label}_backbone_rmsd.png',
        )
        shutil.copy(rmsd_dat, results_dir / f'{label}_backbone_rmsd.dat')
    else:
        print(f'    [WARN] rmsd1.dat missing — cpptraj may have failed.')

    lig_dat = md_dir / 'rmsd_ligand.dat'
    if lig_dat.exists():
        frames, lig_rmsd = parse_two_col_dat(lig_dat)
        summary['lig_rmsd_mean'] = float(np.mean(lig_rmsd))
        summary['lig_rmsd_std']  = float(np.std(lig_rmsd))
        save_rmsd_plot(
            frames, lig_rmsd,
            title=f'Ligand RMSD — {label}',
            ylabel='RMSD (Å)',
            out=results_dir / f'{label}_ligand_rmsd.png',
        )
        shutil.copy(lig_dat, results_dir / f'{label}_ligand_rmsd.dat')

    rmsf_dat = md_dir / 'rmsf_byres.dat'
    if rmsf_dat.exists():
        res, rmsf = parse_two_col_dat(rmsf_dat)
        summary['rmsf_max_res'] = int(res[np.argmax(rmsf)]) if len(res) else None
        summary['rmsf_max_val'] = float(np.max(rmsf))       if len(rmsf) else None
        save_rmsf_plot(
            res, rmsf,
            title=f'Per-residue RMSF — {label}',
            out=results_dir / f'{label}_rmsf.png',
        )
        shutil.copy(rmsf_dat, results_dir / f'{label}_rmsf_byres.dat')

    mmgbsa = collect_mmgbsa(md_dir)
    summary['mmgbsa'] = mmgbsa
    for key in mmgbsa:
        src = md_dir / f'{key}.dat'
        if src.exists():
            shutil.copy(src, results_dir / f'{label}_{key}.dat')

    return summary


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    base        = Path.cwd()
    results_dir = base / 'results'
    results_dir.mkdir(exist_ok=True)

    print(divider())
    print('PROTACtor Analysis')
    print(divider())

    lines = [
        divider(),
        'PROTACtor Analysis Report',
        divider(),
    ]

    # ── Part 1: Protein–Protein Docking ───────────────────────────────────────
    lines.append(section('Protein–Protein Docking'))
    try:
        lines.append(f'  Compatibility score          : {pp_compatibility():.4f}')
    except Exception as e:
        lines.append(f'  Compatibility score          : ERROR ({e})')

    try:
        lo, hi = min_max('lig_distances.txt')
        lines.append(f'  Ligand distance range        : {lo:.2f} – {hi:.2f} Å')
    except Exception as e:
        lines.append(f'  Ligand distance range        : ERROR ({e})')

    # ── Part 2: Linker Generation ──────────────────────────────────────────────
    lines.append(section('Linker Generation'))
    try:
        smiles = get_warhead_smiles()
        lines.append(f'  Warhead SMILES (1)           : {smiles[0].strip()}')
        if len(smiles) > 1:
            lines.append(f'  Warhead SMILES (2)           : {smiles[1].strip()}')
    except Exception as e:
        lines.append(f'  Warhead SMILES               : ERROR ({e})')

    try:
        lines.append(f'  Total linkers sampled        : {get_total_linkers()}')
    except Exception as e:
        lines.append(f'  Total linkers sampled        : ERROR ({e})')

    try:
        lines.append(f'  Best linker score            : {get_max_linker_score():.4f}')
    except Exception as e:
        lines.append(f'  Best linker score            : ERROR ({e})')

    # ── Part 3: PROTAC Docking ─────────────────────────────────────────────────
    lines.append(section('PROTAC Docking'))
    try:
        mean_aff, std_aff = get_warheads_binding_affinity()
        lines.append(
            f'  Warhead binding affinity     : '
            f'{fmt_kcal(mean_aff, std_aff)}'
        )
    except Exception as e:
        lines.append(f'  Warhead binding affinity     : ERROR ({e})')

    try:
        lines.append(
            f'  Highest PROTAC affinity      : '
            f'{get_highest_protac_binding_affinity()}'
        )
    except Exception as e:
        lines.append(f'  Highest PROTAC affinity      : ERROR ({e})')

    try:
        table = display_top_results('top5_complexes.csv')
        lines.append('\n  Top PROTAC complexes:')
        for row in table.splitlines():
            lines.append(f'    {row}')
    except Exception as e:
        lines.append(f'  Top PROTAC complexes         : ERROR ({e})')

    # ── Part 4: MD Analysis ────────────────────────────────────────────────────
    print(section('MD Analysis'))
    dirs = find_completed_md_dirs(base)
    protac_md_dirs = dirs['protac']
    control_md_dir = dirs['control']

    print(f'  PROTAC MD directories  : {len(protac_md_dirs)}')
    print(f'  Control MD directory   : '
          f'{"found" if control_md_dir else "not found"}')

    lines.append(section('MD Analysis'))

    if not protac_md_dirs and control_md_dir is None:
        lines.append('  No completed MD runs found (looking for 07-prod.nc).')
        lines.append('  Re-run this script once MD completes.')
    else:
        protac_summaries = []
        for md_dir in protac_md_dirs:
            label = f'{md_dir.parent.parent.name}_{md_dir.parent.name}'
            s = analyse_md_dir(md_dir, label, results_dir)
            protac_summaries.append(s)

        ctrl_summary = None
        if control_md_dir:
            ctrl_summary = analyse_md_dir(control_md_dir, 'control', results_dir)

        # PROTAC ternary complex summaries
        if protac_summaries:
            lines.append('\n### PROTAC Ternary Complexes')
            for s in protac_summaries:
                lines.append(f"\n  Complex          : {s['label']}")
                if 'bb_rmsd_mean' in s:
                    lines.append(
                        f"  Backbone RMSD    : "
                        f"{fmt_angstrom(s.get('bb_rmsd_mean'), s.get('bb_rmsd_std'))}"
                    )
                if 'lig_rmsd_mean' in s:
                    lines.append(
                        f"  Ligand RMSD      : "
                        f"{fmt_angstrom(s.get('lig_rmsd_mean'), s.get('lig_rmsd_std'))}"
                    )
                if s.get('rmsf_max_res'):
                    lines.append(
                        f"  Max RMSF residue : residue {s['rmsf_max_res']} "
                        f"({s['rmsf_max_val']:.2f} Å)"
                    )
                for key, (m, sd) in s.get('mmgbsa', {}).items():
                    lines.append(f"  {key} ΔG         : {fmt_kcal(m, sd)}")

            all_dg = [
                s['mmgbsa'].get('out_ligand', (None, None))[0]
                for s in protac_summaries
                if 'out_ligand' in s.get('mmgbsa', {})
            ]
            if all_dg:
                lines.append(
                    f"\n  Average PROTAC ΔG (all complexes) : "
                    f"{np.mean(all_dg):.2f} kcal/mol"
                )

        # Control summary
        if ctrl_summary:
            lines.append('\n### Control (Warheads Only)')
            if 'bb_rmsd_mean' in ctrl_summary:
                lines.append(
                    f"  Backbone RMSD    : "
                    f"{fmt_angstrom(ctrl_summary.get('bb_rmsd_mean'), ctrl_summary.get('bb_rmsd_std'))}"
                )
            for key, (m, sd) in ctrl_summary.get('mmgbsa', {}).items():
                lines.append(f"  {key} ΔG         : {fmt_kcal(m, sd)}")

        # ΔΔG
        if protac_summaries and ctrl_summary:
            lines.append('\n### Linker Contribution: ΔΔG = ΔG(PROTAC) − ΔG(control)')
            lines.append(
                '  Interpretation: negative ΔΔG = linker stabilises the ternary\n'
                '  complex relative to the warhead-only control.'
            )
            avg_protac_mmgbsa: dict = {}
            for s in protac_summaries:
                for k, (m, sd) in s.get('mmgbsa', {}).items():
                    avg_protac_mmgbsa.setdefault(k, []).append((m, sd))
            avg_p = {
                k: (np.mean([v[0] for v in vs]),
                    np.sqrt(np.mean([v[1]**2 for v in vs])))
                for k, vs in avg_protac_mmgbsa.items()
            }

            ddg = compute_ddg(avg_p, ctrl_summary.get('mmgbsa', {}))
            if ddg:
                lines.append(f"  PROTAC ΔG                         : {fmt_kcal(*ddg['PROTAC_deltaG'])}")
                lines.append(f"  Control ΔG                        : {fmt_kcal(*ddg['Control_deltaG'])}")
                lines.append(f"  ΔΔG (linker contribution)         : {fmt_kcal(*ddg['ddG_linker'])}")
            else:
                lines.append(
                    '  [WARN] Could not compute ΔΔG — MM/GBSA data missing '
                    'from one or both runs.'
                )

    # ── Output files ────────────────────────────────────────────────────────────
    lines.append(section('Output Files'))
    lines.append(f'  Plots and data : {results_dir}/')
    lines.append('  RMSD plots     : *_backbone_rmsd.png, *_ligand_rmsd.png')
    lines.append('  RMSF plots     : *_rmsf.png')
    lines.append('  MM/GBSA data   : *_out_ligand*.dat')

    # ── Write report ─────────────────────────────────────────────────────────
    report   = '\n'.join(map(str, lines))
    out_file = base / 'output.txt'
    out_file.write_text(report + '\n')

    print(report)
    print(f'\nReport written to {out_file}')
    print(f'All analysis files saved to {results_dir}/')

    # Copy key input/output files into results/
    for fname in ['output.txt', 'top5_complexes.csv', 'lig_distances.txt', 'smiles.smi']:
        src = base / fname
        if src.exists():
            shutil.copy(src, results_dir / fname)


if __name__ == '__main__':
    main()