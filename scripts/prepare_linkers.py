#!/usr/bin/env python3
"""
prepare_tl_dataset.py
=====================
Builds a PROTAC linker SMILES file for Transfer Learning with Link-INVENT.

Two modes:
  1. --protacdb  path/to/PROTAC-DB_export.csv
        Parses a CSV downloaded from https://protacdb.cadd.zju.edu.cn/
        and extracts the linker SMILES column (expects 'Linker Smiles').

  2. (no flag)
        Uses a built-in curated set of ~40 drug-like PROTAC linkers
        covering PEG, alkyl, piperazine, triazole, and heterocyclic chemotypes
        commonly found in published PROTACs.

Output: one linker SMILES per line with `*` attachment points.
        e.g.  *CC(=O)NCCOCCOCCN*
              *c1ccncc1CC(=O)N*

Usage:
    python prepare_tl_dataset.py                       # uses built-in set
    python prepare_tl_dataset.py --protacdb data.csv   # parse PROTAC-DB export
    python prepare_tl_dataset.py --out my_linkers.smi  # custom output path

Author: Jordan Harrison
"""

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# ── Built-in curated linker set ────────────────────────────────────────────────
# Representative linkers from published PROTACs (ARV-110, ARV-471, MZ1,
# dBET6, AT1, ACBI1, and others). Attachment points marked with *.
# These span a range of lengths (4–20 heavy atoms) and chemotypes.

BUILTIN_LINKERS = [
    # ── PEG-based (flexible, polar) ──
    "*CC(=O)NCCOCCOCCN*",
    "*CC(=O)NCCOCCOCCOCCN*",
    "*CNCCOCCOCCN*",
    "*CC(=O)NCCOCCOCC(=O)N*",
    "*CCOCCOCCN*",
    "*CCOCCOCCOCCOCC(=O)N*",
    "*CC(=O)NCCOCCOCCOCCOCCOCC(=O)N*",
    # ── Alkyl chains ──
    "*CCCCN*",
    "*CCCCCN*",
    "*CCCCCCN*",
    "*CC(=O)NCCCCC(=O)N*",
    "*CCNCCN*",
    # ── Piperazine/piperidine ──
    "*CN1CCN(CC1)C(=O)*",
    "*CN1CCCC(C1)N*",
    "*C(=O)N1CCN(CC1)CC(=O)N*",
    "*CN1CCN(CCOCCOCCOCCOCC(=O)N*",
    "*C1CN(CCN1CC(=O)N*)C(=O)*",
    # ── Amide/urea-containing ──
    "*CC(=O)NCCC(=O)N*",
    "*NC(=O)CCNC(=O)CCC(=O)N*",
    "*CC(=O)NCC(=O)N*",
    "*C(=O)NCCCCCNC(=O)*",
    "*C(=O)NCCCCNC(=O)*",
    # ── Triazole-containing ──
    "*Cc1cn(CC(=O)N*)nn1",
    "*CC(=O)Nn1ccnn1*",
    "*CN1C=NN=C1CCC(=O)N*",
    # ── Aromatic spacers ──
    "*Cc1ccc(CC(=O)N*)cc1",
    "*CC(=O)NCc1cccc(CC(=O)N*)c1",
    "*c1cc(CN*)ccc1C(=O)N*",
    "*CC(=O)Nc1cccc(CC(=O)N*)c1",
    # ── Sulfonamide/sulfone ──
    "*CCS(=O)(=O)CCN*",
    "*CCNS(=O)(=O)CCN*",
    # ── Mixed short/medium (common in clinical PROTACs) ──
    "*CC(=O)NCCO*",
    "*CCOCCN*",
    "*CCOCCOCCO*",
    "*CC(=O)NCCN*",
    "*CNCCN*",
    "*CCNCC(=O)N*",
    # ── Longer (for wider anchor distances) ──
    "*CC(=O)NCCOCCOCCOCCOCCN*",
    "*CCOCCOCCOCCOCCOCC(=O)N*",
    "*CC(=O)NCCCCCCCN*",
]

# ── REINVENT Allowed Characters  ───────────────────────────────────────────────

char = ({'3', 's', '[N+]', '[n+]', '#', '6', 'C', '5', '*', '1', '(', 'Cl', '[O]', '=', 'n', 'O', 'Br', '[nH]', '|', '<pad>', 
         '-', 'N', ')', 'o', '$', 'F', '2', 'S', '[S+]', '^', '[O-]', '4', '[s+]', 'c'})

# ── Validation ─────────────────────────────────────────────────────────────────

def has_allowed_characters(smiles: str, allowed_chars: set) -> bool:
    """
    Check if a SMILES string contains only allowed characters.
    Handles both single characters and multi-character tokens like [N+], Cl, Br.
    
    :param smiles: SMILES string to check
    :param allowed_chars: Set of allowed character tokens
    :return: True if all characters are allowed, False otherwise
    """
    i = 0
    while i < len(smiles):
        # Handle bracketed atoms like [N+], [O-], [nH]
        if smiles[i] == '[':
            j = i + 1
            while j < len(smiles) and smiles[j] != ']':
                j += 1
            if j < len(smiles):
                token = smiles[i:j+1]
                if token not in allowed_chars:
                    return False
                i = j + 1
            else:
                return False  # Unmatched bracket
        # Handle two-character tokens like Cl, Br
        elif i + 1 < len(smiles) and smiles[i:i+2] in allowed_chars:
            i += 2
        # Handle single characters
        elif smiles[i] in allowed_chars:
            i += 1
        else:
            return False
    return True

def validate_linker_smiles(smiles_list: list[str]) -> list[str]:
    """
    Filter to valid SMILES with exactly two * attachment points
    and only allowed characters. Returns cleaned list.
    """
    if not HAS_RDKIT:
        print('[WARN] RDKit not available — skipping SMILES validation.')
        return smiles_list

    valid = []
    for smi in smiles_list:
        smi = smi.strip()
        if not smi or smi.startswith('#'):
            continue
            '''
        # Check attachment point count
        if smi.count('*') != 2:
            print(f'  [SKIP] Not exactly 2 attachment points: {smi}')
            continue
            '''
        # Check for forbidden characters
        if not has_allowed_characters(smi, char[0]):
            print(f'  [SKIP] Contains forbidden characters: {smi}')
            continue
        # Check RDKit can parse it
        mol = Chem.MolFromSmiles(smi.replace('*', '[*]'))
        if mol is None:
            print(f'  [SKIP] RDKit could not parse: {smi}')
            continue
        valid.append(smi)
    return valid


# ── PROTAC-DB parser ───────────────────────────────────────────────────────────

def remove_stereochemistry(smiles):
    '''
    Remove stereochemistry from a SMILES string.
    :param smiles: Input SMILES string
    :return: SMILES string without stereochemistry
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)

def parse_protacdb(csv_path: str) -> list[str]:
    """
    Extract linker SMILES from a PROTAC-DB CSV export.
    The CSV column is typically named 'Linker Smiles' or 'linker_smiles'.
    Download from: https://protacdb.cadd.zju.edu.cn/
    """
    if not HAS_RDKIT:
        sys.exit('RDKit required for PROTAC-DB parsing.')

    import pandas as pd

    df = pd.read_csv(csv_path)

    # Find linker column (case-insensitive)
    linker_col = None
    for col in df.columns:
        if 'smiles_r' in col.lower():
            linker_col = col
            break

    if linker_col is None:
        print(f'Available columns: {list(df.columns)}')
        sys.exit(
            'ERROR: Could not find a linker SMILES column in the CSV. '
            'Expected a column containing both "linker" and "smiles" in its name.'
        )

    print(f'Using column: {linker_col}  ({len(df)} rows)')

    raw = df[linker_col].dropna().astype(str).str.strip('/').unique().tolist()
    print(f'Unique raw linker SMILES: {len(raw)}')

    # Ensure attachment points use *
    cleaned = []
    for smi in raw:
        # PROTAC-DB sometimes uses [*] or [*:1] — normalise to *
        smi = smi.replace('[R1]', '*').replace('[R2]', '*').replace('[*]', '*')
        smi = remove_stereochemistry(smi)
        if not has_allowed_characters(smi, char[0]):
            print(f'  [SKIP] Invalid character in {smi}')
            continue
        cleaned.append(smi)

    return cleaned


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Prepare a PROTAC linker SMILES file for Link-INVENT Transfer Learning.'
    )
    parser.add_argument(
        '--protacdb', default=None,
        help='Path to a PROTAC-DB CSV export (optional). If not provided, '
             'a built-in curated set of ~40 linkers is used.'
    )
    parser.add_argument(
        '--out', default='protac_linkers.smi',
        help='Output SMILES file (default: protac_linkers.smi).'
    )
    parser.add_argument(
        '--no-validate', action='store_true',
        help='Skip RDKit validation of SMILES.'
    )
    args = parser.parse_args()

    if args.protacdb:
        print(f'Parsing PROTAC-DB export: {args.protacdb}')
        raw_linkers = parse_protacdb(args.protacdb)
    else:
        print('Using built-in curated PROTAC linker set.')
        raw_linkers = list(BUILTIN_LINKERS)

    print(f'Raw linkers collected : {len(raw_linkers)}')

    if not args.no_validate:
        valid_linkers = validate_linker_smiles(raw_linkers)
    else:
        valid_linkers = [s.strip() for s in raw_linkers if s.strip()]

    print(f'Valid linkers         : {len(valid_linkers)}')

    out_path = Path(args.out)
    out_path.write_text('\n'.join(valid_linkers) + '\n')

    print(f'\nWritten to: {out_path}')
    print(
        f'\nTo use in link_it.py, add:  --tl_dataset {out_path}\n'
        f'Example:\n'
        f'  python link_it.py --smiles_csv smiles.smi --dist_file input.txt \\\n'
        f'                    --tl_dataset {out_path} --tl_epochs 50\n'
    )


if __name__ == '__main__':
    main()