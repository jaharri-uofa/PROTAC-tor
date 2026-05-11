'''
link_it.py
Builds and submits a Link-INVENT staged-learning job using REINVENT 4.
Optionally fine-tunes the prior via Transfer Learning on a PROTAC linker
dataset before running RL, which biases generation toward PROTAC-like
linker chemistry.

New flag:
  --tl_dataset  Path to a SMILES file of known PROTAC linkers (one per line,
                with `*` attachment points, e.g. *CC(=O)NCCOCCN*).
                When supplied, a TL stage runs first to fine-tune the prior,
                and the RL stage uses the resulting model as its agent.

Author: Jordan Harrison
'''

import argparse
import pandas as pd
import subprocess
import os
import toml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops
import numpy as np

# ── Molecular feature extraction ───────────────────────────────────────────────

def molecule_features(smiles: str) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolecularWeight":  rdMolDescriptors.CalcExactMolWt(mol),
        "TPSA":             rdMolDescriptors.CalcTPSA(mol),
        "HBondAcceptors":   rdMolDescriptors.CalcNumHBA(mol),
        "HBondDonors":      rdMolDescriptors.CalcNumHBD(mol),
        "NumRotBond":       rdMolDescriptors.CalcNumRotatableBonds(mol),
        "NumRings":         rdMolDescriptors.CalcNumRings(mol),
        "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "SlogP":            rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
    }


def extract_warhead_smiles(smiles: str) -> tuple[str, str]:
    parts = smiles.split('|')
    if len(parts) != 2:
        raise ValueError(f"Invalid PROTAC SMILES format (expected 'wh1|wh2'): {smiles}")
    return parts[0].strip(), parts[1].strip()


def longest_path_length(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return int(rdmolops.GetDistanceMatrix(mol).max())

def prepare_linkinvent_tl_smiles(raw_smi_file: str,
                                  out_smi_file: str,
                                  warhead1: str,
                                  warhead2: str) -> int:
    """
    Convert a single-column linker SMILES file into the 3-column TSV format
    that REINVENT4 Link-INVENT TL requires:
        warhead1 <TAB> linker <TAB> warhead2

    Lines that are blank or fail RDKit sanitization are skipped.
    Returns the number of valid linkers written.
    """
    from rdkit import Chem
    written = 0
    with open(raw_smi_file) as fin, open(out_smi_file, 'w') as fout:
        for line in fin:
            smi = line.strip().split()[0]  # ignore any trailing comments/names
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"  [WARN] skipping invalid SMILES: {smi}")
                continue
            fout.write(f"{warhead1}\t{smi}\t{warhead2}\n")
            written += 1
    print(f"  Reformatted {written} linkers → {out_smi_file}")
    return written


# ── Transfer Learning TOML ─────────────────────────────────────────────────────

def generate_tl_toml(smiles_csv: str, tl_dataset: str, prior_file: str,
                     output_toml: str, finetuned_prior: str,
                     num_epochs: int = 50):
    """
    Generate a REINVENT4 Transfer Learning TOML for Link-INVENT.

    The TL run fine-tunes `prior_file` on the PROTAC linker SMILES in
    `tl_dataset` and saves the result to `finetuned_prior`.

    Dataset format (one linker per line, attachment points marked with *):
        *CC(=O)NCCOCCN*
        *c1ccncc1CC(=O)NCCN*
        ...

    REINVENT4 TL reference:
        https://github.com/MolecularAI/REINVENT4/tree/main/configs
    Note: parameter names follow REINVENT4 ≥ 4.4. If you are on an older
    version, check configs/transfer_learning/linkinvent.toml in the repo.
    """
    # Reformat single-column linker SMILES → 3-column TSV required by Link-INVENT TL
    with open(smiles_csv) as f:
        smiles1, smiles2 = extract_warhead_smiles(f.readline().strip())

    formatted_dataset = tl_dataset.replace('.smi', '_linkinvent_fmt.smi')
    n = prepare_linkinvent_tl_smiles(tl_dataset, formatted_dataset, smiles1, smiles2)
    if n == 0:
        raise ValueError(f"No valid linker SMILES found in {tl_dataset}. Aborting.")

    config = {
        "run_type":       "transfer_learning",
        "device":         "cuda:0",
        "tb_logdir":      "tb_logs_tl",
        "json_out_config": "tl_out_config.json",
        "parameters": {
            # Input model (the base linkinvent prior)
            "prior_file":        prior_file,
            # Where the fine-tuned model is saved
            "agent_file":        finetuned_prior,
            # Training SMILES — one linker per line with * attachment points
            "smiles_file":       tl_dataset,
            "batch_size":        64,
            "sample_batch_size": 128,
            "num_epochs":        num_epochs,
            # Save a checkpoint every N epochs so you can resume if interrupted
            "save_every_n_epochs": max(1, num_epochs // 5),
            # Clip gradient norm — prevents exploding gradients during fine-tuning
            "clip_gradient_norm": 1.0,
        },
    }

    with open(output_toml, 'w') as f:
        toml.dump(config, f)

    print(f"TL TOML written to {output_toml}")
    print(f"  Prior        : {prior_file}")
    print(f"  Dataset      : {tl_dataset}")
    print(f"  Fine-tuned → : {finetuned_prior}")
    print(f"  Epochs       : {num_epochs}")


# ── RL (staged learning) TOML ──────────────────────────────────────────────────

def generate_rl_toml(smiles_csv: str, dist_file: str,
                     output_toml: str,
                     prior_file: str = 'linkinvent.prior',
                     agent_file: str | None = None):
    """
    Generates the staged RL TOML for Link-INVENT.
    If `agent_file` is provided (i.e. a fine-tuned prior from TL), the RL
    run starts from that checkpoint instead of the raw prior.
    """
    with open(dist_file) as f:
        min_dist, max_dist = map(float, f.readline().strip().split(','))

    with open(smiles_csv) as f:
        smiles1, smiles2 = extract_warhead_smiles(f.readline().strip())

    chem_data = [molecule_features(s) for s in [smiles1, smiles2]
                 if molecule_features(s) is not None]

    carb = 1.3336   # projected C-C bond length (Å)
    peg  = 44.01 / 2  # PEG repeat unit mass / 2 (approximate linker mass per bond)

    base_stage = {
        "termination":  "simple",
        "chkpt_file":   "",
        "max_score":    0.6,
        "min_steps":    1000,
        "max_steps":    10000,
        "scoring": {
            "type": "geometric_mean",
            "component": [
                {"FragmentMolecularWeight": {"endpoint": [{
                    "name": "Linker MW",
                    "weight": 1,
                    "transform": {
                        "type": "reverse_sigmoid",
                        "high": (int(max_dist) / carb) * peg,
                        "low":  (int(min_dist) / carb) * peg,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentGraphLength": {"endpoint": [{
                    "name": "Graph length",
                    "weight": 1,
                    "transform": {
                        "type": "sigmoid",
                        "high": int(max_dist) / carb,
                        "low":  int(min_dist) / carb,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentEffectiveLength": {"endpoint": [{
                    "name": "Effective length (Å)",
                    "weight": 1,
                    "transform": {
                        "type": "sigmoid",
                        "high": int(max_dist),
                        "low":  int(min_dist),
                        "k": 0.5,
                    },
                }]}},
                {"FragmentLengthRatio": {"endpoint": [{
                    "name": "Length ratio",
                    "weight": 1,
                    "transform": {
                        "type": "sigmoid",
                        "high": 1.0,
                        "low":  0.99,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentTPSA": {"endpoint": [{
                    "name": "TPSA",
                    "weight": 1,
                    "transform": {
                        "type": "sigmoid",
                        "high": 80.0,
                        "low":  40.0,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentSlogP": {"endpoint": [{
                    "name": "SlogP",
                    "weight": 1,
                    "transform": {
                        "type": "reverse_sigmoid",
                        "high": 5.0,
                        "low":  1.0,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentHBondAcceptors": {"endpoint": [{
                    "name": "HBA",
                    "weight": 1,
                    "transform": {
                        "type": "reverse_sigmoid",
                        "high": 5.0,
                        "low":  1.0,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentHBondDonors": {"endpoint": [{
                    "name": "HBD",
                    "weight": 1,
                    "transform": {
                        "type": "reverse_sigmoid",
                        "high": 1.0,
                        "low":  0.0,
                        "k": 0.5,
                    },
                }]}},
                {"FragmentNumRotBond": {"endpoint": [{
                    "name": "Rotatable bonds",
                    "weight": 1,
                    "transform": {
                        "type": "sigmoid",
                        "high": int(max_dist) / carb,
                        "low":  int(min_dist) / carb,
                        "k": 0.5,
                    },
                }]}},
                # SA score downweights synthetically inaccessible linkers
                {"SAScore": {"endpoint": [{
                    "name": "SA score",
                    "weight": 3,
                }]}},
            ],
        },
    }

    stages = []
    for i in range(3):
        stage = dict(base_stage)
        stage["chkpt_file"] = f"stage{i + 1}.chkpt"
        stage["min_steps"]  = 1000 + i * 1000
        stage["max_steps"]  = 5000 + i * 2000
        stages.append(stage)

    # If a fine-tuned prior was produced by TL, use it as the starting agent.
    # This means RL exploration begins from a PROTAC-biased distribution.
    effective_agent = agent_file if agent_file else prior_file

    config = {
        "run_type":       "staged_learning",
        "device":         "cuda:0",
        "tb_logdir":      "tb_logs",
        "json_out_config": "staged_linkinvent.json",
        "parameters": {
            "summary_csv_prefix": "linkinvent_stage",
            "use_checkpoint":     False,
            "purge_memories":     False,
            "prior_file":         prior_file,         # always the base prior
            "agent_file":         effective_agent,    # TL output (or base prior)
            "smiles_file":        smiles_csv,
            "batch_size":         64,
            "unique_sequences":   True,
            "randomize_smiles":   True,
            "tb_isim":            False,
        },
        "learning_strategy": {
            "type":  "dap",
            "sigma": 256,
            "rate":  0.0001,
        },
        "diversity_filter": {
            "type":               "ScaffoldSimilarity",
            "bucket_size":        100,
            "minscore":           0.4,
            "minsimilarity":      0.3,
            "penalty_multiplier": 1.0,
        },
        "stage": stages,
    }

    with open(output_toml, 'w') as f:
        toml.dump(config, f)

    print(f"RL TOML written to {output_toml}")
    if agent_file:
        print(f"  Agent (fine-tuned prior) : {agent_file}")
    else:
        print(f"  Agent (base prior)       : {prior_file}")


# ── SLURM scripts ──────────────────────────────────────────────────────────────

MODULE_BLOCK = """module --force purge
module load StdEnv/2023
module load openbabel/3.1.1 gcc/12.3 cmake cuda/12.6
module load python/3.11.5 scipy-stack/2025a rdkit/2024.09.6 python-build-bundle/2025b

source ~/reinvent4/bin/activate
export PATH=$HOME/.local/bin:$PATH
"""


def write_tl_slurm(tl_toml: str, rl_toml: str,
                   slurm_script: str = 'submit_tl_then_rl.sh'):
    """
    Write a SLURM script that:
      1. Runs Transfer Learning to fine-tune the prior.
      2. On success, chains the RL staged-learning job as a dependency.
    """
    contents = f"""#!/bin/bash
#SBATCH --job-name=linkinvent_tl
#SBATCH --output=linkinvent_tl.out
#SBATCH --error=linkinvent_tl.err
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=0-04:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

{MODULE_BLOCK}

echo "=== Step 1: Transfer Learning (fine-tuning prior on PROTAC linkers) ==="
reinvent -l tl.log {tl_toml}

TL_EXIT=$?
if [ $TL_EXIT -ne 0 ]; then
    echo "ERROR: Transfer Learning failed (exit $TL_EXIT). Aborting RL stage."
    exit $TL_EXIT
fi

echo "TL complete. Fine-tuned prior ready."

echo "=== Step 2: Submitting RL staged-learning job ==="
TL_JOBID=$SLURM_JOB_ID
rl_jobid=$(sbatch --parsable \\
    --dependency=afterok:$TL_JOBID \\
    --job-name=linkinvent_rl \\
    --output=linkinvent_rl.out \\
    --error=linkinvent_rl.err \\
    --gres=gpu:1 --mem=4G --cpus-per-task=1 --time=0-02:30 \\
    --account=def-aminpour \\
    --wrap="{MODULE_BLOCK.strip()} && reinvent -l staged.log {rl_toml}")
echo "Submitted RL job as $rl_jobid"

# Chain dock.py after RL finishes
dock_jobid=$(sbatch --parsable \\
    --dependency=afterok:$rl_jobid \\
    --mem=2G --job-name=dockpy --output=dockpy.out --error=dockpy.err \\
    --wrap="module load StdEnv/2023 python/3.11 scipy-stack/2025a rdkit/2024.09.6 \\
            openbabel/3.1.1 gcc/12.3 cmake cuda/12.2 python-build-bundle/2025b \\
            openmpi/4.1.5 ambertools/25.0; python dock.py")
echo "Submitted dock.py as job $dock_jobid"
"""
    with open(slurm_script, 'w') as f:
        f.write(contents)
    print(f"TL+RL SLURM script written to {slurm_script}")


def write_rl_only_slurm(rl_toml: str,
                        slurm_script: str = 'submit_linkinvent.sh'):
    """Original RL-only SLURM script (no TL stage)."""
    contents = f"""#!/bin/bash
#SBATCH --job-name=linkinvent_gpu
#SBATCH --output=linkinvent.out
#SBATCH --error=linkinvent.err
#SBATCH --gres=gpu:1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:30
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca

{MODULE_BLOCK}

echo "Running REINVENT Link-INVENT staged learning..."
reinvent -l staged.log {rl_toml}

echo "Exit code: $?"

link_jobid=$SLURM_JOB_ID

dock_jobid=$(sbatch --parsable \\
    --dependency=afterok:$link_jobid \\
    --mem=2G --job-name=dockpy --output=dockpy.out --error=dockpy.err \\
    --wrap="module load StdEnv/2023 python/3.11 scipy-stack/2025a rdkit/2024.09.6 \\
            openbabel/3.1.1 gcc/12.3 cmake cuda/12.2 python-build-bundle/2025b \\
            openmpi/4.1.5 ambertools/25.0; python dock.py")
echo "Submitted dock.py as job $dock_jobid"
"""
    with open(slurm_script, 'w') as f:
        f.write(contents)
    print(f"RL SLURM script written to {slurm_script}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build and submit a Link-INVENT job (with optional TL fine-tuning).'
    )
    parser.add_argument('--smiles_csv',   required=True,
                        help='Warhead SMILES file (warhead1|warhead2 format).')
    parser.add_argument('--dist_file',    required=True,
                        help='Distance file: single line "min,max" in Å.')
    parser.add_argument('--output_toml',  default='sampling.toml',
                        help='Output RL TOML filename.')
    parser.add_argument('--slurm_script', default='submit_linkinvent.sh',
                        help='Output SLURM script filename.')

    # ── Transfer Learning options ──────────────────────────────────────────
    parser.add_argument(
        '--tl_dataset',
        default=None,
        help=(
            'Path to a SMILES file of known PROTAC linkers for Transfer Learning. '
            'One SMILES per line with * attachment points, e.g. *CC(=O)NCCOCCN*. '
            'When provided, a TL stage fine-tunes the prior before RL begins, '
            'biasing generation toward PROTAC-like linker chemistry. '
            'See prepare_tl_dataset.py to build this file from PROTAC-DB.'
        ),
    )
    parser.add_argument(
        '--tl_epochs',
        type=int,
        default=50,
        help='Number of TL training epochs (default: 50).',
    )
    parser.add_argument(
        '--tl_toml',
        default='tl_config.toml',
        help='Output TL TOML filename (default: tl_config.toml).',
    )
    parser.add_argument(
        '--finetuned_prior',
        default='linkinvent_finetuned.prior',
        help='Filename for the TL output (fine-tuned prior).',
    )
    parser.add_argument(
        '--prior_file',
        default='linkinvent.prior',
        help='Base Link-INVENT prior (default: linkinvent.prior).',
    )
    args = parser.parse_args()

    assert os.path.exists(args.smiles_csv), f'Missing: {args.smiles_csv}'
    assert os.path.exists(args.dist_file),  f'Missing: {args.dist_file}'

    use_tl = args.tl_dataset is not None

    if use_tl:
        assert os.path.exists(args.tl_dataset), f'Missing TL dataset: {args.tl_dataset}'
        print('\n=== Transfer Learning enabled ===')
        print(f'Dataset  : {args.tl_dataset}')
        print(f'Epochs   : {args.tl_epochs}')
        print(f'Output   : {args.finetuned_prior}')
        print()

        # 1. Write TL TOML
        generate_tl_toml(
            smiles_csv=args.smiles_csv,
            tl_dataset=args.tl_dataset,
            prior_file=args.prior_file,
            output_toml=args.tl_toml,
            finetuned_prior=args.finetuned_prior,
            num_epochs=args.tl_epochs,
        )

        # 2. Write RL TOML (agent = fine-tuned prior)
        generate_rl_toml(
            smiles_csv=args.smiles_csv,
            dist_file=args.dist_file,
            output_toml=args.output_toml,
            prior_file=args.prior_file,
            agent_file=args.finetuned_prior,  # ← key: use TL output
        )

        # 3. Write combined TL → RL SLURM script
        write_tl_slurm(
            tl_toml=args.tl_toml,
            rl_toml=args.output_toml,
            slurm_script=args.slurm_script,
        )
        print(f'\nRun:  sbatch {args.slurm_script}')
        print('This will: (1) fine-tune the prior via TL, then (2) run RL from the fine-tuned model.')

    else:
        print('\n=== RL only (no TL fine-tuning) ===')
        print('Tip: supply --tl_dataset to fine-tune on known PROTAC linkers first.')
        print()

        generate_rl_toml(
            smiles_csv=args.smiles_csv,
            dist_file=args.dist_file,
            output_toml=args.output_toml,
            prior_file=args.prior_file,
            agent_file=None,
        )

        write_rl_only_slurm(
            rl_toml=args.output_toml,
            slurm_script=args.slurm_script,
        )
        print(f'\nRun:  sbatch {args.slurm_script}')


if __name__ == '__main__':
    main()
