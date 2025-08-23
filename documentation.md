# PROTAC-tor Computational Workflow Documentation

This document details the computational methods and workflows implemented in the **PROTAC-tor** pipeline. Each script is explained with emphasis on the scientific methodology, data flow, and analytical purpose. The intended audience is researchers working with PROTAC docking, linker analysis, and binding free energy calculations.

---

## 1. `prodock.py`

### Purpose
Automates protein–protein docking using the ZDOCK algorithm. Prepares docking jobs, submits them to a high-performance cluster (HPC) via SLURM, and organizes results for downstream analysis.

### Methods
- Identifies protein complex PDB structures
- Generates docking jobs using ZDOCK
- Submits jobs with appropriate resource requests
- Ensures reproducibility via standardized SLURM scripts

### Functions
1. Scan project directories for complexes (`complexX_complexY`)
2. Construct input paths for receptor and ligand
3. Write SLURM batch script invoking ZDOCK
4. Submit job with `sbatch`

### Input
- Two PDB files: receptor and ligand
- Directory structure:

```
Proteincomplexes/
└── complex1_complex2/
├── receptor.pdb
├── ligand.pdb
```

### Output
- Docking job scripts (`dock.job`) per complex
- ZDOCK raw output files (poses, scores)
- SLURM submission logs

---

## 2. `lig_dist.py`

### Purpose
Evaluates geometric feasibility of linker connection between docked PROTAC components. Computes anchor–anchor distances and filters implausible ternary poses.

### Methods
- Parses PDB docking results
- Identifies linker attachment points (anchors)
- Computes Euclidean distances between anchor atoms
- Filters complexes according to min/max linker thresholds

### Functions
- Distance calculation using Euclidean norm (`sqrt(x²+y²+z²)`)
- Filtering complexes outside the distance window

### Input
- Docked complex PDB files
- User-defined linker distance thresholds

### Filtering
- Complexes are first sorted by the anchor atom distances (smallest to biggest)
- Top N ternary complexes based on protein interactions are then passed as output

### Output
- Filtered list of feasible PROTAC complexes
- Tabular summary:
  - Complex ID
  - Distance between anchors
  - Pass/fail status

---

## 3. `link_it.py` – Automated Linker Design with Link-INVENT

### Overview
Automates de novo linker design using Link-INVENT (REINVENT 4.0). Integrates molecular feature extraction, distance constraints, and multi-stage reinforcement learning (RL). Generates TOML configuration files, SLURM scripts, and submits GPU jobs on HPC.

### Computational Methods

#### 1. Molecular Feature Extraction
Uses RDKit to extract descriptors:

- MolecularWeight
- TPSA
- HBondAcceptors / HBondDonors
- NumRotBond, NumRings, NumAromaticRings
- SlogP
- Longest path length

#### 2. Warhead Handling
- PROTACs represented as `warhead1|warhead2`
- Function `extract_warhead_smiles()` parses and validates input

#### 3. Distance Constraints and Scoring
- Reads min/max linker lengths
- Scoring components:
  - FragmentGraphLength
  - FragmentEffectiveLength
  - FragmentLengthRatio
  - FragmentTPSA, FragmentHBondAcceptors/Donors
  - SlogP, SAScore, Ring Count
- Sigmoid transforms normalize rewards

#### 4. Staged Learning Configuration
- Multi-stage TOML (3 stages)
- Checkpoints allow resuming
- RL strategy: Diversity-augmented policy (DAP)
- Scaffold similarity filter for chemical diversity

#### 5. HPC Integration
- Generates SLURM scripts for GPU jobs (1 GPU, 1 CPU, 4GB RAM, 2.5h runtime)
- Loads required modules
- Executes `reinvent -l <logfile> <TOML>`
- `submit_job()` submits jobs

### Workflow Summary
**Input:** SMILES CSV + distance file  
**Processing:** Feature extraction → geometric constraints → TOML → SLURM  
**Output:** `sampling.toml`, `submit_linkinvent.sh`, logs, checkpoints

---

## 4. `dock.py` – Automated PROTAC Docking with gnina

### Overview
Automates ternary complex docking using PRosettaC-inspired workflows and Gnina.

### Methods

#### 1. Ligand Preparation
- Converts SMILES → 3D SDF
- Adds hydrogens
- ETKDG embedding
- UFF optimization

#### 2. Warhead Parsing
- Enforces `warhead1|warhead2` format

#### 3. Protein Preprocessing
- `remove_ligand()` generates:
  - Ligand-stripped PDB (`*_nolig.pdb`)
  - Ligand-only PDB (`*_lig.pdb`)

#### 4. PROTAC Selection
- `get_PROTAC()` selects top-N candidates based off the linkinvent score
- Writes `top_smiles.txt` and `protac.sdf`

#### 5. Docking Configuration
- Job directory per complex
- Gnina config: receptor, ligand, autobox, num_modes=5, exhaustiveness=32

#### 6. Job Automation
- Writes `job.sh`
- Resources: 16 CPUs, 128 MB/CPU, 2.5h
- Executes:

gnina --config config


---

## 5. `md.py` – Post-Docking MD Preparation

### Overview
Automates ternary PROTAC complex preparation for MD simulations.

### Methods
- Extracts ligands from `docked.sdf.gz`
- Converts SDF → PDB
- Combines receptor + ligand
- Filters top complexes by docking affinity
- Manages directories
- Calls `md_mmgbsa.py`

### Key Functions
- `add_ligand()` combines receptor + ligand
- `distance()` & `lys_dist()` compute distances
- `extract_sdf_gz()`, `sdf_to_smiles_affinity()`
- `get_main_ligand_id()` identifies primary ligand
- Saves top 5 unique complexes → `top5_complexes.csv`

### Output Structure

```
ternary_complex1/
├─ complex1.pdb
├─ receptor1.pdb
├─ ligand1.pdb
├─ ligand_resname.txt
├─ ligand.sdf
├─ ternary.pdb
```

---

## 6. `md_mmgbsa.py` – MD Setup & GPU Submission

### Overview
Prepares ternary PROTAC complexes for Amber MD simulations on GPUs.

### Methods
- Input: complex, receptor, ligand PDB + ligand residue name
- Directory setup: `prep/` and `md/`
- Hydrogen manipulation (`pdb4amber`, `reduce`)
- Ligand parameterization (GAFF)
- Receptor & complex parameterization (solvation, ions)
- MD input files: `01-min1.in` → `06-prod.in`
- SLURM GPU job: `run_md.job`

### Output Structure
```
md/
├─ complex.prmtop
├─ complex.inpcrd
├─ receptor.prmtop
├─ receptor.inpcrd
├─ ligand.prmtop
├─ ligand.inpcrd
├─ 01-min1.in
├─ 02-min2.in
├─ 03-min3.in
├─ 04-heat.in
├─ 05-npt.in
├─ 06-prod.in
├─ run_md.job
└─ [MD output files]
```

---

## 7. `run_mmgbsa.py` – Automated MM/GBSA Job Preparation

### Overview
Automates MM/GBSA calculations for PROTAC ternary complexes.

### Methods
- Loads AmberTools 25:

module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber-pmemd/24.3 ambertools/25.0


- Searches `ternary*` directories → `md/` subdirs
- Copies inputs: `traj.in`, `mmgbsa.in`, `complex_stripped.prmtop`, `receptor.prmtop`, `ligand.prmtop`
- Preprocess trajectories:

cpptraj -i traj.in -o md.dcd


- Generates SLURM job `mmgbsa.job` and runs:

$AMBERHOME/bin/MMPBSA.py -O -i mmgbsa.in -o out1.dat
-cp complex_stripped.prmtop
-rp receptor.prmtop
-lp ligand.prmtop
-y md.dcd


- Submits jobs with `sbatch`

---

## 8. `analysis.py` – PROTACtor Output Analysis

### Overview
Automates analysis of docking, linker generation, and MD/MM-GBSA outputs. Writes `output.txt`.

### Methods
- Protein-Protein Docking: `pp_compatibility()`
- Distance Analysis: `min_max(path)`
- Lysine & Warhead Extraction: `get_lysines()`, `get_warhead_smiles()`
- Linker Metrics: `get_total_linkers()`, `get_max_linker_score()`
- PROTAC Docking: `_extract_affinities()`, `get_warheads_binding_affinity()`, `display_top_results()`, `get_highest_protac_binding_affinity()`
- MD/MM-GBSA placeholders: `get_trajectory_rmsd()`, `get_lysine_accessibility_score()`, `get_mmgbsa_scores()`

### Output
- Aggregates metrics into `output.txt`:
  - Docking compatibility
  - Min/Max distances
  - Warhead SMILES
  - Linker stats
  - Top PROTAC affinities
  - MD/MM-GBSA placeholders

---
