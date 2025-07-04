# PROTAC-tor: Automated PROTAC Generation Pipeline

Author: Jordan Harrison
Environment: Compute Canada (SLURM-based HPC)
Language: Python 3.11, C

---

## Overview

This repository automates the discovery of linker molecules for PROTAC complexes. It includes scripts that:

1. **Process protein-ligand complexes (PDB files)**
2. **Identify ligand pairs and calculate their atomic distances**
3. **Extract SMILES representations** of ligands using RDKit/OpenBabel
4. **Run Link-INVENT** to generate linkers between the fragments
5. **Submit the job to a GPU-enabled SLURM cluster**

---

## File Breakdown

### `lig_dist.py`

This script does the following:

* Iterates over all `.pdb` files in the working directory
* For each file:

  * Identifies all HETATM residue names (ligand IDs)
  * Filters files that contain **exactly 2 ligands**
  * Extracts 3D coordinates of both ligands
  * Calculates the **minimum inter-atomic distance** between the two ligands
  * Stores the closest complex
* Extracts the ligands from the closest complex into individual `.mol` files
* Converts them into **SMILES format**
* Writes the results to:

  * `smiles.csv` — the fragment pair (Kinase\_Ligand, E3\_Ligand)
  * `input.txt` — the min/max linker length bounds for Link-INVENT

This script also defines the `LinkInvent()` function. It:

* Loads SMILES and distance bounds
* Constructs a **JSON config** for Link-INVENT
* Writes a SLURM job script (`submit_linkinvent.sh`) with GPU resources
* Submits the job using `sbatch`

### `prodock.py`

This is the main engine it preprocesses files and submits the main job to the HPC

---


## How to Run the Pipeline on Compute Canada

### Step 1: Upload PDB's to PROTAC-tor directory

The PDB's Must:

* Each have a ligand
* Be prepared in a tool such as MOE
* Be correct PDB files with HETATM headers, and a 3 digit ligand name

### Step 2: Launch and Monitor SLURM Jobs

```bash
python scripts/prodock.py complex1.pdb complex2.pdb
```

This script:

* prepares the proteins
* Sets up the ZDOCK docking folder
* Prepares and submits a SLURM job (`run_docking.sh`) for the pair
* The SLURM script will:

  * Run ZDOCK
  * Extract ligand distances
  * Write SMILES and distance files
  * Launch `link.py` for Link-INVENT

```bash
squeue -u $USER
```

Use `cat zdock.out` or `linkinvent.out` inside complex subdirectories to check output.

---

## Output

Each complex subdirectory will contain:

* `zdock_result.out`: ZDOCK output
* `lig_distances.txt`: ligand distance metrics
* `smiles.csv`: SMILES of extracted ligands
* `linkinvent_output/`: Link-INVENT designed molecules

---

## Notes

* Ligand ID extraction assumes 3-letter residue names from `HETATM` records.
* OpenBabel sometimes fails on charged or malformed ligands.
* You can customize the linker scoring in `link.py` under the `scoring_function` field.

---

## Future Improvements

* Parallelization for ligand distance extraction
* More robust ligand detection (ignore ions or solvents)
* Integration with MOE for pre-processing
* Export top linkers to 3D PDBQT for downstream docking

---

Happy linking!

— Jordan

