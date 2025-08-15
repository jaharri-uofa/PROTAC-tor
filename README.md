# PROTACtor

A pipeline for automated PROTAC docking, linker design, and molecular dynamics analysis.

## Project Structure

```
PROTAC-tor/
    protactor.sh           # Main driver script for the full pipeline (SLURM)
    README.md
    scripts/
        analysis.py        # Final analysis of results
        dock.py            # Post-processing of docking results
        lig_dist.py        # Ligand distance calculation
        link_it.py         # Linker design with REINVENT
        md_mmgbsa.py       # MMGBSA calculation after MD
        md.py              # MD preparation and filtering
        prodock.py         # Docking setup and job submission (config-based)
    ZDOCK/
        zdock, create_lig, create.pl, mark_sur, uniCHARMM, linkinvent.prior, libg2c.so.0
    shell/
        driver.sh          # Main driver in complex diretory
        prodock.sh         # Shell script to launch protein-protein docking
        link_it.sh         # Shell script to launch Link-Invent
```

---

## Quick Start

### 1. Prepare Your Input

Create a `config.txt` file in your project root with **5 lines**:

```
////
e3_ligase.pdb
poi.pdb
E3_LIGAND_SMILES
POI_LIGAND_SMILES
```

- `e3_ligase.pdb`: Path to your E3 ligase PDB file
- `poi.pdb`: Path to your protein of interest PDB file
- `E3_LIGAND_SMILES`: SMILES string for the E3 ligand
- `POI_LIGAND_SMILES`: SMILES string for the POI ligand

- Multiple inputs are also supported, simply duplicate the text above ie:
```
////
e3_ligase1.pdb
poi1.pdb
E3_LIGAND_SMILES1
POI_LIGAND_SMILES1
////
e3_ligase2.pdb
poi2.pdb
E3_LIGAND_SMILES2
POI_LIGAND_SMILES2
```
---

### 2. Run the Pipeline

Submit the main driver script to SLURM:

```bash
sbatch protactor.sh
```

This script will:
- Run all pipeline steps in order (`prodock.py`, `lig_dist.py`, `link_it.py`, `dock.py`, `md.py`, `md_mmgbsa.py`, `analysis.py`)
- Wait for all jobs to finish between steps

---

## Pipeline Steps

1. **prodock.py**  
   Reads `config.txt`, sets up the docking complex, copies required files, and submits the initial docking job.

2. **lig_dist.py**  
   Calculates ligand distances for further filtering.

3. **link_it.py**  
   Designs linkers using REINVENT.

4. **dock.py**  
   Post-processes docking results, extracts top poses, and prepares for MD.

5. **md.py**  
   Prepares molecular dynamics input files and filters top complexes.

6. **md_mmgbsa.py**  
   Runs molecular dynamics and prepares for MM/GBSA.

7. **analysis.py**  
   Final analysis and summary of results.

---

## Notes

- All python scripts are located in the `scripts/` directory.
- All shell scripts are located in the `shell/` directory.
- The pipeline is designed to be run from the project root.
- Make sure all dependencies (RDKit, OpenBabel, gnina, etc.) are available in your environment or loaded via modules as in the SLURM scripts.

---

## Contact

Author: jaharri-uofa, SirTurtle, paolav8ero
Email: jaharri1@ualberta.ca

---

Happy PROTACting!
