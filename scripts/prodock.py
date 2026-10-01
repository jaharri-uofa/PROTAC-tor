'''
Protein Docking Automation Script (multi-config mode), prepares protein complexes and establishes directorys for the rest of the run
Author: Jordan Harrison
'''

#!/usr/bin/env python3
import os
import shutil
from pathlib import Path
import stat
import subprocess
from rdkit import Chem

def remove_stereochemistry(smiles):
    '''
    Removes stereochemistry from a smiles string
    :param smiles: smiles string
    :return: smiles without stereochemistry
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)

print("Starting ZDOCK docking automation (multi-config mode)...")

# get directory paths
base_dir = Path.cwd()
zdock_dir = base_dir / "ZDOCK"
scripts_dir = base_dir / "scripts"
shell_dir = base_dir / "shell"
protein_complexes_dir = base_dir / "complexes"
protein_complexes_dir.mkdir(exist_ok=True)

# includes all the scripts in the github
zdock = ["zdock", "create_lig", "create.pl", "mark_sur", "uniCHARMM", 'linkinvent.prior', 'traj.in', 'mmgbsa.in', 'rmsd.in', 'protac_linkers.smi']
python = ["lig_dist.py", "prodock.py", "link_it.py", "dock.py", "analysis.py", "md.py", "md_mmgbsa.py", "control_md.py", "run_mmgbsa.py", "md_analysis.py", "benchmark.py"]
shell = ["driver.sh", "link_it.sh", "prodock.sh"]
bench = ["benchmark.py"]

# Read config.txt with multiple blocks
config_path = base_dir / "config.txt"
if not config_path.exists():
    print("ERROR: config.txt not found!")
    exit(1)

with open(config_path, "r") as f:
    config_content = f.read()

# Split config into blocks by '////'
blocks = [block.strip() for block in config_content.split('////') if block.strip()]

for block in blocks:
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    benchmark_mode = False
    if lines and lines[0].lower() == "benchmark":
        benchmark_mode = True
        lines = lines[1:]

    if len(lines) != 4:
        print("ERROR: Each config block must have 4 lines: e3_ligase.pdb, poi.pdb, e3ligand_smiles, poi_ligand_smiles")
        if benchmark_mode:
            print("       Or benchmark block can start with 'benchmark' followed by the same 4 lines.")
        continue

    ligase_pdb, poi_pdb, lig1_smiles, lig2_smiles = lines

    # Remove stereochemistry
    lig1_smiles = remove_stereochemistry(lig1_smiles)
    lig2_smiles = remove_stereochemistry(lig2_smiles)

    # Use POI PDB stem as subdirectory name
    complex_name = f"{Path(poi_pdb).stem}"
    complex_dir = protein_complexes_dir / complex_name
    complex_dir.mkdir(parents=True, exist_ok=True)

    # Copy required ZDOCK files
    for filename in zdock:
        src = zdock_dir / filename
        dest = complex_dir / filename
        shutil.copy(src, dest)
        dest.chmod(dest.stat().st_mode | stat.S_IEXEC)

    # Copy scripts
    for script_name in python:
        shutil.copy(scripts_dir / script_name, complex_dir)

    for script_name in shell:
        shutil.copy(shell_dir / script_name, complex_dir)

    if benchmark_mode:
        for script_name in bench:
            shutil.copy(scripts_dir / script_name, complex_dir)

    # Copy receptor and ligand PDBs
    shutil.copy(poi_pdb, complex_dir / "receptor.pdb")
    print(f"Copied {poi_pdb} to {complex_dir / 'receptor.pdb'}")
    shutil.copy(ligase_pdb, complex_dir / "ligand.pdb")
    print(f"Copied {ligase_pdb} to {complex_dir / 'ligand.pdb'}")

    # Add dummy SEQRES
    (complex_dir / "SEQRES").write_text("DUMMYSEQRES\n")

    # Write benchmark mode marker if requested
    if benchmark_mode:
        (complex_dir / "BENCHMARK_MODE").write_text("benchmark\n")

    # Write smiles.smi
    smiles_path = complex_dir / "smiles.smi"
    with open(smiles_path, "w") as f:
        f.write(f"{lig1_smiles}|{lig2_smiles}\n")
    '''
    # Execute SLURM script
    slurm_script_path = complex_dir / "prodock.sh"
    os.chmod(slurm_script_path, 0o755)

    print(f"Submitting protein docking job for {complex_name}...")
    subprocess.run(["sbatch", "prodock.sh"], cwd=complex_dir)
    print(f"Job for {complex_name} submitted.")
    '''