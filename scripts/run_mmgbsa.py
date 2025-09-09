# Script to automatize the preparation of job files for MM/GBSA
# Author: Paola Vottero, Jordan Harrison

import os
from sympy import re
import shutil

# Load amber 25 command string
amber = 'module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber-pmemd/24.3 ambertools/25.0'

def pull_files(file_path_in, file_name, file_path_out):
    os.system(f'cp {file_path_in}/{file_name} {file_path_out}')

def process_trajectory(traj_file):
    os.system(f'cpptraj -i {traj_file} -o md.dcd')


def write_mmgbsa_job_file(md_dir, ligands):
    job_file = 'mmgbsa.job'
    with open(job_file, 'w') as fp:
        fp.write(f"""#!/bin/bash
#SBATCH --time=00-24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --job-name={md_dir}_mmgbsa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca
#SBATCH --account=def-aminpour

module --force purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 ambertools/25.0

cpptraj -i traj.in -o cpptraj.out

""")
        # MMGBSA for each ligand
        for ligand in ligands:
            fp.write(f"""$AMBERHOME/bin/MMPBSA.py -O -i mmgbsa.in -o out_{ligand}.dat \\
    -cp complex_stripped.prmtop \\
    -rp receptor.prmtop \\
    -lp {ligand}.prmtop \\
    -y md.dcd

""")
    return job_file

def main():
    DIR = os.path.dirname(os.path.abspath(__file__))
    complex_dirs = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('ternary')]

    for dir in complex_dirs:
        print(f'Processing directory: {dir}')
        md_dir = os.path.join(dir, "md")

        if not os.path.isdir(md_dir):
            print(f"Warning: {md_dir} does not exist, skipping...")
            continue

        # copy inputs into md/
        pull_files(DIR, 'traj.in', md_dir)
        pull_files(DIR, 'mmgbsa.in', md_dir)
        pull_files(os.path.join(dir, 'prep'), 'complex_stripped.prmtop', md_dir)
        pull_files(os.path.join(dir, 'prep'), 'receptor.prmtop', md_dir)

        # Detect ligand prmtop files
        ligand_files = []
        for i in range(1, 3):
            ligfile = os.path.join(dir, 'prep', f'ligand{i}.prmtop')
            if os.path.exists(ligfile):
                shutil.copy(ligfile, os.path.join(md_dir, f'ligand{i}.prmtop'))
                ligand_files.append(f'ligand{i}')
        # Fallback to single ligand
        if not ligand_files:
            ligfile = os.path.join(dir, 'prep', 'ligand.prmtop')
            if os.path.exists(ligfile):
                shutil.copy(ligfile, os.path.join(md_dir, 'ligand.prmtop'))
                ligand_files.append('ligand')

        cwd = os.getcwd()
        try:
            os.chdir(md_dir)
            job_file = write_mmgbsa_job_file(str(dir) + str(md_dir)[-1], ligand_files)
            os.system(f'sbatch {job_file}')
        finally:
            os.chdir(cwd)

main()