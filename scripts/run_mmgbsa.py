# Script to automatize the preparation of job files for MM/GBSA
# Author: Paola Vottero, Jordan Harrison

import os

# Load amber 25 command string
amber = 'module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber-pmemd/24.3 ambertools/25.0'

def pull_files(file_path_in, file_name, file_path_out):
    os.system(f'cp {file_path_in}/{file_name} {file_path_out}')

def process_trajectory(traj_file):
    os.system(f'cpptraj -i {traj_file} -o md.dcd')


def write_mmgbsa_job_file(dir):
    with open('mmgbsa.job', 'w') as fp:
        fp.write(f"""#!/bin/bash
    #SBATCH --time=00-12:00:00
    #SBATCH --mem-per-cpu=8G
    #SBATCH --cpus-per-task=1
    #SBATCH --job-name={dir}_mmgbsa
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=jaharri1@ualberta.ca
    #SBATCH --account=def-aminpour

    module load {amber}

    $AMBERHOME/bin/MMPBSA.py -O -i mmgbsa.in -o out1.dat -cp complex_stripped.prmtop -rp receptor.prmtop -lp ligand.prmtop -y md1.dcd
    """)
    
def main():
    DIR = os.path.dirname(os.path.abspath(__file__))
    complex_dirs = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith('ternary')]
    for dir in complex_dirs:
        print(f'Processing directory: {dir}')
        pull_files(DIR, 'traj.in', f'{dir}/md')
        pull_files(DIR, 'mmgbsa.in', f'{dir}/md')
        pull_files(f'{dir}/prep', 'complex_stripped.prmtop', f'{dir}/md')
        pull_files(f'{dir}/prep', 'receptor.prmtop', f'{dir}/md')
        pull_files(f'{dir}/prep', 'ligand.prmtop', f'{dir}/md')

        os.chdir(f'{dir}/md')

        process_trajectory('traj.in')

        write_mmgbsa_job_file(DIR)

        os.system('sbatch mmgbsa.job')

main()








    


