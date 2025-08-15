# Updated Script with GPU Bash Script for Running MD Simulations
# Author: Paola Vottero
# Editor: Jordan Harrison
import os
import sys
import math

amber = 'module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber/22.5-23.5'

# Ensure three arguments are passed
if len(sys.argv) != 4 or not sys.argv[1].endswith('.pdb') or not sys.argv[2].endswith('.pdb') or not sys.argv[3].endswith('.pdb'):
    sys.exit('Please provide three PDB files: the complex, receptor, and ligand.')

# Input paths
complex_path = sys.argv[1]
receptor_path = sys.argv[2]
ligand_path = sys.argv[3]

# Create directories
complex_name = os.path.basename(complex_path).replace('.pdb', '')
ligand_filename = os.path.basename(ligand_path)
os.makedirs(f'{complex_name}/prep', exist_ok=True)
os.makedirs(f'{complex_name}/md', exist_ok=True)

# Copy files to prep directory
os.system(f'cp {complex_path} {complex_name}/prep/')
os.system(f'cp {receptor_path} {complex_name}/prep/')
os.system(f'cp {ligand_path} {complex_name}/prep/')

# Move to prep directory
os.chdir(f'{complex_name}/prep')

# Get filenames
complex_filename = os.path.basename(complex_path)
receptor_filename = os.path.basename(receptor_path)

# Remove hydrogens
os.system(f'pdb4amber --nohyd -i {complex_filename} -o complex_noh.pdb')
os.system(f'pdb4amber --nohyd -i {receptor_filename} -o receptor_noh.pdb')
os.system(f'pdb4amber --nohyd -i {ligand_filename} -o ligand_noh.pdb')

# Add hydrogens
os.system('reduce complex_noh.pdb>complex_h.pdb')
os.system('reduce receptor_noh.pdb>receptor_h.pdb')
os.system('reduce ligand_noh.pdb>ligand_h.pdb')


# Generate ligand parameters
print("Generating ligand parameters...")
try:
    print(f"Running antechamber to generate ligand files for {ligand_filename}...")
    os.system(f"antechamber -i ligand_h.pdb -fi pdb -o ligand.mol2 -fo mol2 -at gaff")
    os.system(f"parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod")
    print(f"Ligand files generated successfully: {ligand_filename[:-4]}.mol2 and {ligand_filename[:-4]}.frcmod")
except Exception as e:
    print(f"Error in ligand file generation: {e}")

# Ensure ligand files exist
if not os.path.exists('ligand.mol2') or not os.path.exists('ligand.frcmod'):
    sys.exit("Ligand parameter generation failed. Check ligand PDB file.")

# Prepare receptor with pdb4amber
#print("Preparing receptor structure...")
#os.system(f'pdb4amber -i {os.path.basename(receptor_path)} -o recept.pdb --nohyd')

# Generate lib file for ligand with tleap

print("Generate amber parameters for complex...")
with open('tleap.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.gaff2
UNL = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
saveoff UNL ligand.lib
quit
""")
os.system(f'tleap -f tleap.in')

# Generate Amber files for the receptor
print("Generate amber parameters for receptor...")
with open('tleap_receptor.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.protein.ff19SB
source leaprc.water.opc
RECEPTOR = loadpdb receptor_h.pdb
saveamberparm RECEPTOR receptor.prmtop receptor.inpcrd
savepdb RECEPTOR receptor_solvated.pdb
quit
""")
os.system(f'tleap -f tleap_receptor.in')

# Ensure receptor files exist
if not os.path.exists('receptor.prmtop') or not os.path.exists('receptor.inpcrd'):
    sys.exit("TLEaP failed to generate receptor topology or coordinate files.")

# Generate Amber files for the ligand
print("Generate amber parameters for ligand...")
with open('tleap_ligand.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.gaff2
loadoff ligand.lib
loadamberparams ligand.frcmod
LIGAND = loadpdb ligand_h.pdb
saveamberparm LIGAND ligand.prmtop ligand.inpcrd
quit
""")
os.system(f'tleap -f tleap_ligand.in')

# Ensure ligand files exist
if not os.path.exists('ligand.prmtop') or not os.path.exists('ligand.inpcrd'):
    sys.exit("TLEaP failed to generate ligand topology or coordinate files.")

# Copy receptor and ligand files to MD directory
os.system(f'cp receptor.prmtop receptor.inpcrd ../md/')
os.system(f'cp ligand.prmtop ligand.inpcrd ../md/')

# Run tleap to find charge and box volume
print("Running tleap to find charge and box volume...")
with open('tleap.in', 'w') as tleap_file:
    tleap_file.write(f"""
    source leaprc.protein.ff19SB
    source leaprc.gaff2
    source leaprc.water.opc
    loadoff ligand.lib
    loadamberparams ligand.frcmod
    COMPLEX = loadpdb complex_h.pdb
    charge COMPLEX
    solvateOct COMPLEX OPCBOX 10.0
    quit
    """)

os.system('tleap -f tleap.in')

# Read leap.log file and obtain #waters and charge of the system
with open('leap.log', 'r') as fp:
    for line in fp:
        if 'total atoms in file:' in line:
            natoms = int(line.split()[-1])
        elif line.endswith('residues.\n'):
            nwat = int(line.split()[1])
        elif line.startswith('Total perturbed charge:'):
            charge = int(float(line.split()[-1]))


# Calculate number of Na+ and Cl- ions to add to the solvated system
No = (nwat*0.15)/56  # Expected # of ions
# The method is only valid if No/charge >= 1 - add a check
if charge != 0:
    if No/abs(charge) < 1:
        sys.exit('No/charge < 1')
# Adjust for the charge of the system
Na = math.ceil(No - charge/2)
Cl = math.ceil(No + charge/2)

# Generate amber parameters for complex with tleap
print("Generate amber parameters for complex...")
with open('tleap.in', 'w') as tleap_file:
    tleap_file.write(f"""
    source leaprc.protein.ff19SB
    source leaprc.gaff2
    source leaprc.water.opc
    loadoff ligand.lib
    loadamberparams ligand.frcmod
    COMPLEX = loadpdb complex_h.pdb
    solvateOct COMPLEX OPCBOX 10.0
    addionsrand COMPLEX Na+ {Na} Cl- {Cl}
    saveamberparm COMPLEX complex.prmtop complex.inpcrd
    savepdb COMPLEX complex_solvated.pdb
    quit
    """)

os.system(f'tleap -f tleap.in')

# Strip solvent and ions from the solvated complex
print("Stripping solvent and ions from the solvated complex...")
with open('strip.in', 'w') as strip_file:
    strip_file.write("""
        parm complex.prmtop
        trajin complex_solvated.pdb
        strip :WAT
        strip :Na+
        strip :Cl-
        trajout complex_stripped.pdb
         """)
os.system('cpptraj -i strip.in')

# Ensure the stripped PDB file exists
if not os.path.exists('complex_stripped.pdb'):
    sys.exit("Stripping solvent and ions failed. Check cpptraj input.")

# Generate Amber parameters for the stripped complex
print("Generating non-PBC Amber parameters for the stripped complex...")
with open('tleap_strip.in', 'w') as tleap_strip_file:
    tleap_strip_file.write(f"""
source leaprc.protein.ff19SB
source leaprc.gaff2
loadoff ligand.lib
loadamberparams ligand.frcmod
COMPLEX = loadpdb complex_stripped.pdb
saveamberparm COMPLEX complex_stripped.prmtop complex_stripped.inpcrd
quit
""")
os.system('tleap -f tleap_strip.in')

# Ensure stripped topology files exist
if not os.path.exists('complex_stripped.prmtop') or not os.path.exists('complex_stripped.inpcrd'):
    sys.exit("TLEaP failed to generate stripped topology or coordinate files.")

# Copy files to MD directory
os.system(f'cp complex.prmtop complex.inpcrd ../md/')

# Prepare input files for minimization, heating, equilibration, and production
os.chdir('../md')

# Input templates for MD
input_templates = {
    "01-min1.in": """Minimization 1 (solute heavy atoms restrained)
&cntrl
imin=1,
maxcyc=10000,
ncyc=5000,
ntpr=1000,
cut=10.0,
ntr=1,
restraint_wt=10.0,
restraintmask='!:WAT&!:Na+&!:Cl-&!:H'
/
""",
    "02-min2.in": """Minimization 2 (no restraints)
&cntrl
imin=1,
maxcyc=10000,
ncyc=5000,
ntpr=1000,
cut=10.0,
ntr=0
/
""",
    "03-heat.in": """Heating in NVT ensemble
&cntrl
nstlim=250000,
dt=0.002,
ntpr=500,
tempi=100.0,
temp0=310.0,
ntt=3,
gamma_ln=1.0,
cut=10.0,
ntb=1,
ntc=2,
ntf=2
/
""",
    "04-npt.in": """Equilibration in NPT ensemble
&cntrl
nstlim=250000,
dt=0.002,
ntpr=500,
tempi=310.0,
temp0=310.0,
ntt=3,
gamma_ln=1.0,
ntp=1,
ntb=2,
cut=10.0,
irest=1,
ntx=5,
ntc=2,
ntf=2
/
""",
    "05-prod.in": """Production (NPT ensemble)
&cntrl
nstlim=50000000,
dt=0.002,
ntpr=1000,
ntwx=1000,
tempi=310.0,
temp0=310.0,
ntt=3,
gamma_ln=1.0,
ntp=1,
ntb=2,
cut=10.0,
irest=1,
ntx=5,
ntc=2,
ntf=2
/
"""
}

# Write input files
for filename, content in input_templates.items():
    with open(filename, 'w') as file:
        file.write(content)

# GPU Run Script
with open('run_md.job', 'w') as job_file:
    job_file.write(f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca
#SBATCH --job-name=md_{complex_name}

{amber}

# Minimization with restraints
pmemd.cuda -O -i 01-min1.in -p complex.prmtop -c complex.inpcrd -o 01-min1.out -r 01-min1.rst -ref complex.inpcrd -inf 01-min1.info

# Full minimization
pmemd.cuda -O -i 02-min2.in -p complex.prmtop -c 01-min1.rst -o 02-min2.out -r 02-min2.rst -inf 02-min2.info

# Heating
pmemd.cuda -O -i 03-heat.in -p complex.prmtop -c 02-min2.rst -o 03-heat.out -r 03-heat.rst -x 03-heat.nc -inf 03-heat.info

# Equilibration
pmemd.cuda -O -i 04-npt.in -p complex.prmtop -c 03-heat.rst -o 04-npt.out -r 04-npt.rst -x 04-npt.nc -inf 04-npt.info

# Production
pmemd.cuda -O -i 05-prod.in -p complex.prmtop -c 04-npt.rst -o 05-prod.out -r 05-prod.rst -x 05-prod.nc -inf 05-prod.info
""")

print("Preparation complete. Move to the MD directory to start simulations.")