#!/usr/bin/env python3
# md_mmgbsa_2lig.py
# Adapted from md_mmgbsa.py to accept two ligands.
# Author: (adapted by ChatGPT for Jordan Harrison)
# Keeps behavior and flow similar to original; runs antechamber once per ligand.

import os
import sys
import math
import subprocess
import shutil
import re

amber = 'module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber-pmemd/24.3'

# Ensure four arguments are passed (complex, receptor, ligand1, ligand2)
if len(sys.argv) != 5 or not all(arg.endswith('.pdb') for arg in sys.argv[1:5]):
    sys.exit('Please provide four PDB files: the complex (with two ligands), receptor, ligand1, and ligand2.')

# locate ligand_resname.txt (we expect two resnames, one per line: ligand1 then ligand2)
if os.path.exists('ligand_resname.txt'):
    fname = 'ligand_resname.txt'
elif os.path.exists('../ligand_resname.txt'):
    fname = '../ligand_resname.txt'
else:
    raise FileNotFoundError("ligand_resname.txt not found in current or parent directory. Expected two lines: ligand1_resname, ligand2_resname")

with open(fname) as f:
    resnames = [l.strip() for l in f.readlines() if l.strip()]

if len(resnames) < 2:
    sys.exit("ligand_resname.txt must contain two residue names (one per line): ligand1_resname then ligand2_resname")

lig_resname1, lig_resname2 = resnames[0], resnames[1]

# Input paths
complex_path = sys.argv[1]
receptor_path = sys.argv[2]
lig1_path = sys.argv[3]
lig2_path = sys.argv[4]

DIR = os.path.dirname(os.path.abspath(__file__))

# Create directories
complex_filename = os.path.basename(complex_path)
receptor_filename = os.path.basename(receptor_path)
lig1_filename = os.path.basename(lig1_path)
lig2_filename = os.path.basename(lig2_path)

os.makedirs('prep', exist_ok=True)
os.makedirs('md', exist_ok=True)

# Copy input files and the ligand_resname file into prep/
for fn in (complex_filename, receptor_filename, lig1_filename, lig2_filename):
    shutil.copy(fn, 'prep/')

shutil.copy(fname, 'prep/')

print("Current working directory:", os.getcwd())
print("Files in directory:", os.listdir())

# Move into prep directory
os.chdir('prep')

# (re)read resnames from the copy inside prep
with open('ligand_resname.txt') as f:
    resnames = [l.strip() for l in f.readlines() if l.strip()]

lig_resname1, lig_resname2 = resnames[0], resnames[1]

# Remove hydrogens from all inputs
os.system(f'pdb4amber --nohyd -i {complex_filename} -o complex_noh.pdb')
os.system(f'pdb4amber --nohyd -i {receptor_filename} -o receptor_noh.pdb')
os.system(f'pdb4amber --nohyd -i {lig1_filename} -o ligand1_noh.pdb')
os.system(f'pdb4amber --nohyd -i {lig2_filename} -o ligand2_noh.pdb')

# Add hydrogens with reduce
os.system('reduce complex_noh.pdb>complex_h.pdb')
os.system('reduce receptor_noh.pdb>receptor_h.pdb')
os.system('reduce ligand1_noh.pdb>ligand1_h.pdb')
os.system('reduce ligand2_noh.pdb>ligand2_h.pdb')

# Generate ligand parameters for ligand1
print("Generating ligand1 parameters...")
try:
    print(f"Running antechamber to generate ligand files for {lig1_filename}...")
    os.system(f"antechamber -i ligand1_h.pdb -fi pdb -o ligand1.mol2 -fo mol2 -at gaff")
    os.system(f"parmchk2 -i ligand1.mol2 -f mol2 -o ligand1.frcmod")
    print("Ligand1 files generated: ligand1.mol2, ligand1.frcmod")
except Exception as e:
    print(f"Error in ligand1 file generation: {e}")

# Generate ligand parameters for ligand2
print("Generating ligand2 parameters...")
try:
    print(f"Running antechamber to generate ligand files for {lig2_filename}...")
    os.system(f"antechamber -i ligand2_h.pdb -fi pdb -o ligand2.mol2 -fo mol2 -at gaff")
    os.system(f"parmchk2 -i ligand2.mol2 -f mol2 -o ligand2.frcmod")
    print("Ligand2 files generated: ligand2.mol2, ligand2.frcmod")
except Exception as e:
    print(f"Error in ligand2 file generation: {e}")

# Ensure ligand files exist
if not os.path.exists('ligand1.mol2') or not os.path.exists('ligand1.frcmod') or not os.path.exists('ligand2.mol2') or not os.path.exists('ligand2.frcmod'):
    sys.exit("Ligand parameter generation failed for one or both ligands. Check ligand PDB files and antechamber output.")

# Create lib files for both ligands with tleap
print("Creating ligand libraries...")
with open('tleap_liglibs.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.gaff2
{lig_resname1} = loadmol2 ligand1.mol2
{lig_resname2} = loadmol2 ligand2.mol2
loadamberparams ligand1.frcmod
loadamberparams ligand2.frcmod
saveoff {lig_resname1} ligand1.lib
saveoff {lig_resname2} ligand2.lib
quit
""")
os.system('tleap -f tleap_liglibs.in')

# Generate Amber files for receptor
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
os.system('tleap -f tleap_receptor.in')

if not os.path.exists('receptor.prmtop') or not os.path.exists('receptor.inpcrd'):
    sys.exit("TLEaP failed to generate receptor topology or coordinate files.")

# Generate Amber files for ligand1
print("Generate amber parameters for ligand1...")
with open('tleap_lig1.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.gaff2
loadoff ligand1.lib
loadamberparams ligand1.frcmod
LIG1 = loadpdb ligand1_h.pdb
saveamberparm LIG1 ligand1.prmtop ligand1.inpcrd
quit
""")
os.system('tleap -f tleap_lig1.in')

# Generate Amber files for ligand2
print("Generate amber parameters for ligand2...")
with open('tleap_lig2.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.gaff2
loadoff ligand2.lib
loadamberparams ligand2.frcmod
LIG2 = loadpdb ligand2_h.pdb
saveamberparm LIG2 ligand2.prmtop ligand2.inpcrd
quit
""")
os.system('tleap -f tleap_lig2.in')

if not (os.path.exists('lig1.prmtop') and os.path.exists('lig1.inpcrd') and os.path.exists('lig2.prmtop') and os.path.exists('lig2.inpcrd')):
    sys.exit("TLEaP failed to generate ligand topology or coordinate files for ligand1 and/or ligand2.")

# Copy receptor and ligand files to MD directory
os.system('cp receptor.prmtop receptor.inpcrd ../md/')
os.system('cp lig1.prmtop lig1.inpcrd ../md/')
os.system('cp lig2.prmtop lig2.inpcrd ../md/')

# Run tleap to find charge and box volume for the full complex (both ligands present)
print("Running tleap to find charge and box volume for complex (with both ligands)...")
with open('tleap_complex.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc
loadoff ligand1.lib
loadoff ligand2.lib
loadamberparams ligand1.frcmod
loadamberparams ligand2.frcmod
COMPLEX = loadpdb complex_h.pdb
charge COMPLEX
solvateOct COMPLEX OPCBOX 10.0
quit
""")
os.system('tleap -f tleap_complex.in')

# Parse leap.log for natoms, nwat, charge
natoms = None
nwat = None
charge = None
with open('leap.log', 'r') as fp:
    for line in fp:
        if 'total atoms in file:' in line:
            try:
                natoms = int(line.split()[-1])
            except:
                pass
        elif line.rstrip().endswith('residues.'):
            parts = line.split()
            # line like "   100 residues."
            try:
                nwat = int(parts[1])
            except:
                pass
        elif line.strip().startswith('Total perturbed charge:'):
            try:
                charge = int(float(line.split()[-1]))
            except:
                pass

if nwat is None or charge is None:
    # try alternative patterns if needed (keep original behavior tolerant)
    print("Warning: could not parse nwat or charge from leap.log; defaulting to charge=0 and nwat estimate")
    if nwat is None:
        nwat = 0
    if charge is None:
        charge = 0

# Calculate number of Na+ and Cl- ions to add
No = (nwat * 0.15) / 56.0
if charge != 0 and No / abs(charge) < 1:
    sys.exit('No/charge < 1')

Na = math.ceil(No - charge / 2.0)
Cl = math.ceil(No + charge / 2.0)

# Generate Amber parameters for complex (with both ligands) and solvate/add ions
print("Generate amber parameters for complex (with both ligands) and solvate/add ions...")
with open('tleap_complex_build.in', 'w') as tleap_file:
    tleap_file.write(f"""
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.opc
loadoff ligand1.lib
loadoff ligand2.lib
loadamberparams ligand1.frcmod
loadamberparams ligand2.frcmod
COMPLEX = loadpdb complex_h.pdb
solvateOct COMPLEX OPCBOX 10.0
addionsrand COMPLEX Na+ {Na} Cl- {Cl}
saveamberparm COMPLEX complex.prmtop complex.inpcrd
savepdb COMPLEX complex_solvated.pdb
quit
""")
os.system('tleap -f tleap_complex_build.in')

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

if not os.path.exists('complex_stripped.pdb'):
    sys.exit("Stripping solvent and ions failed. Check cpptraj input.")

# Generate non-PBC Amber parameters for the stripped complex (with both ligands)
print("Generating non-PBC Amber parameters for the stripped complex...")
with open('tleap_strip.in', 'w') as tleap_strip_file:
    tleap_strip_file.write(f"""
source leaprc.protein.ff19SB
source leaprc.gaff2
loadoff ligand1.lib
loadoff ligand2.lib
loadamberparams ligand1.frcmod
loadamberparams ligand2.frcmod
COMPLEX = loadpdb complex_stripped.pdb
saveamberparm COMPLEX complex_stripped.prmtop complex_stripped.inpcrd
quit
""")
os.system('tleap -f tleap_strip.in')

if not os.path.exists('complex_stripped.prmtop') or not os.path.exists('complex_stripped.inpcrd'):
    sys.exit("TLEaP failed to generate stripped topology or coordinate files.")

# Copy complex topology/coords to md directory
os.system('cp complex.prmtop complex.inpcrd ../md/')

# Move into md directory and create MD input files similar to original script
os.chdir('../md')

input_templates = {
    "01-min1.in": """Minimization 1 (solute heavy atoms restrained)
&cntrl
imin=1,
ntmin=2,
maxcyc=500000,
ntpr=1000,
cut=10.0,
ntr=1,
drms=0.05,
restraint_wt=10.0,
restraintmask='!:WAT&!:Na+&!:Cl-&!:H'
/
""",
    "02-min2.in": """Minimization 2 (Backbone heavy atoms restrained)
&cntrl
imin=1,
ntmin=2,
drms=0.05,
maxcyc=500000,
ntpr=500,
cut=10.0,
ntr=1,
restraint_wt=2.0,
restraintmask='@CA,C,N,O'
/
""",
    "03-min3.in": """Minimization 3 (no restraints)
&cntrl
imin=1,
ntmin=2,
maxcyc=500000,
ntpr=1000,
drms=0.05,
cut=10.0,
ntr=0
/
""",
    "04-heat.in": """Heating in NVT ensemble
&cntrl
nstlim=500000,
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
    "05-npt.in": """Equilibration in NPT ensemble
&cntrl
nstlim=500000,
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
    "06-prod.in": """Production (NPT ensemble)
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

for filename, content in input_templates.items():
    with open(filename, 'w') as file:
        file.write(content)

md_dir = os.getcwd()
protein_name = os.path.basename(os.path.dirname(os.path.dirname(md_dir)))
complex_dir = os.path.basename(os.path.dirname(md_dir))
match = re.search(r'ternary_complex(\d+)', complex_dir)
complex_number = match.group(1) if match else 'X'
job_name = f"{protein_name}_{complex_number}"

# GPU Run Script (very similar to original)
with open('run_md.job', 'w') as job_file:
    job_file.write(f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --gres=gpu:h100:1
#SBATCH --time=0-72:00:00
#SBATCH --account=def-aminpour
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaharri1@ualberta.ca
#SBATCH --job-name={job_name}

module --force purge
{amber}
module load amber-pmemd/24.3
module list

for i in {{1..20}}; do
    if which pmemd.cuda > /dev/null 2>&1; then
        echo "Amber CUDA path: $(which pmemd.cuda)"
        break
    else
        echo "pmemd.cuda not found in PATH, retry $i/20..."
        sleep 5
        module --force purge
        {amber}
        module load amber-pmemd/24.3
    fi
    if [ $i -eq 20 ]; then
        echo "pmemd.cuda could not be found after 20 attempts. Exiting."
        exit 1
    fi
done

echo "Amber CUDA path: $(which pmemd.cuda)"
module list

# Minimization 1
pmemd.cuda -O -i 01-min1.in -p complex.prmtop -c complex.inpcrd -o 01-min1.out -r 01-min1.rst -ref complex.inpcrd -inf 01-min1.info

# Minimization 2
pmemd.cuda -O -i 02-min2.in -p complex.prmtop -c 01-min1.rst -o 02-min2.out -r 02-min2.rst -ref complex.inpcrd -inf 02-min2.info

# Minimization 3
pmemd.cuda -O -i 03-min3.in -p complex.prmtop -c 02-min2.rst -o 03-min3.out -r 03-min3.rst -inf 03-min3.info

# Heating
pmemd.cuda -O -i 04-heat.in -p complex.prmtop -c 03-min3.rst -o 04-heat.out -r 04-heat.rst -x 04-heat.nc -inf 04-heat.info

# Equilibration
pmemd.cuda -O -i 05-npt.in -p complex.prmtop -c 04-heat.rst -o 05-npt.out -r 05-npt.rst -x 05-npt.nc -inf 05-npt.info

# Production
pmemd.cuda -O -i 06-prod.in -p complex.prmtop -c 05-npt.rst -o 06-prod.out -r 06-prod.rst -x 06-prod.nc -inf 06-prod.info

# === run_mmgbsa.py ===
echo "Submitting run_mmgbsa.py to slurm after md_mmgbsa.py completes..."
mmgbsa_jobid=$(sbatch --parsable --dependency=afterok:$md_jobid --job-name=mmgbsa --output=mmgbsa.out --error=mmgbsa.err --wrap="python run_mmgbsa.py")
echo "Submitted run_mmgbsa.py as job $mmgbsa_jobid (after md.py)"

# === analysis.py ===
echo "Submitting analysis.py to SLURM after md.py completes..."
analysis_jobid=$(sbatch --parsable --dependency=afterok:$md_jobid --job-name=analyzpy --output=analyzpy.out --error=analyzpy.err --wrap="python analysis.py")
echo "Submitted analysis.py as job $analysis_jobid (after md.py)"
""")

print("Preparation complete. Submitting job now...")
subprocess.run(['sbatch', 'run_md.job'])
