#!/Users/tjams20/opt/anaconda3/bin/python

import os
from ase import units
from ase.calculators.ccs import CCS
from ase.calculators.dftb import Dftb
from ase.calculators.mixing import LinearCombinationCalculator
from ase import Atoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import force_temperature, MaxwellBoltzmannDistribution
from ase.io import Trajectory, read, write
import ase.db as db
import json
import numpy as np
import sys
from pymatgen import Structure
from ase.constraints import FixAtoms

with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)
with open('CCS_input.json', 'r') as f:
   CCS_input = json.load(f)

eps_Cs = np.float(CCS_params['eps']['Cs'])
eps_Pb = np.float(CCS_params['eps']['Pb'])
eps_Br = np.float(CCS_params['eps']['Br'])

charge={'Cs': 0.0, 'Pb': 0.0, 'Br': 0.0} # REPEAT CP2K tjams20
eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
CCS_calc = CCS(charge=True,q=charge,eps=eps,rc=10)

os.environ["DFTB_PREFIX"] = "/Users/tjams20/Documents/UBath/DFTB/rhea_rsync/python_workflow/SKF/Cs_w3_00_c9_00_Pb_w3_00_c9_00_Br_w3_00_c9_00/"
DFTB_calc = Dftb(label='opt',
                         Hamiltonian_SCC='Yes',
                         Hamiltonian_ShellResolvedSCC = 'Yes',
                         Hamiltonian_SCCTolerance = "1.0E-006",
                         Hamiltonian_ReadInitialCharges = 'No',
                         Hamiltonian_MaxSCCIterations = "1500",
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_Cs='"s"',
                         Hamiltonian_MaxAngularMomentum_Pb='"p"',
                         Hamiltonian_MaxAngularMomentum_Br='"p"',
                         Hamiltonian_KPointsAndWeights = 'SupercellFolding {2 0 0 0 2 0 0 0 2 0.0 0.0 0.0}')

calcs =[DFTB_calc, CCS_calc]
weights=[1, 1]
calc=LinearCombinationCalculator(calcs, weights)

# structure = read("/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/PSTRESS_0/pris_uc_opt_high_acc/OUTCAR", index=-1)

# structure = read("/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN/PSTRESS_0/NEB_opt/21_25/02/OUTCAR", index=-1)
structure = read("/Users/tjams20/Documents/UBath/isambard_rsync/strain_study_to_isam_SCAN_320/pris_run_v2/PSTRESS_0/sup_40.POSCAR")

# structure.rattle(stdev=0.1)
Temp = int(sys.argv[1])
MaxwellBoltzmannDistribution(structure, temperature_K=500)
force_temperature(structure, Temp)
structure.calc = calc
c = FixAtoms(indices=[0])
structure.set_constraint(c)

nrg = structure.get_potential_energy()
print(nrg)

dyn = VelocityVerlet(structure, dt=5.0 * units.fs,
                     trajectory='md.traj', logfile='md.log')
dyn.run(500)  # take 1000 steps
