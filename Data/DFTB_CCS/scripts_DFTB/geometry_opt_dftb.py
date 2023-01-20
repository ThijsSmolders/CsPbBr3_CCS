#!/Users/tjams20/opt/anaconda3/bin/python

import os 
from ase import Atoms
from ase.calculators.dftb import Dftb
from ase.calculators.ccs import CCS
from ase.calculators.mixing import LinearCombinationCalculator 
from ase import io

from ase.optimize import QuasiNewton
from ase.optimize import FIRE
from ase.io import Trajectory, read, write
import ase.db as db 
from ase.units import *
import numpy as np
import json
from ase.calculators.mixing import LinearCombinationCalculator 

structure = read("/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/PSTRESS_0/pris_uc_opt_high_acc/OUTCAR", index=-1)

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
                         Hamiltonian_KPointsAndWeights = 'SupercellFolding {6 0 0 0 4 0 0 0 6 0.0 0.0 0.0}') 

with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)

eps_Cs = np.float(CCS_params['eps']['Cs'])
eps_Pb = np.float(CCS_params['eps']['Pb'])
eps_Br = np.float(CCS_params['eps']['Br'])

charge={'Cs': 0.0, 'Pb': 0.0, 'Br':0.0} 
eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
CCS_calc = CCS(charge=True, q=charge, eps=eps, rc=10)

calcs =[DFTB_calc, CCS_calc]
weights=[1, 1]
calc=LinearCombinationCalculator(calcs, weights) 

DFT_data="DFT.db"

DFT_DB=db.connect(DFT_data)

structure.calc = calc

nrg = structure.get_potential_energy()
print(nrg)

dyn = QuasiNewton(atoms=structure, trajectory='QN_DFTB.traj')
# dyn = BFGS(structure)
dyn.run(fmax=0.05)
