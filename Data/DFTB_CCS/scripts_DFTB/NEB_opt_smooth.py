#!/Users/tjams20/opt/anaconda3/bin/python

import os 
from ase import Atoms
from ase.calculators.dftb import Dftb
from ase.calculators.ccs import CCS
from ase.calculators.mixing import LinearCombinationCalculator 
from ase import io
from ase.neb import NEB
from ase.dyneb import DyNEB
from ase.optimize import MDMin
from ase.optimize import QuasiNewton
from ase.optimize import FIRE
from ase.io import Trajectory, read, write
import ase.db as db 
from ase.units import *
import numpy as np
import json
from ase.calculators.mixing import LinearCombinationCalculator 
from ase.constraints import FixAtoms
import matplotlib.pyplot as plt
from scipy import interpolate

def gen_spline(x,y):
    f = interpolate.interp1d(x, y, kind='quadratic') # ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’

    xnew = np.linspace(x[0], x[-1], num=100)
    ynew = f(xnew)
    return [xnew, ynew]

NEB_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN/PSTRESS_0/NEB_opt/"
NEB_folder_Cs = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/PSTRESS_0/NEB_opt_Cs/"

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
                         Hamiltonian_Filling = 'Fermi { Temperature [Kelvin] = 20.0}',
                         Hamiltonian_KPointsAndWeights = 'SupercellFolding {2 0 0 0 2 0 0 0 2 0.0 0.0 0.0}') 

with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)

eps_Cs = float(CCS_params['eps']['Cs'])
eps_Pb = float(CCS_params['eps']['Pb'])
eps_Br = float(CCS_params['eps']['Br'])

charge={'Cs': 0.0, 'Pb': 0.0, 'Br':0.0} 
eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
CCS_calc = CCS(charge=True, q=charge, eps=eps, rc=10)

calcs =[DFTB_calc, CCS_calc]
weights=[1, 1]
calc=LinearCombinationCalculator(calcs, weights) 

NEB_paths = ['21_25', '21_26', '21_37', '21_38', '25_26', '25_36', '25_38']
NEB_paths_Cs = ['2_5', '2_6', '2_7']

for NEB_path in NEB_paths_Cs:

    images = io.read('neb_' + NEB_path + '_rough.traj@-7:')

    neb = DyNEB(images, k=0.25, climb=False, fmax=0.02, dynamic_relaxation=True, scale_fmax=3.)
    for image in images[1:6]:
        calc = LinearCombinationCalculator(calcs, weights)
        image.calc = calc

    opt = MDMin(neb, trajectory='neb_' + NEB_path + '_smooth.traj')
    opt.run(fmax=0.02)

    os.system("ase nebplot neb_" + NEB_path + "_smooth.traj")
