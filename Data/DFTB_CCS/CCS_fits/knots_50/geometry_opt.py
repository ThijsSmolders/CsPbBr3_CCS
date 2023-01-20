#!/Users/tjams20/opt/anaconda3/bin/python

from ase.calculators.ccs import CCS
from ase import Atoms
from ase.optimize import QuasiNewton
from ase.io import Trajectory, read, write
import json
from pymatgen import Structure
import numpy as np
import matplotlib.pyplot as plt
import ase.db as db

with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)

q = float(CCS_params["Charge scaling factor"])
eps_Cs = float(CCS_params["eps"]["Cs"])
eps_Pb = float(CCS_params["eps"]["Pb"])
eps_Br = float(CCS_params["eps"]["Br"])
# charge = {'Cs':q, 'Pb':2*q, 'Br':-q}
charge_dict={'Cs': 0.812*q, 'Pb': 1.330*q, 'Br':-0.714*q} # REPEAT CP2K tjams20
eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
calc = CCS(charge=True, q=charge_dict, eps=eps)

structure = read("/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/PSTRESS_0/pris_uc_opt_high_acc/OUTCAR", index=-1)

nrgs = []
for d in np.linspace(3.0,5.5,26):
    br2 = Atoms('Br2', positions=[(0, 0, 0), (0, 0, d)])
    br2.set_cell(10 * np.identity(3))
    print(br2)
    br2.calc = calc
    nrg = br2.get_potential_energy()
    nrgs = np.append(nrgs, nrg)

plt.plot(np.linspace(3.0,5.5,26), nrgs, marker='o')
plt.show()
