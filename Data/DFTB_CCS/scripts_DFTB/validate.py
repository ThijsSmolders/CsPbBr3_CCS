#!/Users/tjams20/opt/anaconda3/bin/python

from ase.db import connect
from ase.io import Trajectory,read, write
import re
from ase.calculators.neighborlist  import *
import numpy as np
from pymatgen import Lattice, Structure
from pymatgen.analysis import ewald
from ase.visualize import view
from ase.calculators.ccs import CCS
import ase.db as db
import json

with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)
with open('CCS_input.json', 'r') as f:
    CCS_input = json.load(f)

q = np.float(CCS_params['Charge scaling factor'])
#q = np.sqrt(np.float(CCS_input["General"]["ewald_scaling"]))
eps_Cs = np.float(CCS_params['eps']['Cs'])
eps_Pb = np.float(CCS_params['eps']['Pb'])
eps_Br = np.float(CCS_params['eps']['Br'])

# charge={'Cs':1.0*q, 'Pb': 2.0*q, 'Br': -1.0*q}
# charge={'Cs': 0.835*q, 'Pb': 0.920*q, 'Br': -0.585*q} # Freddie values
# charge={'Cs': 1.0520*q, 'Pb': 0.9199*q, 'Br': -0.6573*q} # Balestra_2020
# charge = {'Cs':q, 'Pb':2*q, 'Br':-q}
charge={'Cs': 0.812*q, 'Pb': 1.330*q, 'Br':-0.714*q} # REPEAT CP2K tjams20
eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
calc = CCS(charge=True, q=charge, eps=eps, rc=10)
DFT_data="DFT.db"
DFT_DB=db.connect(DFT_data)

for row in DFT_DB.select():
    structure=row.toatoms()
    EDFT=structure.get_total_energy()
    structure.calc=calc
    print(EDFT,structure.get_potential_energy(), structure.get_potential_energy()-EDFT)
