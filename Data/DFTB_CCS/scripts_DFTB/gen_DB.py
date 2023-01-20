#!/Users/tjams20/opt/anaconda3/bin/python

import os
import ase.db as db
from ase.io import Trajectory, read, write
from ase.calculators.dftb import Dftb
from pymatgen.io.vasp.outputs import Vasprun
from ase.calculators.ccs import CCS
from ase.calculators.mixing import LinearCombinationCalculator
import json

DFT_DB = db.connect("DFT.db")
DFTB_DB = db.connect("DFTB.db")

def write_to_DB(DFT_DB, structures_DFT, data): 
    try:
        for structure_DFT in structures_DFT:
            EDFT = structure_DFT.get_potential_energy()
            DFT_DB.write(structure_DFT, key=counter, data=data) 
            tmp = up_counter()
    
    except:
        structure_DFT = structures_DFT
        EDFT = structure_DFT.get_potential_energy()
        DFT_DB.write(structure_DFT, key=counter, data=data)
        tmp = up_counter()

os.environ["DFTB_PREFIX"] = "/Users/tjams20/Documents/UBath/DFTB/rhea_rsync/python_workflow/SKF/Cs_w3_00_c9_00_Pb_w3_00_c9_00_Br_w3_00_c9_00/"


with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)

eps_Cs = float(CCS_params['eps']['Cs'])
eps_Pb = float(CCS_params['eps']['Pb'])
eps_Br = float(CCS_params['eps']['Br'])

charge={'Cs': 0.0, 'Pb': 0.0, 'Br':0.0}
eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
CCS_calc = CCS(charge=True, q=charge, eps=eps, rc=10)


for row in DFT_DB.select():
    key_int=int(row.key)
    try:
        kpts = row.data.kpts
    except:
        kpts = [2,2,2]

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
                     Hamiltonian_KPointsAndWeights = 'SupercellFolding {{ {} 0 0 0 {} 0 0 0 {}     0.0 0.0 0.0}}'.format(kpts[0], kpts[1], kpts[2]))

    calcs =[DFTB_calc, CCS_calc]
    weights=[1, 1]
    calc=LinearCombinationCalculator(calcs, weights)

    struct = row.toatoms()
    struct.calc = calc
    DFTB_CCS_nrg = struct.get_potential_energy()
    print(DFTB_CCS_nrg)
    DFTB_DB.write(struct, key=key_int)
