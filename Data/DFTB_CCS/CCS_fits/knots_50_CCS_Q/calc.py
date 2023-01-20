#!/Users/tjams20/opt/anaconda3/bin/python

import os
from ase.calculators.ccs import CCS
from ase.calculators.dftb import Dftb
from ase.calculators.mixing import LinearCombinationCalculator

import json

os.environ["DFTB_PREFIX"] = "/Users/tjams20/Documents/UBath/DFTB/DFTB_CCS/DFTB_PARAMS/Cs_w3_00_c9_00_Pb_w3_00_c9_00_Br_w3_00_c9_00/"

DFTB_calc_40 = Dftb(label='opt',
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
                     Hamiltonian_KPointsAndWeights = 'SupercellFolding { 2 0 0 0 2 0 0 0 2 0.0 0.0 0.0}')

DFTB_calc_320 = Dftb(label='opt',
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
                     Hamiltonian_KPointsAndWeights = 'SupercellFolding { 1 0 0 0 1 0 0 0 1 0.0 0.0 0.0}')

with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)

eps_Cs = float(CCS_params['eps']['Cs'])
eps_Pb = float(CCS_params['eps']['Pb'])
eps_Br = float(CCS_params['eps']['Br'])
charge={'Cs': 0.0, 'Pb': 0.0, 'Br':0.0}

eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
CCS_calc = CCS(charge=True, q=charge, eps=eps, rc=10)

calcs_40 =[DFTB_calc_40, CCS_calc]
weights=[1, 1]
calc_40=LinearCombinationCalculator(calcs_40, weights)

calcs_320 =[DFTB_calc_40, CCS_calc]
weights=[1, 1]
calc_320=LinearCombinationCalculator(calcs_320, weights)
