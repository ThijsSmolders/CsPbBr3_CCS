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

import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate

plt.style.use('classic')
matplotlib.rcParams.update({'font.size': 30, 'xtick.major.pad': 14, 'ytick.major.pad': 14, 'axes.labelpad': 16}) # s, 'text.usetex': True}) #, 'usetex':True})


NEB_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN/PSTRESS_0/NEB_opt/"
NEB_folder_Cs = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/"

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
                         Hamiltonian_KPointsAndWeights = 'SupercellFolding {1 0 0 0 1 0 0 0 1 0.0 0.0 0.0}') 

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
pressures = np.linspace(0,20,6)

Ea_fors = np.zeros((len(pressures), len(NEB_paths_Cs)))
Ea_revs = np.zeros((len(pressures), len(NEB_paths_Cs)))

for press_id, press in enumerate(pressures):
    for NEB_id, NEB_path in enumerate(NEB_paths_Cs):
        struc_ini = read(NEB_folder_Cs + "PSTRESS_" + str(int(press)) + "/vac_opt_Cs_large/vac_" + NEB_path[0] + "/CONTCAR")
        struc_fin = read(NEB_folder_Cs + "PSTRESS_" + str(int(press)) + "/vac_opt_Cs_large/vac_" + NEB_path[-1] + "/CONTCAR")
        struc_dim = read(NEB_folder_Cs + "PSTRESS_" + str(int(press)) + "/dimer_opt_Cs_large/" + NEB_path + "/CENTCAR.POSCAR")
       
        structures = [struc_ini, struc_dim, struc_fin]
    
        DFTB_CCS_NEB_nrgs = []
        for structure in structures:
            structure.calc = calc
        
            nrg = structure.get_potential_energy()
            DFTB_CCS_NEB_nrgs = np.append(DFTB_CCS_NEB_nrgs, nrg)
       
        Ea_for = DFTB_CCS_NEB_nrgs[1]-DFTB_CCS_NEB_nrgs[0]
        Ea_rev = DFTB_CCS_NEB_nrgs[1]-DFTB_CCS_NEB_nrgs[2]

        Ea_fors[press_id, NEB_id] = Ea_for
        Ea_revs[press_id, NEB_id] = Ea_rev

        print(NEB_path + " at press " + str(int(press)) + " done.")

### Functions
def gen_spline(x,y):
    f = interpolate.interp1d(x, y, kind='quadratic') # ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’

    xnew = np.linspace(x[0], x[-1], num=100)
    ynew = f(xnew)
    return [xnew, ynew]

def color_dict(label):
    cur_dict = {
     '2_5': 'xkcd:green',
     '2_6': 'xkcd:blue',
     '2_7': 'xkcd:red'}
    return cur_dict[label]

def marker_dict(label):
    cur_dict = {
     '2_5': "o",
     '2_6': "s",
     '2_7': "D"}
    return cur_dict[label]

outline_thk = 4.05
line_thk = 4
mrk_size = 15
mrk_outline = 15.1
plt_frame_thk = 10

fig, ax = plt.subplots() # 1, figsize=(8,9))
plt.locator_params(axis='y', nbins=7)
for NEB_i, NEB_path in enumerate(NEB_paths_Cs):
    NEB_nrgs = Ea_fors[:,NEB_i]
    [spl_x, spl_y] = gen_spline(pressures/10, NEB_nrgs)
    ax.plot(spl_x, spl_y, linewidth=outline_thk, color='black')
    ax.plot(spl_x, spl_y, linewidth=line_thk, label=NEB_path, color=color_dict(NEB_path))
    ax.plot(pressures/10, NEB_nrgs, marker_dict(NEB_path), markersize=mrk_outline, color='k')
    ax.plot(pressures/10, NEB_nrgs, marker_dict(NEB_path), markersize=mrk_size, color=color_dict(NEB_path))

ax.set_xlim(-0.1, 2.1)
ax.set_ylim(0.6, 1.8)
ax.set_xlabel("Pressure (GPa)", fontsize=30)
ax.set_ylabel('Activation Ent. (eV)', fontsize=30)
# fig.show()
fig.savefig("mig_bar_Cs_Hydro_dimer_static.pdf", dpi=300, bbox_inches = 'tight', pad_inches=0.35)
