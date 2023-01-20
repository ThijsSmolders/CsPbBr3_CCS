#!/Users/tjams20/opt/anaconda3/bin/python

from ase import units
from ase.calculators.ccs import CCS
from ase.calculators.dftb import Dftb
from ase import Atoms
from ase.optimize import QuasiNewton
from ase.io import Trajectory, read, write
from pymatgen.io.vasp.outputs import Vasprun
from ase.calculators.mixing import LinearCombinationCalculator
import numpy as np
import matplotlib.pyplot as plt
import json
import pickle
import os

tilt_folder =  "/Users/tjams20/Documents/UBath/DFTB/Training_Data/Tilting_Patterns/"
with open('/Users/tjams20/Documents/UBath/DFTB/df_phases.pkl', 'rb') as handle:
    df_phases = pickle.load(handle)
with open('CCS_params.json', 'r') as f:
    CCS_params = json.load(f)
os.environ["DFTB_PREFIX"] = "/Users/tjams20/Documents/UBath/DFTB/DFTB_CCS/DFTB_PARAMS/Cs_w3_00_c9_00_Pb_w3_00_c9_00_Br_w3_00_c9_00/"


def sg_dict(sg):
    sg_dict = {
        'C2/c': r'C2/c',
        'C2/m': r'C2/m',
        'Cmcm': r'Cmcm',
        'I4/mcm': r'I4/mcm',
        'I4/mmm': r'I4/mmm',
        'Im-3': r'Im$\bar{3}$',
        'Imma': r'Imma',
        'Immm': r'Immm',
        'P-1': r'P$\bar{1}$',
        'P2_1/m': r'P2$_1$/m',
        'P4/mbm': r'P4/mbm',
        'P4_2/nmc': r'P4$_2$/nmc',
        'Pm-3m': r'Pm$\bar{3}$m',
        'Pnma': r'Pnma',
        'R-3c': r'R$\bar{3}$c'
    }
    return sg_dict[sg]

sg_entries = ['Pm-3m', 'P4/mbm', 'I4/mmm', 'Im-3', 'Immm', 'I4/mcm', 'Imma', 'R-3c', 'C2/m', 'C2/c', 'P-1', 'Cmcm', 'Pnma', 'P2_1/m', 'P4_2/nmc']
sg_entries = np.sort(sg_entries)

eps_Cs = float(CCS_params['eps']['Cs'])
eps_Pb = float(CCS_params['eps']['Pb'])
eps_Br = float(CCS_params['eps']['Br'])
charge={'Cs': 0.0, 'Pb': 0.0, 'Br':0.0}

eps = {'Cs':eps_Cs, 'Pb':eps_Pb, 'Br':eps_Br}
CCS_calc = CCS(charge=True, q=charge, eps=eps, rc=10)

atoms_cub = read(tilt_folder + "Pm3m/OUTCAR", index=-1)
DFT_nrg_cub = atoms_cub.get_potential_energy()

DFT_nrgs = []
DFTB_CCS_nrgs = []

for sg_entry in sg_entries:
    atoms = read(tilt_folder + sg_entry.replace('/','_').replace('-','') + "/OUTCAR", index=-1)
    vrun = Vasprun(tilt_folder + sg_entry.replace('/','_').replace('-','') + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
    kpts = vrun.kpoints.kpts[0]
    DFT_nrg = atoms.get_potential_energy()
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
    atoms.calc = calc
    DFTB_CCS_nrg = atoms.get_potential_energy()
    fu = atoms.get_number_of_atoms()/5
    DFT_nrgs = np.append(DFT_nrgs, DFT_nrg/fu)
    
    qn = QuasiNewton(atoms, trajectory = sg_entry.replace('/','_') + '.traj')
    qn.run(fmax=0.01)
    
    opt_sg = read(sg_entry.replace('/','_') + '.traj')
    DFTB_CCS_nrg = opt_sg.get_potential_energy()

    DFTB_CCS_nrgs = np.append(DFTB_CCS_nrgs, DFTB_CCS_nrg/fu)
    print("DFT", sg_entry, DFT_nrg/fu)
    print(sg_entry, DFTB_CCS_nrg/fu)
    
plot_x = 0
plot_sgs = []

fig, ax = plt.subplots()

for index, DFT_nrg in enumerate(DFT_nrgs):
    ax.plot(plot_x, (DFT_nrg-DFT_nrg_cub)*1000, 's', markersize=15, color='k')
    ax.plot(plot_x, (DFTB_CCS_nrgs[index]-np.max(DFTB_CCS_nrgs))*1000, 's', markersize=12, color='r')
    plot_x += 1

plot_sgs = sg_entries
plot_sgs = [sg_dict(sg) for sg in plot_sgs]
ax.set_xlim(-1, plot_x) 
# ax.set_ylim(-150,10)
ax.set_xticks(np.arange(0, plot_x, 1))
ax.set_xticklabels(plot_sgs, rotation=270, fontsize=24)
ax.tick_params(axis='y', which='both', labelsize=24)
ax.hlines(0, -1, plot_x, linestyle='--', color='k', linewidth=2.2, zorder=0)
plt.ylabel("Energy (meV/f.u.)", fontsize=27)
fig.savefig("comp_tilt_DFTB_CCS_opt.pdf", dpi=300, bbox_inches = 'tight', pad_inches=0.35)
