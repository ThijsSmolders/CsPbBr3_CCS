#!/Users/tjams20/opt/anaconda3/bin/python

from calc import *

import os 
from ase import Atoms
from ase import io

from ase.optimize import QuasiNewton
from ase.optimize import FIRE
from ase.io import Trajectory, read, write
import ase.db as db 
from ase.units import *
import numpy as np

import matplotlib.pyplot as plt
from scipy import interpolate

def gen_spline(x,y):
    f = interpolate.interp1d(x, y, kind='quadratic') # ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’

    xnew = np.linspace(x[0], x[-1], num=100)
    ynew = f(xnew)
    return [xnew, ynew]

NEB_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN/PSTRESS_0/NEB_opt/"
NEB_folder_Cs = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/PSTRESS_0/NEB_opt_Cs/"

NEB_paths = ['21_25', '21_26', '21_37', '21_38', '25_26', '25_36', '25_38']
NEB_paths_Cs = ['2_5', '2_6', '2_7']

fig, axs = plt.subplots(3, 4, figsize=(16,8))

for NEB_id, NEB_path in enumerate(NEB_paths):
    structures = [read(NEB_folder + NEB_path + "/0" + str(i) + "/CONTCAR") for i in range(6)]
    
    cur_nebef = open(NEB_folder + NEB_path + "/nebef.dat")
    cur_path_DFT = []
    for line in cur_nebef.readlines():
        cur_path_DFT = np.append(cur_path_DFT, float(line.split(' ')[-2]))
    
    DFTB_CCS_NEB_nrgs = []
    for structure in structures:
        structure.calc = calc_40
    
        nrg = structure.get_potential_energy()
        DFTB_CCS_NEB_nrgs = np.append(DFTB_CCS_NEB_nrgs, nrg)
    
    DFTB_CCS_NEB_nrgs = DFTB_CCS_NEB_nrgs-DFTB_CCS_NEB_nrgs[0]

    [spl_DFT_x, spl_DFT_y] = gen_spline(range(6), cur_path_DFT)
    [spl_DFTB_x, spl_DFTB_y] = gen_spline(range(6), DFTB_CCS_NEB_nrgs)
    axs[int(NEB_id/4), int(NEB_id%4)].plot(spl_DFT_x, spl_DFT_y, 'k')
    axs[int(NEB_id/4), int(NEB_id%4)].plot(spl_DFTB_x, spl_DFTB_y, 'r')
    axs[int(NEB_id/4), int(NEB_id%4)].plot(range(6), cur_path_DFT, 'ks')
    axs[int(NEB_id/4), int(NEB_id%4)].plot(range(6), DFTB_CCS_NEB_nrgs, 'rs')
    axs[int(NEB_id/4), int(NEB_id%4)].set_title(NEB_path)

fig.savefig("NEB_paths_Br.pdf", dpi=300, bbox_inches = 'tight', pad_inches=0.35)

fig, axs = plt.subplots(1, 3, figsize=(9,3))

for NEB_id, NEB_path in enumerate(NEB_paths_Cs):
    structures = [read(NEB_folder_Cs + NEB_path + "/0" + str(i) + "/CONTCAR") for i in range(6)]
    
    cur_nebef = open(NEB_folder_Cs + NEB_path + "/nebef.dat")
    cur_path_DFT = []
    for line in cur_nebef.readlines():
        cur_path_DFT = np.append(cur_path_DFT, float(line.split(' ')[-2]))
    
    DFTB_CCS_NEB_nrgs = []
    for structure in structures:
        structure.calc = calc_320
    
        nrg = structure.get_potential_energy()
        DFTB_CCS_NEB_nrgs = np.append(DFTB_CCS_NEB_nrgs, nrg)
    
    DFTB_CCS_NEB_nrgs = DFTB_CCS_NEB_nrgs-DFTB_CCS_NEB_nrgs[0]

    [spl_DFT_x, spl_DFT_y] = gen_spline(range(6), cur_path_DFT)
    [spl_DFTB_x, spl_DFTB_y] = gen_spline(range(6), DFTB_CCS_NEB_nrgs)
    axs[NEB_id].plot(spl_DFT_x, spl_DFT_y, 'k')
    axs[NEB_id].plot(spl_DFTB_x, spl_DFTB_y, 'r')
    axs[NEB_id].plot(range(6), cur_path_DFT, 'ks')
    axs[NEB_id].plot(range(6), DFTB_CCS_NEB_nrgs, 'rs')
    axs[NEB_id].set_title(NEB_path)

axs[0].set_ylabel("Migration enth. (eV)")

fig.savefig("NEB_paths_Cs.pdf", dpi=300, bbox_inches = 'tight', pad_inches=0.35)
