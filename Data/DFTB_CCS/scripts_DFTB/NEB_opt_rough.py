#!/Users/tjams20/opt/anaconda3/bin/python

from calc import *

import os 
from ase import Atoms
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
import pymatgen as pmg

def gen_spline(x,y):
    f = interpolate.interp1d(x, y, kind='quadratic') # ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’

    xnew = np.linspace(x[0], x[-1], num=100)
    ynew = f(xnew)
    return [xnew, ynew]

NEB_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN/PSTRESS_0/NEB_opt/"
NEB_folder_Cs = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/"

pressures = [0] # np.linspace(0,20,6)
NEB_paths = ['21_25', '21_26', '21_37', '21_38', '25_26', '25_36', '25_38']
NEB_paths_Cs = ['2_5', '2_6', '2_7']

for press in pressures:
    for NEB_path in NEB_paths_Cs:
        NEB_ini_id = NEB_path[0]
        NEB_fin_id = NEB_path[2]

        sup_40_SPC = pmg.Structure.from_file("/Users/tjams20/Documents/UBath/isambard_rsync/strain_study_to_isam_SCAN_320/PSTRESS_" + str(int(press)) + "/sup_40.POSCAR")

        # struc_ini = pmg.Structure.from_file(NEB_folder_Cs + "PSTRESS_" + str(int(press)) + "/vac_opt_Cs_large/vac_" + NEB_ini_id + "/CONTCAR")
        # struc_final = pmg.Structure.from_file(NEB_folder_Cs + "PSTRESS_" + str(int(press)) + "/vac_opt_Cs_large/vac_" + NEB_fin_id + "/CONTCAR")
        
        struc_ini = sup_40_SPC.copy()
        struc_final = sup_40_SPC.copy()

        struc_ini.remove_sites([int(NEB_ini_id)])
        struc_final.remove_sites([int(NEB_fin_id)])

        struc_ini.append(struc_ini.sites[int(NEB_fin_id)-1].species, struc_ini.sites[int(NEB_fin_id)-1].frac_coords)
        struc_final.append(struc_final.sites[int(NEB_ini_id)].species, struc_final.sites[int(NEB_ini_id)].frac_coords)

        struc_ini.remove_sites([int(NEB_fin_id)-1])
        struc_final.remove_sites([int(NEB_ini_id)])

        struc_ini.to(filename = 'initial_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.POSCAR')
        struc_final.to(filename = 'final_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.POSCAR')

        struc_ini = read('initial_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.POSCAR')
        struc_final = read('final_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.POSCAR')

        struc_ini.calc = calc_40
        struc_final.calc = calc_40
    
        # Initial state:
        if not os.path.isfile('initial_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.traj'):
            qn = QuasiNewton(struc_ini, trajectory='initial_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.traj')
            qn.run(fmax=0.01)
    
        # Final state:
        if not os.path.isfile('final_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.traj'):
            qn = QuasiNewton(struc_final, trajectory='final_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.traj')
            qn.run(fmax=0.01)
    
        initial = read('initial_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.traj')
        final = read('final_PSTRESS_' + str(int(press)) + '_' + NEB_path + '.traj')
    
        images = [initial]
        images += [initial.copy() for i in range(5)]
        images.append(final)
    
        neb = DyNEB(images, k=0.5, climb=False, fmax=0.50, dynamic_relaxation=True, scale_fmax=3.)
        neb.interpolate()
        for image in images[1:6]:
            calc_40 = LinearCombinationCalculator(calcs_40, weights)
            c = FixAtoms(indices=range(len(image)-1))
            c_at = FixAtoms(indices=[atom.index for atom in image if ((atom.symbol == 'Cs') or (atom.symbol == 'Pb'))])
            # image.set_constraint([c_at])
            red_weights=[1, 0.1]
            calc_red_weights=LinearCombinationCalculator(calcs_40, red_weights)
            image.calc = calc_red_weights
    
        opt = MDMin(neb, trajectory='neb_PSTRESS_' + str(int(press)) + '_' + NEB_path + '_rough.traj')
        for i in [4]:
            opt.attach(io.Trajectory('A2B-%d.traj' % i, 'w', images[i]))
        opt.run(fmax=0.50)
    
        os.system("ase nebplot neb_PSTRESS_" + str(int(press)) + '_' + NEB_path + "_rough.traj")
        os.system("mv nebplots.pdf nebplots_PSTRESS_" + str(int(press)) + "_" + NEB_path + "_rough.pdf")
