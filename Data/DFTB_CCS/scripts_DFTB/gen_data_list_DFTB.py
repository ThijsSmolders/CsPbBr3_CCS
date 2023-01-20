#!/Users/tjams20/opt/anaconda3/bin/python

import numpy as np
import os
import ase.db as db
from ase.io import Trajectory, read, write
from ase.calculators.dftb import Dftb
from pymatgen.io.vasp.outputs import Vasprun

Hydro = True
Biax = True
UniAx = True
Phase_Diag = False
Tilts = True
MD = False
Mode_Following = True
NEB_40 = False
NEB_320 = True

load_indices = -1

master_DB = db.connect("master.db")
DFT_DB = db.connect("DFT.db")

counter = 1

def up_counter(add=1):
    global counter
    counter += add
    return counter

def write_to_DB(master_DB, DFT_DB, structures_DFT, data): 
    try:
        for structure_DFT in structures_DFT:
            EDFT = structure_DFT.get_potential_energy()
            master_DB.write(structure_DFT)
            DFT_DB.write(structure_DFT, key=counter, data=data) 
            tmp = up_counter()
    
    except:
        structure_DFT = structures_DFT
        EDFT = structure_DFT.get_potential_energy()
        master_DB.write(structure_DFT)
        DFT_DB.write(structure_DFT, key=counter, data=data)
        tmp = up_counter()

f = open("data.list", "w")

Training_Data_folder = "/Users/tjams20/Documents/UBath/DFTB/Training_Data/"

### Hydrostatic pressure optimised unit cells

if Hydro:
    Hydro_folder_uc = Training_Data_folder + "Hydro/"

    pressures = np.linspace(0,20,6)

    for press in pressures:
        tmp_folder = Hydro_folder_uc + "PSTRESS_" + str(int(press)) + "/uc_opt/"
        f.write(tmp_folder + "\n")
        structures_DFT = read(tmp_folder + "/OUTCAR", index=load_indices) 
        vrun = Vasprun(tmp_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
        kpts = vrun.kpoints
        data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
        write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Biaxial strains optimised unit cells
if Biax:
    Biax_folder_uc = Training_Data_folder + "Biax/"
    strain_mags = np.linspace(1,5,5) #np.concatenate((np.linspace(-5,-1,5,np.linspace(1,5,5)))
    strain_planes = ["xy_opt_z", "xz_opt_y", "yz_opt_x"]

    strain_strs = []
    for strain_plane in strain_planes:
        for strain_mag in strain_mags:
            strain_mag_str = "{:.2f}".format(strain_mag).replace(".","_").replace("-","")
            if strain_mag == 0:
                strain_str = "pris_uc_opt_high_acc"
            elif strain_mag < 0:
                strain_str = "strain_" + strain_plane[:2] + "_m" + strain_mag_str + strain_plane[2:]
            elif strain_mag > 0:
                strain_str = "strain_" + strain_plane[:2] + "_p" + strain_mag_str + strain_plane[2:]
            if strain_str not in strain_strs:
                strain_strs.append(strain_str)
                tmp_folder = Biax_folder_uc + strain_str + "/"
                f.write(tmp_folder + "\n")
                structures_DFT = read(tmp_folder + "/OUTCAR", index=load_indices)
                vrun = Vasprun(tmp_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
                kpts = vrun.kpoints
                data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
                write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Biaxial strains optimised unit cells
if UniAx:
    UniAx_folder_uc = Training_Data_folder + "UniAx/"
    strain_mags = np.linspace(1,5,5)
    strain_planes = ["x_opt_yz", "y_opt_xz", "z_opt_xy"]

    strain_strs = []
    for strain_plane in strain_planes:
        for strain_mag in strain_mags:
            strain_mag_str = "{:.2f}".format(strain_mag).replace(".","_").replace("-","")
            if strain_mag == 0:
                strain_str = "pris_uc_opt_high_acc"
            elif strain_mag < 0:
                strain_str = "strain_" + strain_plane[:1] + "_m" + strain_mag_str + strain_plane[1:]
            elif strain_mag > 0:
                strain_str = "strain_" + strain_plane[:1] + "_p" + strain_mag_str + strain_plane[1:]
            if strain_str not in strain_strs:
                strain_strs.append(strain_str)
                tmp_folder = UniAx_folder_uc + strain_str + "/"
                f.write(tmp_folder + "\n")
                structures_DFT = read(tmp_folder + "/OUTCAR", index=load_indices)
                vrun = Vasprun(tmp_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
                kpts = vrun.kpoints
                data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
                write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Different stoichiometries
if Phase_Diag:
    Phase_diag_folder = Training_Data_folder + "Phase_Diagram/"
    Phases = ["Cs4PbBr6", "CsBr", "CsPb2Br5", "PbBr2"]
    for phase in Phases:
        tmp_folder = Phase_diag_folder + phase + "/"
        f.write(tmp_folder + "\n")
        structures_DFT = read(tmp_folder + "/OUTCAR", index=load_indices)
        vrun = Vasprun(tmp_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
        kpts = vrun.kpoints
        data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
        write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Different Tilting Patterns
if Tilts:
    Tilts_folder = Training_Data_folder + "Tilting_Patterns/"
    Tilts = ["C2_c", "C2_m", "Cmcm", "I4_mcm", "I4_mmm", "Im3", "Imma", "Immm", "P1", "P2_1_m", "P4_2_nmc", "P4_mbm", "Pm3m", "R3c"]
    for tilt in Tilts:
        tmp_folder = Tilts_folder + tilt + "/"
        f.write(tmp_folder + "\n")
        structures_DFT = read(tmp_folder + "/OUTCAR", index=load_indices)
        vrun = Vasprun(tmp_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
        kpts = vrun.kpoints
        data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
        write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Structures from MD runs
if MD:
    MD_folder = Training_Data_folder + "MD_runs/"
    runs = [1,2,3,4,5,6,7,8,9,10]
    for run in runs:
        run_folder = MD_folder + "run_" + str(int(run)) + "/"
        MD_dirs = [name for name in os.listdir(run_folder)]
        for MD_dir in MD_dirs:
            f.write(run_folder + MD_dir + "\n")
            structures_DFT = read(run_folder + MD_dir + "/OUTCAR", index=load_indices)
            vrun = Vasprun(run_folder + MD_dir + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Structures from mode-following
if Mode_Following:
    Mode_folder = Training_Data_folder + "Mode_Following/"
    phases = ["Pm3m", "I4_mcm", "P4_mbm", "Pnma"]
    for phase in phases:
        phase_folder = Mode_folder + phase + "/"
        step_folders = [name for name in os.listdir(phase_folder)]
        for step_folder in step_folders:
            f.write(phase_folder + step_folder + "\n")
            structures_DFT = read(phase_folder + step_folder + "/OUTCAR", index=load_indices)
            vrun = Vasprun(phase_folder + step_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)

### Structures from 40-atom study
pressures = np.linspace(0,20,6)
vacs_Br = [17, 21, 25, 26, 33, 36, 37, 38]
vacs_Cs = [2, 5, 6, 7]

NEB_paths_Br = ["21_25","21_26","21_37","21_38","25_26","25_36","25_38","21_17","21_27","25_33"]
NEB_paths_Cs = ["2_5","2_6","2_7"]

Hydro_strain_40_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN/"
Hydro_strain_320_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/"

if NEB_40 == True:
    for press in pressures:
        press_folder = Hydro_strain_40_folder + "PSTRESS_" + str(int(press)) + "/"
        for vac in vacs_Br:
            vac_folder = press_folder + "vac_opt/vac_" + str(vac)
            f.write(vac_folder + "\n")
            structures_DFT = read(vac_folder + "/OUTCAR", index=load_indices)
            vrun = Vasprun(vac_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)

        for NEB_path in NEB_paths_Br:
            NEB_folder = press_folder + "NEB_opt/" + NEB_path + "/"
            for img in range(4):
                f.write(NEB_folder + "0" + str(img+1) + "\n")
                structures_DFT = read(NEB_folder + "0" + str(img+1) + "/OUTCAR", index=load_indices)
                print(NEB_folder + "/vasprun.xml")
                vrun = Vasprun(NEB_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
                kpts = vrun.kpoints
                data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
                write_to_DB(master_DB, DFT_DB, structures_DFT, data)

if NEB_320 == True:
    for press in pressures:
        press_folder = Hydro_strain_320_folder = "/Users/tjams20/Documents/UBath/isambard_rsync/CsPbBr3_strain/strain_study_SCAN_320/pris_run_v2/PSTRESS_" + str(int(press)) + "/"
#         for vac in vacs_Cs:
#             vac_folder = press_folder + "vac_opt_Cs_large/vac_" + str(vac)
#             f.write(vac_folder + "\n")
#             structures_DFT = read(vac_folder + "/OUTCAR", index=load_indices)
#             write_to_DB(master_DB, DFT_DB, structures_DFT, counter, "PSTRESS_" + str(int(press)), "Cs_320", vac)
#         
#         for NEB_path in NEB_paths_Cs:
#             dim_folder = press_folder + "dimer_opt_Cs_large/" + NEB_path + "/"
#             f.write(dim_folder + NEB_path + "\n")
#                 structures_DFT = read(dim_folder + NEB_path + "/OUTCAR", index=load_indices)
#                 write_to_DB(master_DB, DFT_DB, structures_DFT, counter, "PSTRESS_" + str(int(press)), "Cs_320", "Cs_" + NEB_path)

        for vac in vacs_Br:
            vac_folder = press_folder + "vac_opt/vac_" + str(vac)
            f.write(vac_folder + "\n")
            structures_DFT = read(vac_folder + "/OUTCAR", index=load_indices)
            vrun = Vasprun(vac_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)
        
        for NEB_path in NEB_paths_Br:
            NEB_folder = press_folder + "NEB_opt_4/" + NEB_path + "/"
            for img in range(4):
                f.write(NEB_folder + "0" + str(img+1) + "\n")
            dim_folder = press_folder + "dimer_opt/" + NEB_path + "/dim/"
            f.write(dim_folder + "\n")
            structures_DFT = read(dim_folder + "/OUTCAR", index=load_indices)
            vrun = Vasprun(dim_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)

        for vac in vacs_Cs:
            vac_folder = press_folder + "vac_opt_Cs_large/vac_" + str(vac) + "/"
            f.write(vac_folder + "\n")
            structures_DFT = read(vac_folder + "/OUTCAR", index=load_indices)
            vrun = Vasprun(vac_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)

        for NEB_path in NEB_paths_Cs:
            dim_folder = press_folder + "dimer_opt_Cs_large/" + NEB_path
            f.write(dim_folder + "\n")
            structures_DFT = read(dim_folder + "/OUTCAR", index=load_indices)
            vrun = Vasprun(dim_folder + "/vasprun.xml", ionic_step_skip=None, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            kpts = vrun.kpoints
            data = {'kpts': kpts.kpts[0]} #{"KPOINTS": KPOINTS}
            write_to_DB(master_DB, DFT_DB, structures_DFT, data)

f.close()

### Build DFTB_db
DFT_DB = db.connect("DFT.db")
DFTB_DB = db.connect("DFTB.db")

os.environ["DFTB_PREFIX"] = "/Users/tjams20/Documents/UBath/DFTB/rhea_rsync/python_workflow/SKF/Cs_w3_00_c9_00_Pb_w3_00_c9_00_Br_w3_00_c9_00/"

for row in DFT_DB.select():
    key_int=int(row.key)
    # path_key=str(row.Path)
    # spec_key=str(row.Spec)
    # strain_key=str(row.Strain)
    # print(key_int, path_key, spec_key, strain_key)
    kpts = row.data.kpts

    print(key_int)
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
                         Hamiltonian_KPointsAndWeights = 'SupercellFolding {{ {} 0 0 0 {} 0 0 0 {}     0.0 0.0 0.0}}'.format(kpts[0], kpts[1], kpts[2]))
    
    struct = row.toatoms()
    struct.calc = DFTB_calc
    DFTB_nrg = struct.get_potential_energy()
    DFTB_DB.write(struct, key=key_int)  #, Spec=spec_key, Path=path_key, key=key_int)
