# coding: utf-8
# Copyright (c) Thijs JAM Smolders.

try:
    # New Py>=3.5 import
    from math import gcd
except ImportError:
    # Deprecated import from Py3.5 onwards.
    from fractions import gcd

import six

import numpy as np
import os
import pymatgen as pmg
import time

from monty.io import zopen

"""
This module provides classes used work with perovskite structures.
"""

__author__ = "Thijs JAM Smolders"
__copyright__ = "Copyright to author"
__version__ = "1.0"
__maintainer__ = "Thijs JAM Smolders"
__email__ = "smolders.thijs@gmail.com"
__status__ = "Development"
__date__ = "Dec 11, 2019"

def init_xml_non_disp(f, init_struc_inp):
    f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    f.write("<modeling>\n")
    f.write(" <generator>\n")
    f.write('  <i name="program" type="string">vasp </i>\n')
    f.write(" </generator>\n")
    f.write(" <atominfo>\n")
    f.write("  <atoms>   " + str(len(init_struc_inp.species)) + " </atoms>\n")
    f.write("  <types>   " + str(len(init_struc_inp.types_of_specie)) + " </types>\n")
    f.write('  <array name="atoms" >\n')
    f.write('   <dimension dim="1">ion</dimension>\n')
    f.write('   <field type="string">element</field>\n')
    f.write('   <field type="int">atomtype</field>\n')
    f.write('   <set>\n')
    spec_cnt = 0
    cur_spec = "INIT"
    for spec in init_struc_inp.species:
        if spec.name != cur_spec:
            spec_cnt += 1
            cur_spec = spec.name
        f.write('    <rc><c>' + spec.name + "</c><c>   " + str(spec_cnt) + '</c></rc>\n')
    f.write('   </set>\n')
    f.write('  </array>\n')
    f.write('  </atominfo>\n')

def write_xml_non_disp(f, struc):
    f.write(' <calculation>\n')
    f.write('  <structure>\n')
    f.write('  <crystal>\n')
    f.write('   <varray name="basis" >\n')
    latt = struc.lattice.matrix
    f.write('    <v> ' + str(latt[0,0]) + '   ' + str(latt[0,1]) + '   ' + str(latt[0,2]) + ' </v>\n')
    f.write('    <v> ' + str(latt[1,0]) + '   ' + str(latt[1,1]) + '   ' + str(latt[1,2]) + ' </v>\n')
    f.write('    <v> ' + str(latt[2,0]) + '   ' + str(latt[2,1]) + '   ' + str(latt[2,2]) + ' </v>\n')
    f.write('   </varray>\n')
    f.write('  </crystal>\n')
    f.write('  <varray name="positions" >\n')
    coords = struc.frac_coords
    for i in range(len(coords[:,0])):
        f.write('    <v> ' + str(coords[i,0]) + '   ' + str(coords[i,1]) + '   ' + str(coords[i,2]) + ' </v>\n')
    f.write('   </varray>\n')
    f.write('  </structure>\n')
    f.write(' </calculation>\n')

def init_xml_disp(f, init_struc_inp):
    f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    f.write("<modeling>\n")
    f.write(" <generator>\n")
    f.write('  <i name="program" type="string">vasp </i>\n')
    f.write(" </generator>\n")
    f.write(" <atominfo>\n")
    f.write("  <atoms>   " + str(len(init_struc_inp.species)) + " </atoms>\n")
    f.write("  <types>   " + str(len(init_struc_inp.types_of_specie)) + " </types>\n")
    f.write('  <array name="atoms" >\n')
    f.write('   <dimension dim="1">ion</dimension>\n')
    f.write('   <field type="string">element</field>\n')
    f.write('   <field type="int">atomtype</field>\n')
    f.write('   <set>\n')
    spec_cnt = 0
    cur_spec = "INIT"
    for spec in init_struc_inp.species:
        if spec.name != cur_spec:
            spec_cnt += 1
            cur_spec = spec.name
        f.write('    <rc><c>' + spec.name + "</c><c>   " + str(spec_cnt) + '</c></rc>\n')
    f.write('   </set>\n')
    f.write('  </array>\n')
    f.write('  </atominfo>\n')

def write_xml_disp(f, struc):
    f.write(' <calculation>\n')
    f.write('  <structure>\n')
    f.write('  <crystal>\n')
    f.write('   <varray name="basis" >\n')
    latt = struc.lattice.matrix
    f.write('    <v> ' + str(latt[0,0]) + '   ' + str(latt[0,1]) + '   ' + str(latt[0,2]) + ' </v>\n')
    f.write('    <v> ' + str(latt[1,0]) + '   ' + str(latt[1,1]) + '   ' + str(latt[1,2]) + ' </v>\n')
    f.write('    <v> ' + str(latt[2,0]) + '   ' + str(latt[2,1]) + '   ' + str(latt[2,2]) + ' </v>\n')
    f.write('   </varray>\n')
    f.write('  </crystal>\n')
    f.write('  <varray name="positions" >\n')
    coords = struc.frac_coords
    for i in range(len(coords[:,0])):
        f.write('    <v> ' + str(np.mod(coords[i,0]+0.0625, 1.0)) + '   ' + str(np.mod(coords[i,1]+0.0625, 1.0)) + '   ' + str(np.mod(coords[i,2]+0.0625, 1.0)) + ' </v>\n')
        # f.write('    <v> ' + str(coords[i,0]) + '   ' + str(coords[i,1]) + '   ' + str(coords[i,2]) + ' </v>\n')
    f.write('   </varray>\n')
    f.write('  </structure>\n')
    f.write(' </calculation>\n')

def init_xml(f, struc, trans):
    if trans==False:
        init_xml_non_disp(f, struc)
    elif trans==True:
        init_xml_disp(f, struc)

def write_xml(f, struc, trans):
    if trans==False:
        write_xml_non_disp(f, struc)
    elif trans==True:
        write_xml_disp(f, struc)

def write_NEBrun(nr_imgs=0, interpolate_rate=6, FILE_TYPE='CONTCAR', trans=False):
    NEB_folder = os.getcwd() + "/"
    NEBrun_path = NEB_folder + "NEBrun.xml"
    with open(NEB_folder + 'INCAR', 'r') as INC:
        for line in INC.readlines():
            if 'IMAGES' in line:
                dat = line.split(' ')
                nr_imgs = int(dat[2])
    init_struc = pmg.Structure.from_file(NEB_folder + "00/" + FILE_TYPE)
    with open(NEBrun_path, "w") as fle:
        init_xml(fle, init_struc, trans)
        for img in range(nr_imgs + 2):
            #if img == 1:
            #    continue
            if img < 10:
                img_str = '0' + str(img)
                print(img_str)
            else:
                img_str = str(img)
            if interpolate_rate > 0:
                tmp_interpolate_rate = interpolate_rate + 2
                if img != nr_imgs + 1:
                    next_img = img + 1
                    if next_img < 10:
                        next_img_str = '0' + str(next_img)
                    else:
                        next_img_str = str(next_img)
                    img_struc = pmg.Structure.from_file(NEB_folder + img_str + '/' + FILE_TYPE)
                    next_img_struc = pmg.Structure.from_file(NEB_folder + next_img_str + '/' + FILE_TYPE)
                    diff_frac_unnorm = next_img_struc.frac_coords - img_struc.frac_coords
                    diff_frac_unnorm = np.where(diff_frac_unnorm>0.5, diff_frac_unnorm-1, diff_frac_unnorm)
                    diff_frac = np.where(diff_frac_unnorm<-0.5, diff_frac_unnorm+1, diff_frac_unnorm) / tmp_interpolate_rate
                    diff_latt = (next_img_struc.lattice.matrix - img_struc.lattice.matrix) / tmp_interpolate_rate
                    for cur_int in range(tmp_interpolate_rate-1): # To avoid having duplicate images
                        tmp_struc = pmg.Structure(pmg.Lattice(img_struc.lattice.matrix + cur_int*diff_latt), img_struc.species, img_struc.frac_coords + cur_int*diff_frac)
                        write_xml(fle, tmp_struc, trans)
                else:
                    img_struc = pmg.Structure.from_file(NEB_folder + img_str + '/' + FILE_TYPE)
                    write_xml(fle, img_struc, trans)
            else:
                img_struc = pmg.Structure.from_file(NEB_folder + img_str + '/' + FILE_TYPE)
                write_xml(fle, img_struc, trans)

        fle.write('</modeling>\n')

write_NEBrun(FILE_TYPE='POSCAR', trans=True)
