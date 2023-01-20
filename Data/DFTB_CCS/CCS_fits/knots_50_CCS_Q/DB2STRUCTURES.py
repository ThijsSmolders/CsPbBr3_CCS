#!/Users/tjams20/opt/anaconda3/bin/python

import sys
import json
import itertools as it
from collections import OrderedDict, defaultdict
import numpy as np
from ase import Atoms
from ase import io
from ase.calculators.neighborlist import NeighborList
import ase.db as db
from pymatgen import Lattice, Structure
from pymatgen.analysis import ewald

# charge_dict={'Cs':1.0, 'Pb': 2.0, 'Br': -1.0}
# charge_dict={'Cs': 0.835, 'Pb': 0.920, 'Br': -0.585} # Freddie values
# charge_dict={'Cs': 1.0520, 'Pb': 0.9199, 'Br': -0.6573} # Balestra_2020
# charge_dict={'Cs': 0.9495, 'Pb': 1.1547, 'Br':-0.7014} # Bader tjams20
charge_dict={'Cs': 0.812, 'Pb': 1.330, 'Br':-0.714} # Bader tjams20

def pair_dist(atoms, R_c, ch1, ch2):
    ''' This function returns pairwise distances between two types of atoms within a certain cuttoff
    Args:
        R_c (float): Cut off distance(6. Å)
        ch1 (str): Atom species 1
        ch2 (str): Atoms species 2

    Returns:
        A list of distances in bohr
    '''
    if ch1 == ch2:

        nl = NeighborList(atoms.get_number_of_atoms() *
                          [R_c], self_interaction=False, bothways=False)
    else:
        nl = NeighborList(atoms.get_number_of_atoms() *
                          [R_c], self_interaction=False, bothways=True)

    nl.update(atoms)
    distances = []
    for j in range(atoms.get_number_of_atoms()):
        if (atoms.get_chemical_symbols()[j] == ch1):
            indices, offsets = nl.get_neighbors(j)
            for i, offset in zip(indices, offsets):
                if(atoms.get_chemical_symbols()[i] == ch2):
                    distances.append((np.linalg.norm(
                        atoms.positions[i] + np.dot(offset, atoms.get_cell()) - atoms.positions[j])))

    distances.sort()
    r_distances = [round(elem, 6)
                   for elem in distances]
    return r_distances

def main(DFT_data="DFT.db", R_c=5.0, charge_dict=charge_dict):
    """  Function to read files and output structures.json
    
    Args:
        args(list): list of filenames
        R_c (float, optional): Distance cut-off. Defaults to 5.0.
    """
    DFT_DB=db.connect(DFT_data)   

    species = []
    counter=0
    c = OrderedDict()
    d = OrderedDict()
    for row in DFT_DB.select():
        EDFT=row.energy
        key=str(row.key)
        struct=row.toatoms()
        EDFT=DFT_DB.get('key='+key).energy
        counter = counter +1
        print("Reading configuration: "+ str(counter))
        c = OrderedDict()
        dict_species = defaultdict(int)
        struct.charges=[]
        for elem in struct.get_chemical_symbols():
            dict_species[elem] += 1
            struct.charges.append(charge_dict[elem])
        atom_pair = it.combinations_with_replacement(dict_species.keys(), 2)

        lattice= Lattice(struct.get_cell())
        coords= struct.get_scaled_positions()
        ew_struct= Structure(lattice, struct.get_chemical_symbols(), coords, site_properties={"charge" : struct.charges})

        Ew = ewald.EwaldSummation(ew_struct)
        ES_energy=Ew.total_energy

        c['energy_dft'] = EDFT
        c['ewald'] = ES_energy
        c['atoms'] = dict_species
        for (x, y) in atom_pair:
            pair_distances = pair_dist(struct, R_c, x, y)
            c[str(x)+'-'+str(y)] = pair_distances
        d['S'+str(counter+1)] = c
    with open('structures.json', 'w') as f:
        json.dump(d, f, indent=8)

main("DFT.db", charge_dict=charge_dict)


