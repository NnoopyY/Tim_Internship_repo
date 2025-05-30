#!/usr/bin/env python3

import copy
import os
import sys

import ase
import ase.io
import numpy as np
import spglib


element_symbols = [
    'X',
    'H', 'He',  # Period 1
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',  # Period 2
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',  # Period 3
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  # Period 4
    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',  # Period 4
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',  # Period 5
    'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',  # Period 5
    'Cs', 'Ba',  # Period 6
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',  # Lanthanides
    'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',  # Lanthanides
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  # Period 6
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',  # Period 6
    'Fr', 'Ra',  # Period 7
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',  # Actinides
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',  # Actinides
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',  # Period 7
    'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']  # Period 7


def angle_between(v1, v2):
    p = np.dot(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
    return np.rad2deg(np.arccos(np.clip(p, -1.0, 1.0)))


def cellParameters(lattice):
    return (np.linalg.norm(lattice[0]),
            np.linalg.norm(lattice[1]),
            np.linalg.norm(lattice[2]),
            angle_between(lattice[1], lattice[2]),
            angle_between(lattice[0], lattice[2]),
            angle_between(lattice[0], lattice[1]))


def writeCIF(cell, prec, basename):
    sg = spglib.get_spacegroup(cell, symprec=prec)
    sg, sgid = sg.split(' (')
    sgid = sgid.replace(')', '')

    lattice, scaled_positions, numbers = spglib.refine_cell(cell, prec)
    ncell = (lattice, scaled_positions, numbers)
    sym = spglib.get_symmetry(ncell, prec)
    uniques = np.unique(sym['equivalent_atoms'])
    a, b, c, alpha, beta, gamma = cellParameters(lattice)

    f = open((basename + '_' + sgid + '.cif').replace('/', '|'), 'w')

    f.write(f'# Symmetrized structure with precision = {prec}\n')
    f.write(f'data_{basename}{sg}\n'.replace(' ', '_'))
    f.write('_cell_length_a                  %.7g\n' % a)
    f.write('_cell_length_b                  %.7g\n' % b)
    f.write('_cell_length_c                  %.7g\n' % c)
    f.write('_cell_angle_alpha               %.7g\n' % alpha)
    f.write('_cell_angle_beta                %.7g\n' % beta)
    f.write('_cell_angle_gamma               %.7g\n' % gamma)
    f.write("_symmetry_space_group_name_H-M  '" + sg + "'\n")
    f.write('_symmetry_Int_Tables_number     ' + str(sgid) + '\n')

    f.write('''
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
''')

    for pos, at in zip(scaled_positions[uniques], numbers[uniques]):
        sym = element_symbols[at]
        f.write('%s %s 1.0 %.7f %.7f %.7f\n' % (sym, sym, pos[0], pos[1], pos[2]))

    f.close()


# Main
if len(sys.argv) != 2:
    print(f'Usage: {os.path.basename(sys.argv[0])} <file.cif|POSCAR>')
    sys.exit(1)

filename = sys.argv[1]

try:
    structure = ase.io.read(filename)
    lattice = structure.get_cell()
    positions = structure.get_scaled_positions()
    numbers = structure.get_atomic_numbers()
    cell = (lattice, positions, numbers)

except Exception as e:
    print('Could not read structure from file: ' + filename)
    print('Error message is: ' + str(e))
    sys.exit(1)

basename = os.path.splitext(os.path.basename(filename))[0]

t = [1, 2, 3, 5, 7]
precs = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6] + [i * j for i in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1] for j in t]

old = ''
print('# Tolerance\tSpace group')

for prec in precs:
    s = spglib.get_spacegroup(cell, symprec=prec)
    if s != old:
        print(f'{prec}\t\t{s}')
        writeCIF(cell, prec, basename)
        old = s

sys.exit(0)
