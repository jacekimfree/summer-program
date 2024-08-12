#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../../extra-files')
from masses import *

import numpy as np
from copy import deepcopy

class Molecule():
    
    def __init__(self, filename):
        self.units = []
        self.natom = []
        self.labels = []
        self.masses = []
        self.charges = []
        self.geom = []
        
        lines = []
        with open(filename,'r') as fn:
            for line in fn:
                lines.append(line.strip())
        # print(lines)
        
        # append values to corresponding variables
        self.units.append(lines[1])             # units
        self.natom.append(lines[0])             # natoms
        
        geom = []                               # labels, geom
        for line in lines[2:]:
            atom, x, y, z = line.split()
            self.labels.append(atom)
            geom.append([float(x), float(y), float(z)])
        self.geom = np.array(geom)
        # print(geom)
        
        for atom in self.labels:                # masses, charges
            self.masses.append(get_mass(atom))
            self.charges.append(get_charge(atom))
    
    # convert unit to Bohr
    def to_bohr(self):
        if self.units[0] != 'Bohr':
            self.geom *= 1.88971616463
            self.units[0] = 'Bohr'
        else:
            print('Already in Bohr')
    
    # convert unit to Angstrom
    def to_angstrom(self):
        if self.units[0] != 'Angstrom':
            self.geom /= 1.88971616463
            self.units[0] = 'Angstrom'
        else:
            print('Already in Angstrom')
    
    # print Molecule object
    def xyz_string(self, option=0):
        if option == 0:
            string = str(self.natom) + '\n'
        else:
            string = '{}\n{}'.format(str(self.natom[0]),str(self.units[0]))
        for i in range(len(self.labels)):
            string += '\n{:<10}{:<20}{:<20}{:<20}'.format(str(self.labels[i]), str(self.geom[i][0]), str(self.geom[i][1]), str(self.geom[i][2]))
        return string
        
'''
Test code
'''

if __name__ == "__main__":
    m = Molecule('../../extra-files/molecule.xyz')
    
    print(m.xyz_string())
    
    m.to_angstrom()
    print(m.xyz_string())
    
    m.to_bohr()
    print(m.xyz_string())