#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
sys.path.append('../../extra-files')
from masses import *
sys.path.append('../../0/jacekimfree')
from molecule import Molecule
sys.path.append('../../1/jacekimfree')
from freq import Frequencies

import numpy as np

class Hessian(object):
    def __init__(self, mol, disp_size=0.005, directory='DISPS', command='psi4', energy_prefix='@DF-RHF Final Energy:'):
        self.mol = mol
        self.disp_size = disp_size
        self.directory = directory
        self.command = command
        self.energy_prefix = energy_prefix
        self.cwd = os.getcwd()
        self.hes_directory = '{}/{}'.format(self.cwd,self.directory)
        
        # unify units to angstrom
        if self.mol.units != 'Angstrom':
            mol.to_angstrom()
        
        with open('template.dat','r') as fn:
            self.template = fn.read()
        
    def write_file(self, mol, filename):
        with option(filename,'w') as fn:
            fn.write(self.template.format(mol.xyz_string(option=1)))
            

'''
test
'''

if __name__ =='__main__':
    