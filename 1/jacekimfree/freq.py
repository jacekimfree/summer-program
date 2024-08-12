#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../../extra-files')
from masses import *
sys.path.append('../../0/jacekimfree')
from molecule import Molecule

import numpy as np
import math

def Frequencies(filename, units='Angstrom'):
    
    # read in the Hessian matrix
    hess_list = []
    with open(filename,'r') as fn:
        for line in fn:
            hess_list.append(line.split())      # split and append each line
    hess = np.array(hess_list).astype(float)
    # print(hess)
    
    # build Molecule object from molecule.xyz
    m = Molecule('../../extra-files/molecule.xyz')
    m.to_angstrom()
    # print(m.xyz_string())
    
    # build mass-weighted Hessian matrix
    mw_hess = hess.copy()
    for i in range(len(m.labels)):
        for j in range(len(m.labels)):
            m1 = get_mass(m.labels[i])
            m2 = get_mass(m.labels[j])
            for k in range(3):
                for l in range(3):
                    mw_hess[((3*i+k))][(3*j+l)] = hess[((3*i+k))][(3*j+l)] / ((m1 * m2)**0.5)
    # print(mw_hess)
    
    # compute evals(k) and evecs(mw_q) of mw_hess
    k, mw_q = np.linalg.eigh(mw_hess)
    # print(k)
    # print(q)
    
    # un-mass-weight the evecs(mw_q) to get normal coordinates
    mass_mat = np.zeros(mw_hess.shape)          # mass_mat should end up becoming diagonal
    mass_list = []
    for i in m.labels:
        for j in range(3):
            mass_list.append(get_mass(i))
    for i in range(len(mass_list)):
        mass_mat[i][i] = mass_list[i]**(-0.5)
    # print(mass_mat)
    q = np.dot(mass_mat,mw_q)
    q = np.transpose(q)
    print(q)
    
    # define unit conversion constants
    n1 = 4.35974417e-18/((5.2917721092e-11**2)*1.6605389e-27)
    n2 = (2 * math.pi)
    c = 2.99792458e10
    
    # convert frequencies from a.u. to cm-1
    k = [complex(i * n1,0) for i in k]
    nu = np.sqrt(k)
    nu = [i / n2 for i in nu]
    nu = np.array([i / c for i in nu])
    # print(nu)
    
    # formatting output
    output_string = ''
    for i in range(len(nu)):
        output_string += str(len(m.labels)) + '\n'
        if nu[i] != 0:
            if nu[i].imag != 0:
                output_string += str(np.round(nu[i].imag,2)) + 'i cm^-1\n'
            elif nu[i].real != 0:
                output_string += str(np.round(nu[i].real,2)) + '  cm^-1\n'
            for j in range(len(m.labels)):
                output_string += '{:2}{:>15.10f}{:>15.10f}{:>15.10f}{:>15.10f}{:>15.10f}{:>15.10f}\n'.format(str(m.labels[j]), m.geom[j][0], m.geom[j][1], m.geom[j][2], q[i,3*j], q[i,3*j+1], q[i,3*j+2])
        output_string += '\n'
    return output_string
'''
Test code
'''

if __name__ == "__main__":
    m = Frequencies('../../extra-files/hessian.dat')
    f = open('normal_modes.xyz','w')
    f.write(m)
    f.close()