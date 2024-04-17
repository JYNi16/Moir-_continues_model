# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:32:23 2024

@author: 26526
"""

from numpy import *
import numpy as np
from config import *

#define the interlayer hopping term matrix Tth
def T_mat(i):    
    return w0*s0 + w1*(cos((2*pi/3)*(i-1))*sx + valley*sin((2*pi/3)*(i-1))*sy)


def taylor_exp(ax, ay, theta):
    
    return cos(theta/2)*ax + sin(theta/2) * ay

#define the data array with N cutoff in k-space
def Lattice(n):
    # n is the cutoff label
    Lmn = []
    idx_L = np.zeros((2*N+1, 2*N+1), int)
    count = 0
    for i in np.arange(-n, n+1):
        for j in np.arange(-n, n+1):
            Lmn.append([i, j])
            idx_L[i+n, j+n] = count
            count = count + 1
    
    return idx_L, np.array(Lmn)
