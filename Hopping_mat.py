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
