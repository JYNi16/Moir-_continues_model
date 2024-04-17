# -*- coding: utf-8 -*-
"""
The Hamiltonian of the twist trilayer Graphene with magic angle xxx

@author: Curry
"""

from numpy import *
import numpy as np
from config import *
from mat_ import * 

T1, T1_h = T_mat(1), np.array(np.matrix(T_mat(1)).H)
T2, T2_h = T_mat(2), np.array(np.matrix(T_mat(2)).H)
T3, T3_h = T_mat(3), np.array(np.matrix(T_mat(3)).H)

lattN = (2*N+1)*(2*N+1)
invL, L = Lattice(N)

print("invL is:", invL)
print("L is:", L)

def Hamiltonian(k):
    kx, ky = k
    
    H = array(zeros((6*lattN, 6*lattN)), dtype=complex)
    
    for i in np.arange(lattN):
        #diagonal term
        m = L[i, 0]
        n = L[i, 1]
        
        ax = kx - valley*K1[0] + m*b1m[0] + n*b2m[0]
        ay = ky - valley*K1[1] + m*b1m[1] + n*b2m[1]     
        qx, qy = taylor_exp(ax, ay, theta), taylor_exp(ay, ax, -theta)
        
        #diagonal term of the first layer 
        H[2*i, 2*i+1] = hv * (valley*qx - 1.j*qy)
        H[2*i+1, 2*i] = hv * (valley*qx + 1.j*qy)
        
        #diagonal term of the third layer
        H[2*i+4*lattN, 2*i+4*lattN+1] = hv * (valley*qx - 1.j*qy)
        H[2*i+4*lattN+1, 2*i+4*lattN] = hv * (valley*qx + 1.j*qy)
        
        ax2 = kx - valley*K2[0] + m*b1m[0] + n*b2m[0] 
        ay2 = ky - valley*K2[1] + m*b1m[1] + n*b2m[1] 
        qx2, qy2 = taylor_exp(ax2, ay2, -theta), taylor_exp(ay2, ax2, theta)
        
        #diagonal term of the second layer 
        H[2*(i+lattN), 2*(i+lattN)+1] = hv * (valley*qx2 - 1.j*qy2)
        H[2*(i+lattN)+1, 2*(i+lattN)] = hv * (valley*qx2 + 1.j*qy2)

        #off-diagonal term when m1-m2==0 and n1-n2==0
        #lower left part (1st-layer -> 2nd-layer)
        j = i + lattN
        H[2*j, 2*i]     = T1_h[0, 0]
        H[2*j, 2*i+1]   = T1_h[0, 1]
        H[2*j+1, 2*i]   = T1_h[1, 0]
        H[2*j+1, 2*i+1] = T1_h[1, 1]
        
        #upper right part
        H[2*i, 2*j]     = T1[0, 0]
        H[2*i, 2*j+1]   = T1[0, 1]
        H[2*i+1, 2*j]   = T1[1, 0]
        H[2*i+1, 2*j+1] = T1[1, 1]
        
        #upper right part (2nd-layer -> 3rd-layer)
        H[2*j, 2*i+4*lattN]     = T1[0, 0]
        H[2*j, 2*i+4*lattN+1]   = T1[0, 1]
        H[2*j+1, 2*i+4*lattN]   = T1[1, 0]
        H[2*j+1, 2*i+4*lattN+1] = T1[1, 1]     
            
        #lower left part
        H[2*i+4*lattN, 2*j]   = T1[0, 0]
        H[2*i+4*lattN, 2*j+1]   = T1[0, 1]
        H[2*i+4*lattN+1, 2*j]   = T1[1, 0]
        H[2*i+4*lattN+1, 2*j+1]   = T1[1, 1]    
        
        #off-diagonal term when m1-m2==0, n1-n2==-1*valley
        if (n != valley*N):
            j = invL[m+N, n+valley*1+N] + lattN
            
            #lower left part (1st-layer -> 2nd-layer)
            H[2*j, 2*i]     = T2_h[0, 0]
            H[2*j, 2*i+1]   = T2_h[0, 1]
            H[2*j+1, 2*i]   = T2_h[1, 0]
            H[2*j+1, 2*i+1] = T2_h[1, 1]
            
            #upper right part
            H[2*i, 2*j]     = T2[0, 0]
            H[2*i, 2*j+1]   = T2[0, 1]
            H[2*i+1, 2*j]   = T2[1, 0]
            H[2*i+1, 2*j+1] = T2[1, 1]
            
            #upper right part (2nd-layer -> 3rd-layer)
            H[2*j, 2*i+4*lattN]     = T2_h[0, 0]
            H[2*j, 2*i+4*lattN+1]   = T2_h[0, 1]
            H[2*j+1, 2*i+4*lattN]   = T2_h[1, 0]
            H[2*j+1, 2*i+4*lattN+1] = T2_h[1, 1]          
                
            #lower left part
            H[2*i+4*lattN, 2*j]   = T2[0, 0]
            H[2*i+4*lattN, 2*j+1]   = T2[0, 1]
            H[2*i+4*lattN+1, 2*j]   = T2[1, 0]
            H[2*i+4*lattN+1, 2*j+1]   = T2[1, 1]        
            
        #off-diagonal term when m1-m2==1*valley, n1-n2==0
        if (m != -valley*N):
            j = invL[m-valley*1+N, n+N] + lattN
            
            #lower left part 1st-layer -> 2nd-layer)
            H[2*j, 2*i]     = T3_h[0, 0]
            H[2*j, 2*i+1]   = T3_h[0, 1]
            H[2*j+1, 2*i]   = T3_h[1, 0]
            H[2*j+1, 2*i+1] = T3_h[1, 1]
            
            #upper right part
            H[2*i, 2*j]     = T3[0, 0]
            H[2*i, 2*j+1]   = T3[0, 1]
            H[2*i+1, 2*j]   = T3[1, 0]
            H[2*i+1, 2*j+1] = T3[1, 1]
            
            #upper right part (2nd-layer -> 3rd-layer)
            H[2*j, 2*i+4*lattN]     = T3_h[0, 0]
            H[2*j, 2*i+4*lattN+1]   = T3_h[0, 1]
            H[2*j+1, 2*i+4*lattN]   = T3_h[1, 0]
            H[2*j+1, 2*i+4*lattN+1] = T3_h[1, 1]           
                
            #lower left part
            H[2*i+4*lattN, 2*j]   = T3[0, 0]
            H[2*i+4*lattN, 2*j+1]   = T3[0, 1]
            H[2*i+4*lattN+1, 2*j]   = T3[1, 0]
            H[2*i+4*lattN+1, 2*j+1]   = T3[1, 1]
            
    return H

#set_printoptiopns in printing
#np.set_printoptions(threshold=np.inf)
#H = Hamiltonian(G)
#print(H[0][:])
