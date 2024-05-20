# -*- coding: utf-8 -*-
"""
script to calculate density of states in Graphene

@author: Curry
"""

import numpy as np 
from math import pi
import matplotlib.pyplot as plt

import numba 
from numba import jit 
import config as cf


t1 = -1

sq3 = np.sqrt(3)
eplison = 0.0035

numk_dos = cf.numk

font = {'family': "Times New Roman", "weight":"normal", "size":20,}

#We should transform the honeycomb area of the BZ into rectangle lattice
xxx = np.linspace(0, (4*sq3*np.pi)/3, numk_dos, endpoint=True)
yyy = np.linspace((-2*np.pi)/(3), (4*np.pi)/(3),numk_dos, endpoint=True)

#xxx = np.linspace(-np.pi, np.pi, numk)
#yyy = np.linspace(-np.pi, np.pi, numk)

dky = (2*np.pi)/(cf.numk-1)
dkx = (4*np.pi*sq3/3)/(cf.numk-1)

kp = []
for i in xxx:
    for j in yyy:
        kp.append(np.array([i,j]))

#general method
@jit
def Dos():
    #H = np.zeros((2,2), dtype=complex)
    npoints = 500
    E = np.linspace(-3.5, 3.5, npoints, endpoint=True)
    E_dos = np.linspace(-3.5,3.5, npoints, endpoint=True)
    for i in range(npoints):
        E_d = 0
        print("ith is:", i)
        for j in range(numk_dos):
            for w in range(numk_dos):
                k = np.array([xxx[j], yyy[w]]) 
                gk = t1*(np.exp(1.j*k.dot(cf.a1)) + np.exp(1.j*k.dot(cf.a2)) + np.exp(1.j*k.dot(cf.a3)))
                
                e0 = np.real(np.sqrt(gk * gk.conjugate()))
                e1 = np.real(-np.sqrt(gk * gk.conjugate()))
                E_d += eplison/(((E[i]-e0)**2)+eplison**2) + eplison/(((E[i]-e1)**2)+eplison**2) 
        
        E_dos[i] = E_d/(pi)/(numk_dos*numk_dos)
    
    return E, E_dos

#Using Green function to calculate Density of States
@jit
def Dos_green():
    nume = 200 
    
    E = np.linspace(-3.5,3.5,nume, endpoint=True)
    GR = np.zeros((nume,numk_dos*numk_dos),dtype=complex)
    GA = np.zeros((nume,numk_dos*numk_dos),dtype=complex)
   
    for i in range(nume):
        print("step is:", i)
        for j in range(numk_dos*numk_dos):
            e = E[i]
            k = kp[j]
            gk = t1 * (np.exp(1.j*k.dot(cf.a1)) + np.exp(1.j*k.dot(cf.a2)) + np.exp(1.j*k.dot(cf.a3)))
            
            H =np.zeros((2,2), dtype=complex)
            
            H[0,0] = 0 
            H[1,1] = 0
            H[0,1] = gk 
            H[1,0] = gk.conj()
            
            Htmp_r = (e + 1.j*0.01) * np.eye(2, dtype=complex) - H
            Htmp_a = (e - 1.j*0.01) * np.eye(2, dtype=complex) - H
            Hinv_r = np.linalg.inv(Htmp_r)
            Hinv_a = np.linalg.inv(Htmp_a)
            GR[i,j] = np.trace(Hinv_r)
            GA[i,j] = np.trace(Hinv_a)
            
    
    #print("GR is:", GR)
    A = 1.j*(GR-GA)
    dosA = np.real(A)/(2*np.pi) 
    #print("dosA is:", dosA)
    dos = list(map(sum, dosA))
    dos_uniform = [x/(numk_dos*numk_dos) for x in dos]
    
    return E, dos_uniform
    
    
def plot_dos():
    plt.figure(figsize=(10, 8))
    E, E_dos = Dos()
    plt.plot(E, E_dos)
    plt.xlabel("Energy($meV$)", font)
    plt.ylabel("DOS", font)
    plt.xticks(fontproperties = "Times New Roman", fontsize=20)
    plt.yticks(fontproperties = "Times New Roman", fontsize=20)
    plt.grid(True)
    plt.show()
    
if __name__=="__main__":
    plot_dos()