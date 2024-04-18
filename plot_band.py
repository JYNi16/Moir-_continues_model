# -*- coding: utf-8 -*-
"""
Module to plot moire band for MATBG/MATTG/MATTMD ...

@author: Curry
"""

from config import * 
import numpy as np 
import matplotlib.pyplot as plt
from MATBG_Ham import Hamiltonian
from k_sym_gen import *

def H(k):
    e = np.sort(np.real(np.linalg.eigh(Hamiltonian(k))[0]))
    return e

def band_post():
    syms = [K2, K1, G, M, K2]
    #syms = [K2, K1, G, M, K2]
    k_point_path, k_path, Node = k_path_sym_gen(syms)
    E_band = []
    
    for i in range(len(k_point_path)):
        E_values = np.array(list(map(H, k_point_path[i])))
        #print("E_values is:", E_values)
        if (len(E_values.shape) < 2):
            E_band.append((np.reshape(E_values,[E_values.shape[0], -1])))
        else:
            E_band.append(E_values)
    
    return np.array(E_band), k_point_path, k_path, Node

def plot_band(): 
    k_sym_label =  [r"$K^{\prime}_{m}$",r"$K_{m}$", r"$\Gamma_{m}$", r"$M_{m}$",  r"$K^{\prime}_{m}$"]
    #k_sym_label =  [r"$K_{m}$", r"$K^{\prime}_{m}$",  r"$\Gamma_{m}$", r"$M_{m}$", r"$K_{m}$"]
    #k_sym_label =  [r"$K_{2}$",r"$\Gamma$", r"$K_{1}$", r"$K_{2}$", r"$K_{2}^{\prime}$"]
    font = {'family': "Times New Roman", "weight":"normal", "size":28,}
    E_band, k_point_path, k_path, Node= band_post()
    shape = E_band.shape
    print("E_band.shape is:", shape) 
        
    plt.figure(1, figsize=(10,8))
    
    for i in range(shape[-1]):
        eig_test = [] 
        for j in range(shape[0]):
            eig_test.append(E_band[j][:,i])
        eig = np.hstack(tuple(eig_test))
        plt.plot(k_path, eig, linewidth=3)   
    
    plt.xlim(0, k_path[-1])
    plt.ylim(-1200,1200)
    plt.xticks(Node, k_sym_label) 
    #plt.xlabel("$K$-points", font)
    plt.ylabel("Energy($meV$)", font)
    font_txt = {'style': "normal", "weight":"normal", "size":20, 'family': "Times New Roman"}
    plt.xticks(fontproperties = "Times New Roman", fontsize=20)
    plt.yticks(fontproperties = "Times New Roman", fontsize=20)
    #plt.text(0.2,6.1, "(a)", fontsize=20, style= "Times New Roman")
    #plt.text(4.5,5, "$\mathregular{\Delta J_2 / D = 0.1}$", fontdict = font_txt)
    title = r"Band of TBG with magic angle of {}$^\degree$ $w_1$ = {} and $w_1 / w_0 = {}$".format(theta_v, w1, r1)
    plt.title(title,loc = "center",fontdict={"size":"xx-large","color":"black", "family":"Times New Roman"})
    
    plt.savefig("./figure/MATBG_{}_{}_{}.png".format(theta_v, r1, w1), dpi=500)
    
    plt.show()


if __name__=="__main__":
    plot_band()