# -*- coding: utf-8 -*-
"""
Script for ploting band for MATBG ...

@author: Curry
"""

from config import * 
import numpy as np 
import matplotlib.pyplot as plt
from MATBG_Ham import Hamiltonian

def H(k):
    e = np.sort(np.real(np.linalg.eig(Hamiltonian(k))[0]))
    return e

def Dist(r1, r2):
    return np.linalg.norm(r1-r2)

#define the k-path
def k_sym_path():
    
    kmk = np.linspace(K1,K2,npoints)
    kkg = np.linspace(K2,G,npoints)
    kgm = np.linspace(G,M,npoints)
    kmkm = np.linspace(M,K1,npoints)
    
    k_point_path = [kmk, kkg, kgm, kmkm]
    
    lkmk = Dist(K1,K2)
    lkg = Dist(K2,G)
    lgm = Dist(G,M)
    lmkm = Dist(M,K1)

    lk = np.linspace(0, 1, npoints)
    xkmk = lkmk * lk 
    xkg = lkg * lk + xkmk[-1]
    xgm = lgm * lk + xkg[-1]
    xmkm = lmkm * lk + xgm[-1]
    

    kpath = np.concatenate((xkmk, xkg, xgm, xmkm), axis=0)
    
    Node = [0, xkmk[-1], xkg[-1], xgm[-1], xmkm[-1]]
    k_path = np.concatenate((xkmk, xkg, xgm, xmkm), axis=0)
    
    return k_point_path, k_path, Node

def band_post():
    k_point_path, k_path, Node = k_sym_path()
    E_band = []
    
    for i in range(len(k_point_path)):
        E_values = np.array(list(map(H, k_point_path[i])))
        #print("E_values is:", E_values)
        E_band.append(E_values)
    
    return E_band

def plot_band(): 
    k_sym_label =  [r"$K_{m}$",  r"$K^{\prime}_{m}$", r"$\Gamma_{m}$", r"$M_{m}$",  r"$K_{m}$"]
    font = {'family': "Times New Roman", "weight":"normal", "size":20,}
    k_point_path, k_path, Node = k_sym_path()
    E_band = np.array(band_post())
    shape = E_band.shape
    print("E_band.shape is:", shape) 
        
    plt.figure(1, figsize=(8,8))
    if len(shape) < 2:
        eig = np.hstack(tuple(E_band))
        plt.plot(k_path, eig)
        #plt.xticks(Node,Node_label)
        plt.show()   
        return 
    
    for i in range(shape[-1]):
        eig_test = [] 
        for j in range(shape[0]):
            eig_test.append(E_band[j][:,i])
        eig = np.hstack(tuple(eig_test))
        plt.plot(k_path, eig, "black", linewidth=2)   
    
    plt.xlim(0, k_path[-1])
    plt.ylim(-250,250)
    plt.xticks(Node, k_sym_label) 
    #plt.xlabel("$K$-points", font)
    plt.ylabel("Energy($meV$)", font)
    font_txt = {'style': "normal", "weight":"normal", "size":20, 'family': "Times New Roman"}
    plt.xticks(fontproperties = "Times New Roman", fontsize=20)
    plt.yticks(fontproperties = "Times New Roman", fontsize=20)
    #plt.text(0.2,6.1, "(a)", fontsize=20, style= "Times New Roman")
    #plt.text(4.5,5, "$\mathregular{\Delta J_2 / D = 0.1}$", fontdict = font_txt)
    title = r"Band of TBG with magic angle of {}$^\degree$ and $w_1 / w_0 = {}$".format(theta_v, r1)
    plt.title(title,loc = "center",fontdict={"size":"xx-large","color":"black", "family":"Times New Roman"})
    
    plt.savefig("./figure/MATBG_{}_{}.png".format(theta_v, r1), dpi=500)
    
    plt.show()


if __name__=="__main__":
    plot_band()