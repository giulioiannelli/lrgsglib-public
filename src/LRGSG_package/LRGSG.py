import networkx as nx
import numpy as np
#
import matplotlib.colors as mplc
import matplotlib.pyplot as plt
#
from farrow_and_ball import *
from numpy.linalg import eigvals
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from scipy.linalg import expm
from scipy.spatial.distance import squareform

from tqdm import tqdm

def get_graph_laplacian(G, is_signed=True):
    if is_signed:
        A = nx.adjacency_matrix(G).toarray()
        D = np.diag(np.abs(A).sum(axis=1))
        L = D-A
    else:
        L = nx.laplacian_matrix(G).todense()
    return L
def entropy(G, steps=1000, is_signed=False, wTresh=1e-10):
    if is_signed:
        L = get_graph_laplacian(G)
        w = eigvals(L)
    else:
        w = nx.laplacian_spectrum(G)
    wSig = w[w>wTresh]
    t1 = np.log10(1. / np.max(wSig))
    t2 = np.log10(10. / np.min(wSig))
    t = np.logspace(t1,t2, int(steps))
    S = np.zeros(len(t))
    VarL = np.zeros(len(t))
    N = G.number_of_nodes()
    
    for i, tau in enumerate(tqdm(t)):      
        Tr = np.nansum(np.exp(-tau * w))
        T1 = np.divide(np.exp(-w*tau), Tr)
        S[i] = -np.nansum(T1 * np.log(T1)) / np.log(N)
        avgRho = np.nansum(np.multiply(w,np.exp(-tau*w))) / Tr
        av2Rho = np.nansum(np.multiply(np.multiply(w,w),np.exp(-tau*w))) / Tr
        VarL[i] = (av2Rho-avgRho*avgRho)
        
    dS = np.log(N) * np.diff(1-S)/np.diff(np.log(t))
    return 1-S, dS, VarL, t