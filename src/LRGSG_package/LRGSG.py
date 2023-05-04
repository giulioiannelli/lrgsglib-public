#
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import random
#
import networkx as nx
import numpy as np
#
import matplotlib.colors as mplc
import matplotlib.pyplot as plt
import scipy.special
#
from farrow_and_ball import *
from numpy.linalg import eigvals
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from scipy.linalg import expm, fractional_matrix_power, ishermitian
from scipy.sparse.linalg import eigs, eigsh, ArpackNoConvergence
from scipy.spatial.distance import squareform
#
from tqdm import tqdm
#
ePDF = ".pdf"
eTXT = ".txt"
#
datPath_lmin = "data/lmin/"
setPath_ERp = "conf/ERp/"
pltPath_Sm1C = "plot/Sm1_and_C/"
#
def smallest_prob_for_erconn(N, pstart=0.1, halving=0.8, testset=100):
    p = pstart
    while True:
        areallconn = []
        for _ in range(testset):
            G = nx.erdos_renyi_graph(N, p, seed=None, directed=False)
            areallconn.append(nx.is_connected(G))
        if not all(areallconn):
            lowest_prob = p
            break
        p *= halving
    return lowest_prob
def logpsace_prob_erconn(N, pHigh=0.5, pstart=0.1, halving=0.8, testset=100):
    lowest_prob = smallest_prob_for_erconn(N, pstart, halving, testset)
    lsp = np.logspace(np.log10(lowest_prob), np.log10(pHigh), num=10)
    return lsp
def get_graph_lapl(G, is_signed=True):
    if is_signed:
        lapl = slaplacian_matrix(G).asfptype()
    else:
        lapl = nx.laplacian_matrix(G).asfptype()
    return lapl

def get_graph_lspectrum(G, is_signed=False):
    if is_signed:
        A = nx.adjacency_matrix(G).toarray()
        D = np.diag(np.abs(A).sum(axis=1))
        L = D-A
        w = eigvals(L)
    else:
        L = nx.laplacian_matrix(G).todense()
        w = nx.laplacian_spectrum(G)
    return L, w

def slaplacian_matrix(G, nodelist=None, weight="weight"):
    """Returns the signed Laplacian matrix of G.

    The graph Laplacian is the matrix L = D - A, where
    A is the adjacency matrix and D is the diagonal matrix of node degrees.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : list, optional
       The rows and columns are ordered according to the nodes in nodelist.
       If nodelist is None, then the ordering is produced by G.nodes().

    weight : string or None, optional (default='weight')
       The edge data key used to compute each value in the matrix.
       If None, then each edge has weight 1.

    Returns
    -------
    L : SciPy sparse array
      The Laplacian matrix of G.

    Notes
    -----
    For MultiGraph, the edges weights are summed.

    See Also
    --------
    to_numpy_array
    normalized_laplacian_matrix
    laplacian_spectrum

    Examples
    --------
    For graphs with multiple connected components, L is permutation-similar
    to a block diagonal matrix where each block is the respective Laplacian
    matrix for each component.

    >>> G = nx.Graph([(1, 2), (2, 3), (4, 5)])
    >>> print(nx.laplacian_matrix(G).toarray())
    [[ 1 -1  0  0  0]
     [-1  2 -1  0  0]
     [ 0 -1  1  0  0]
     [ 0  0  0  1 -1]
     [ 0  0  0 -1  1]]

    """
    import scipy as sp
    import scipy.sparse  # call as sp.sparse

    if nodelist is None:
        nodelist = list(G)
    A = nx.to_scipy_sparse_array(G, nodelist=nodelist, weight=weight, format="csr")
    n, m = A.shape
    D = sp.sparse.csr_array(sp.sparse.spdiags(np.abs(A.sum(axis=1)), 0, m, n, format="csr"))
    return D - A

def get_graph_lspectrum_rw(G, is_signed=False):
    A = nx.adjacency_matrix(G).toarray()
    D = np.diag(np.abs(A).sum(axis=1))
    L = np.eye(D.shape[0]) - fractional_matrix_power(D, -.5)@A@fractional_matrix_power(D, -.5)
    if is_signed:
        w = eigvals(L)
    else:
        w = nx.laplacian_spectrum(G)
    return L, w

def entropy(G, steps=600, is_signed=False, wTresh=1e-15, t1=-2, t2=5):
    N = G.number_of_nodes()

    L, w = get_graph_lspectrum_rw(G, is_signed=is_signed)
    wSig = w[w>wTresh]
    
    t = np.logspace(t1,t2, int(steps))
    S = np.zeros(len(t))
    VarL = np.zeros(len(t))
    
    for i, tau in enumerate(t):
        rhoTr = np.exp(-tau * w)
        Tr = np.nansum(rhoTr)
        rho = np.divide(rhoTr, Tr)
        S[i] = -np.nansum(rho * np.log(rho)) / np.log(N)
        avgrho = np.nansum(np.multiply(w, rhoTr)) / Tr
        av2rho = np.nansum(np.multiply(np.multiply(w,w), rhoTr)) / Tr
        VarL[i] = (av2rho - avgrho**2)
        
    dS = np.log(N) * np.diff(1-S)/np.diff(np.log(t))
    return 1-S, dS, VarL, t

def averaged_Sm1(t1Sm1Avg):
    lenLs = t1Sm1Avg[0][0].__len__()
    avgSm1 = np.zeros(lenLs)
    all_t1 = np.concatenate([i[0] for i in t1Sm1Avg])
    commonLs = np.logspace(np.log10(min(all_t1)), np.log10(max(all_t1)), lenLs)
    digitizedLs = []
    for tS in t1Sm1Avg:
        digiTmp = np.digitize(tS[0], bins=commonLs)-1
        np.add.at(avgSm1, digiTmp, tS[1])
        digitizedLs.extend(digiTmp)
    unique, counts = np.unique(digitizedLs, return_counts=True)
    # if avgSm1.shape != counts.shape:
    #     counts = np.delete(counts, 0)
    #     print(avgSm1, unique, counts)
    try:
        avgSm1 /= counts
    except ValueError:
        print("counts has not the same dimensions of avgSm1.")
        print(counts.shape, avgSm1.shape)
    return commonLs, avgSm1

def lapl_dists(L, tau=1e-2, is_signed=False):
    num = expm((-tau*L))
    rho = num/np.trace(num)
    Trho = np.copy(1. / rho)
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0)
    if is_signed:
        old_d = squareform(Trho)
        dists = np.sqrt(np.max(old_d) - old_d)
    else:
        dists = squareform(Trho)
    return dists

def MakeLinkageMatrix(G, tau=1e-2, is_signed=False, method="ward"):
    L, w = get_graph_lspectrum(G, is_signed)
    dists = lapl_dists(L, tau, is_signed)
    linkage_matrix1 = linkage(dists, method=method)
    tmax = linkage_matrix1[::, 2][-1]
    linkage_matrix = linkage(dists/tmax, method=method)
    return linkage_matrix, w