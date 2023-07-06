#
import argparse
import random
import re
import os
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#
import networkx as nx
import numpy as np
# import pyfftlog
#
import matplotlib.animation as animation
import matplotlib.colors as mplc
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import scipy.special
import scipy.sparse as scsp
#
from itertools import product
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Circle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from networkx.classes.graph import Graph
from numpy import ndarray
from numpy.linalg import eigvals, eigvalsh
from scipy import stats
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from scipy.interpolate import make_interp_spline
from scipy.linalg import expm, fractional_matrix_power, ishermitian
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
from scipy.sparse import csr_array
from scipy.sparse.linalg import eigs, eigsh, ArpackNoConvergence
from scipy.spatial.distance import squareform
#
from tqdm import tqdm
#
MAX_DIGITS_ROUND_SIGFIG = 18
DEFAULT_ENTROPY_STEPS = 1000
DEFAULT_ENTROPY_LEXPONENT = -3
DEFAULT_ENTROPY_HEXPONENT = 6
DEFAULT_NUNMBER_AVERAGES = 1000
DEFAULT_SPIKE_THRESHOLD = 0.05
DEFAULT_MAX_THRESHOLD = 2 * DEFAULT_SPIKE_THRESHOLD
#
ePDF = ".pdf"
eTXT = ".txt"
eBIN = ".bin"
#
datPath_lminl2d = "data/lmin_l2d/"
datPath_l2d_sq = "data/l2d_sq/"
lambdaPath_l2d = lambda geometry : f"l2d_{geometry}/"
pltPath_l2d = lambda geometry : f"data/plot/{lambdaPath_l2d(geometry)}"
datPath_l2d = lambda geometry : f"data/{lambdaPath_l2d(geometry)}"
setPath_ERp = "conf/ERp/"
pltPath_Sm1C = "plot/Sm1_and_C/"
#
pflip_fmt = '.3g'
#
def line(x, a, b):
    return a * x + b
# basic math functions
def dv(f_x: ndarray, x: ndarray = None) -> ndarray:
    """Compute the computational derivative of an array with respect to another.

    Parameters
    ----------
    f_x : ndarray
        A N dimensional NumPy array where the outer dimension is the one where
        the `np.diff` method is applied.
    x : ndarray, optional
        An array with the same dimensions of the outer axis of `f_x`. If the
        independent variable array is not passed the one assumed is the
        `range(0, len(f(x)))`.
    
    Returns
    -------
    df_dx : ndarray 

    """
    if x is None:
        x = np.linspace(0, f_x.shape[-1], num=f_x.shape[-1])
    df_dx = np.diff(f_x, axis=-1) / np.diff(x)
    return df_dx
# networkx graph related functions
def slaplacian_matrix(G: Graph, nodelist: list = None, weight: str = "weight"
                      ) -> csr_array:
    """Returns the signed Laplacian matrix of G.

    The graph Laplacian is the matrix L = |D| - A, where
    A is the adjacency matrix and |D| is the diagonal matrix of absolute values
    of node degrees.

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
    """
    if nodelist is None:
        nodelist = list(G)
    A = nx.to_scipy_sparse_array(G, nodelist=nodelist, weight=weight,
                                 format="csr")
    D = scsp.csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 0, *A.shape,
                                    format="csr"))
    return D - A
#
def flip_random_fract_edges(G: Graph, p: float):
    """Flips a fraction p of edges (+1 to -1) of a graph G.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    p : float
       The fraction in [0, 1] of edges to be flipped.
    """
    eset = G.edges()
    nedges = len(eset)
    noflip = int(p*nedges)
    #
    if noflip < 1:
        print('p too small')
        exit()
    rndsmpl = random.sample(range(nedges), noflip)
    #
    # all_weights = {e: 1 for e in eset}
    neg_weights = {e: -1 for i,e in enumerate(eset) if i in rndsmpl}
    #
    nx.set_edge_attributes(G, values=1, name='weight')
    nx.set_edge_attributes(G, values=neg_weights, name='weight')
#
def flip_one_2dgraph(G, coord1, coord2):
    neg_weights = {(coord1, coord2): -1}
    #
    nx.set_edge_attributes(G, values=neg_weights, name='weight')
def get_kth_order_neighbours(G, node, order):
    node = 0
    # Set the desired order
    order = 1
    # Get the kth-order neighbors
    neighbor_dict = nx.single_source_shortest_path_length(G, node, cutoff=order)
    # Extract the nodes at the kth order
    kth_order_neighbors = [n for n, d in neighbor_dict.items() if d == order]
    return kth_order_neighbors
#
class NflipError(Exception):
    pass
# renormalization group for heterogenous network functions
class SignedLaplacianAnalysis:
    slspectrum = None
    Sm1 = None
    VarL = None
    Cspe = None
    taumax = None
    taumax0 = None
    #
    def __init__(self, system, pflip: float = None, is_signed: bool = True,
                 steps: int = DEFAULT_ENTROPY_STEPS, t1: float = -2,
                 t2: float = 5, maxThresh: float = DEFAULT_MAX_THRESHOLD, 
                 nreplica: int = 0, initCond: str = 'gauss_1') -> None:
        self.system = system
        self.is_signed = is_signed
        #
        self.nreplica = nreplica
        self.steps = steps
        self.t1 = t1
        self.t2 = t2
        self.maxThresh = maxThresh
        self.initCond = initCond
        #
        self.pflip = pflip
        if pflip is not None:
            self.nflip = int(self.pflip*self.system.Ne)
            self.randsample = random.sample(range(self.system.Ne), 
                                            self.nflip)
        else:
            self.nflip = 0
            self.randsample = None
        
    #
    def init_weights(self):
        nx.set_edge_attributes(self.system.G, values=1, name='weight')
    #
    def upd_graph_matrices(self):
        self.Adj = self.adjacency_matrix(self.system.G)
        self.Deg = self.degree_matrix(self.Adj)
        self.sDeg = self.absolute_degree_matrix(self.Adj)
        self.Lap = self.laplacian_csr()
        self.sLp = self.signed_laplacian_csr()
    #
    def check_pflip(self):
        if self.nflip < 1:
            raise NflipError("""The probability of flipping an edge times the 
                             number of edges is < 1, then no edges would be
                             flipped. Skipping the analysis for this value.""")
    #
    def flip_random_fract_edges(self):
        """Flips a fraction p of edges (+1 to -1) of a graph G.
        """
        
        #
        try:
            self.check_pflip()
        except NflipError:
            return None
        neg_weights = {e: -1 for i,e in enumerate(self.system.eset)
                       if i in self.randsample}
        nx.set_edge_attributes(self.system.G, values=neg_weights, name='weight')
        self.upd_graph_matrices()
    #
    def flip_sel_edges(self, neg_weights_dict):
        """Flips a specific edges of a graph G.
        """
        #
        nx.set_edge_attributes(self.system.G, values=neg_weights_dict, 
                               name='weight')
        self.upd_graph_matrices()
    #
    def adjacency_matrix(self, G: Graph,
                         nodelist: list = None, weight: str = "weight"):
        if nodelist is None: nodelist = list(G)
        return nx.to_scipy_sparse_array(G, nodelist=nodelist, 
                                        weight=weight, format="csr")
    #
    def degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(scsp.spdiags(A.sum(axis=1), 
                                        0, *A.shape, format="csr"))
    #
    def absolute_degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 
                                        0, *A.shape, format="csr"))
    #
    def laplacian_csr(self) -> csr_array:
        """Returns the signed Laplacian matrix of G.
        The graph Laplacian is the matrix L = D - A, where
        A is the adjacency matrix and D is the diagonal matrix of node degrees

        Returns
        -------
        L : SciPy sparse array
        The Laplacian matrix of G.
        """
        return self.Deg - self.Adj
    #
    def signed_laplacian_csr(self) -> csr_array:
        """Returns the signed Laplacian matrix of G.
        The graph Laplacian is the matrix L = |D| - A, where
        A is the adjacency matrix and |D| is the diagonal matrix of absolute 
        values of node degrees

        Returns
        -------
        L : SciPy sparse array
        The Laplacian matrix of G.
        """
        return self.sDeg - self.Adj
    #
    def compute_laplacian_spectrum(self, MODE_lapspec: str = 'numpy') -> None:
        if MODE_lapspec == 'networkx':
            self.slspectrum = nx.laplacian_spectrum(self.system.G)
        elif MODE_lapspec == 'numpy':
            self.slspectrum = eigvals(self.sLp.toarray())
    #
    def timescale_for_S(self) -> ndarray:
        return np.logspace(self.t1, self.t2, self.steps)
    #
    def timescale_for_C(self) -> ndarray:
        t = self.timescale_for_S()
        return .5 * (t[1:] + t[:-1])
    #
    def compute_entropy(self) -> None:# tuple[ndarray, ndarray, ndarray, ndarray]
        if self.slspectrum is None:
            self.compute_laplacian_spectrum()
        t = self.timescale_for_S()
        w = self.slspectrum
        S = np.zeros(len(t))
        VarL = np.zeros(len(t))
    
        for i, tau in enumerate(t):
            rhoTr = np.exp(-tau * w)
            Tr = np.nansum(rhoTr)
            rho = np.divide(rhoTr, Tr)
            avgrho = np.nansum(np.multiply(w, rhoTr)) / Tr
            av2rho = np.nansum(np.multiply(np.multiply(w,w), rhoTr)) / Tr
            S[i] = -np.nansum(rho * np.log(rho)) / np.log(self.system.N)
            VarL[i] = (av2rho - avgrho**2)
        self.Sm1  = 1 - S
        self.VarL = VarL
    #
    def compute_Cspe(self) -> None:
        if self.Sm1 is None:
            self.compute_entropy()
        N = self.system.N
        Sm1 = self.Sm1
        t = self.timescale_for_S()
        self.Cspe = np.log(N) * dv(Sm1, np.log(t))
    #
    def compute_taumax_array(self) -> None:
        if self.Cspe is None:
            self.compute_Cspe()
        t = self.timescale_for_C()
        maxIdx = argrelextrema(self.Cspe, np.greater)[0]
        maxIdxCondition = self.Cspe[maxIdx] > self.maxThresh
        self.taumax = t[maxIdx[maxIdxCondition]]
        self.taumax0 = t[maxIdx[maxIdxCondition][0]]
    #
    def laplacian_dynamics_init(self):
        if self.initCond == 'uniform_1':
            self.status_array = np.random.uniform(-1, 1, self.system.N)
        elif self.initCond == 'delta_1':
            self.status_array = np.zeros(self.system.N)
            self.status_array[self.system.N//2] = self.system.N
        elif self.initCond == 'gauss_1':
            self.status_array = np.random.normal(-1, 1, self.system.N)
        elif self.initCond == 'all_1':
            self.status_array = np.ones(self.system.N)
        elif self.initCond.startswith('ground_state'):
            self.eigenModeInit = int(self.initCond.split('_')[-1])
        #
        if self.system.pbc is False:
            L = int(np.sqrt(self.system.N))
            self.fixed_border_idxs = np.array(sorted([i for i in range(L)] + \
                        [(L-1) * L + i for i in range(L)] + \
                        [i * L for i in range(1, L-1)] + 
                        [(i+1) * L - 1 for i in range(1, L-1)]))
            self.status_array[self.fixed_border_idxs] = self.fbc_val
#
class Lattice2D(Graph):
    p_c = None
    lsp = None
    def __init__(self, side1: int, geometry: str = 'squared', side2: int = 0,
                 lsp_mode: str = 'intervals', incoming_graph_data=None, pbc = True, 
                 **attr) -> None:
        super().__init__(incoming_graph_data, **attr)
        self.side1 = side1
        self.side2 = side2 if side2 else side1
        self.geometry = geometry
        self.pbc = pbc
        self.G = self.lattice_selection()
        self.H = nx.convert_node_labels_to_integers(self.G)
        self.node_map = dict(zip(self.G, self.H))
        self.edge_map = dict(zip(self.G.edges(), self.H.edges()))
        self.N = self.G.number_of_nodes()
        self.eset = self.G.edges()
        self.Ne = self.G.number_of_edges()
        self.lsp_mode = lsp_mode

    def lattice_selection(self) -> Graph:
        if self.geometry == 'triangular':
            nxfunc = nx.triangular_lattice_graph
            self.p_c = 0.146
        elif self.geometry == 'squared':
            nxfunc = nx.grid_2d_graph
            self.p_c = 0.103
        elif self.geometry == 'hexagonal':
            nxfunc = nx.hexagonal_lattice_graph
            self.p_c = 0.065
        return nxfunc(self.side1, self.side2, periodic=self.pbc)
    
    def lsp_selection(self, custom_list):
        if self.lsp_mode == 'custom':
                self.lsp = np.array(custom_list)
        elif self.lsp_mode == 'intervals':
            intervals = []
            tmp = max([vset['rsf'] for vset in custom_list])
            for vset in custom_list:
                if vset['kind'] == 'log':
                        spacing_f = np.logspace
                        vset['start'] = np.log10(vset['start'])
                        vset['stop'] = np.log10(vset['stop'])
                elif vset['kind'] == 'lin':
                        spacing_f = np.linspace
                intervals.append(#
                    round_sigfig_n(#
                        spacing_f(vset['start'], vset['stop'],
                                    num=vset['num'],
                                    endpoint=False),
                    vset['rsf'])
                )
            self.lsp = (intervals := np.concatenate(intervals))
            while set(self.lsp).__len__() == intervals.__len__():
                tmp = tmp - 1
                self.lsp = np.round(self.lsp , tmp)
            tmp = tmp + 1
            self.lsp = np.round(intervals, tmp)
    
    def default_dict_lsp(self, num_low = 3, num_at = 6, num_high = 3):
        d = (#
            {'kind': 'lin', 'start': 0.001, 'stop': self.p_c-self.p_c*num_at/100, 
              'num': num_low, 'rsf': 1}, 
             {'kind': 'lin', 'start': self.p_c-self.p_c*num_at/100, 
              'stop': self.p_c+self.p_c*num_at/100, 'num': num_at, 'rsf': 3},
              {'kind': 'lin', 'start': self.p_c+self.p_c*num_at/100, 
              'stop': 1, 'num': num_high, 'rsf': 1},
        )
        return d
    


def lsp_read_values(folder_path):
    #file_pattern = r"p=(\d+\.\d+)_Sm1.bin"
    file_pattern = r"Sm1_avg_p=(\d+\.\d+)"
    value_pattern = r"p=(\d+\.\d+)"
    #
    # Get all files in the folder
    #
    files = os.listdir(folder_path)
    #
    # Filter files based on the pattern
    file_names = [file_name for file_name in files if re.match(file_pattern, file_name)]
    #
    # Extract values from file names
    values = []
    for file_name in file_names:
        match = re.search(value_pattern, file_name)
        if match:
            value = float(match.group(1))
            values.append(value)
    #
    # Sort the values if needed
    values.sort()
    return np.array(values) 
def round_sigfig_n(num, n: int = 1):
    if n not in range(1, MAX_DIGITS_ROUND_SIGFIG):
        raise ValueError("Significant figures number not in [1, 15].")
    expn = -np.floor(np.log10(np.abs(num))).astype('int')
    if hasattr(num, "__len__"):
        rr = np.array([np.round(nn, ee+n-1) for nn,ee in zip(num, expn)])
    else:
        rr = np.round(num, expn+n-1)
    return rr




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
        lapl = slaplacian_matrix(G)
    else:
        lapl = nx.laplacian_matrix(G)
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


def Cspe_plot_ax(ax):
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\log(N)\langle{C}\rangle$")
    ax.set_xscale('log')

def log_binning(data, binnum=20):
    log_data = np.log10(data)
    min_val = int(np.floor(np.min(log_data)))
    max_val = int(np.ceil(np.max(log_data)))
    bins = np.logspace(min_val, max_val, num=binnum)
    hist, _ = np.histogram(data, bins=bins)
    bin_w = (bins[1:] + bins[:-1])
    bin_centers = (bins[1:] + bins[:-1]) / 2.0
    return bin_centers, hist, bin_w

def neglog_binning(data, binnum=20):
    abs_data = np.abs(data)
    log_data = np.log10(abs_data)
    min_val = int(np.floor(np.min(log_data)))
    max_val = int(np.ceil(np.max(log_data)))
    bins = -np.logspace(min_val, max_val, num=binnum)[::-1]  # Reverse the bins to have a negative x-scale
    hist, _ = np.histogram(data, bins=bins)
    bin_w = (bins[1:] + bins[:-1])
    bin_centers = (bins[1:] + bins[:-1]) / 2.0
    return bin_centers, hist, bin_w


def symlog_binning(full_data, binnum=20):
    datap = full_data[full_data > 0]
    datam = full_data[full_data < 0]
    if datap.size:
        outp = log_binning(datap)
    else:
        outp = None
    if datam.size:
        outm = neglog_binning(datam)
    else:
        outm = None
    return outp, outm

def create_custom_colormap(c1="#0000ff", c2="#fc0303"):
    colors = [c1, c2]  # Red to black
    cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)
    return cmap

def radial_correlation(field, neighbor_dict, central_index):
    # Flatten the field array
    flattened_field = field.flatten()
    # Compute the mean of the field values
    mean_field = np.mean(flattened_field)
    # Compute the correlation function
    correlation = []
    for ord in range(1, int(np.sqrt(len(flattened_field)))):
        kth_order_neighbors = [n for n, d in neighbor_dict.items() if d == ord]
        for neighbor_index in kth_order_neighbors:
            correlations_node = []
            correlation_node = np.mean((flattened_field[central_index] - mean_field) * flattened_field[neighbor_index]) / np.std(flattened_field) ** 2
            correlations_node.append(correlation_node)
        correlation.append([ord, np.mean(correlations_node)])

    return correlation
# Define a function to calculate the Local Moran's I for a given cell
def local_moran_i(i, field_array, adj):
    N = len(field_array)
    x_avg = np.mean(field_array)
    x_i = field_array[i]
    x_j = field_array
    w_ij = adj[i]
    m2 = np.sum((x_i-x_avg)**2)/N
    I_i = (x_i-x_avg)/m2 * np.sum(w_ij*(x_j-x_avg))
    # # Extract the values of the cell and its neighbors
    # cell_value = field_array[row, col]
    # cell_neighbors = field_array[max(0, row - window_size):min(n_rows, row + window_size + 1),
    #                  max(0, col - window_size):min(n_cols, col + window_size + 1)].flatten()
    # # Calculate the Moran's I numerator (sum of (x_i - x_mean) * (x_j - x_mean) for all i and j)
    # numerator = (cell_value - field_mean) * np.sum(cell_neighbors - field_mean)
    # # Calculate the Moran's I denominator (sum of (x_i - x_mean)^2 for all i)
    # denominator = np.sum((cell_value - field_mean) ** 2)/(n_rows * n_cols)
    # # Calculate Moran's I statistic
    # moran_i = (window_size + 1) * numerator / denominator
    return I_i


