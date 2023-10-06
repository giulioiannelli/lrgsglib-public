#
import argparse
import random
import re
import os
import sys
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
#
import networkx as nx
import numpy as np
import random as rd

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
from matplotlib.patches import Circle, Rectangle
from matplotlib.ticker import ScalarFormatter
from networkx.classes.graph import Graph
from numba import jit
from numpy import inf, ndarray
from numpy.linalg import eigvals, eigvalsh
from numpy.typing import NDArray
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
from subprocess import call
from typing import Any, Optional, Union

#
from tqdm import tqdm

#
from .nx_patches import *
from .nx_objects import *
from .LRGSG_const import *
from .LRGSG_errwar import *
from .LRGSG_plots import *
from .LRGSG_utils import *

np.random.seed(0)
random.seed(0)

#
__all__ = [
    "np",
    "nx",
    "plt",
    "tqdm",
    "flip_random_fract_edges",
    "SignedLaplacianAnalysis",
    "lsp_read_values",
    "FullyConnected",
    "Lattice2D",
]


#
# renormalization group for heterogenous network functions
class SignedLaplacianAnalysis:
    slspectrum = None
    Sm1 = None
    VarL = None
    Cspe = None
    taumax = None
    taumax0 = None
    frames_dynsys = []
    Ising_clusters = []
    eigv = None
    ACCERR_LAPL_DYN = 1e-10
    MAXVAL_LAPL_DYN = 200

    #
    def __init__(
        self,
        system,
        is_signed: bool = True,
        steps: int = DEFAULT_ENTROPY_STEPS,
        t1: float = -2,
        t2: float = 5,
        maxThresh: float = DEFAULT_MAX_THRESHOLD,
        nreplica: int = 0,
        initCond: str = "gauss_1",
        t_steps=10,
        no_obs=1,
    ) -> None:
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
        self.t_steps = t_steps
        self.no_obs = no_obs


    #
    def timescale_for_S(self) -> ndarray:
        return np.logspace(self.t1, self.t2, self.steps)
    #
    def timescale_for_C(self) -> ndarray:
        return 0.5 * (self.tTsS[1:] + self.tTsS[:-1])
    #
    def init_computation(self):
        if self.system.slspectrum is None:
            self.system.compute_laplacian_spectrum()
        self.tTsS = self.timescale_for_S()
        self.tTsC = self.timescale_for_C()
    #
    def compute_entropy(self) -> None:
        w = self.system.slspectrum
        S = np.zeros(len(self.tTsS))
        VarL = np.zeros(len(self.tTsS))
        for i, tau in enumerate(self.tTsS):
            rhoTr = np.exp(-tau * w)
            Tr = np.nansum(rhoTr)
            rho = np.divide(rhoTr, Tr)
            avgrho = np.nansum(np.multiply(w, rhoTr)) / Tr
            av2rho = np.nansum(np.multiply(np.multiply(w, w), rhoTr)) / Tr
            S[i] = -np.nansum(rho * np.log(rho)) / np.log(self.system.N)
            VarL[i] = av2rho - avgrho**2
        self.Sm1 = 1 - S
        self.VarL = VarL
    #
    def compute_Cspe(self) -> None:
        if self.Sm1 is None:
            self.compute_entropy()
        self.Cspe = np.log(self.system.N) * dv(self.Sm1, np.log(self.tTsS))
    #
    def compute_taumax_array(self) -> None:
        if self.Cspe is None:
            self.compute_Cspe()
        maxIdx = argrelextrema(self.Cspe, np.greater)[0]
        maxIdxCondition = self.Cspe[maxIdx] > self.maxThresh
        self.taumax = self.tTsC[maxIdx[maxIdxCondition]]
        self.taumax0 = self.tTsC[maxIdx[maxIdxCondition][0]]
    #
    def laplacian_dynamics_init(
        self,
        t_stepsMultiplier: int = 1,
        window_size=0,
        window_shift_x=0,
        window_shift_y=0,
        win_val=1,
    ):
        N = self.system.N
        #
        self.Deltat = 1.0 / self.t_steps
        self.simulationTime = N * self.t_steps * t_stepsMultiplier
        self.sampling = np.logspace(
            0, np.log10(self.simulationTime), num=self.no_obs, dtype=int
        )
        #
        if not "ground_state" in self.initCond:
            self.compute_k_eigvV()
        if self.initCond.startswith("ground_state"):
            self.eigenModeInit = int(self.initCond.split("_")[-1])
            self.compute_k_eigvV(howmany=self.eigenModeInit + 1)
            self.field = self.system.eigV.T[self.eigenModeInit]
        elif self.initCond == "uniform_1":
            self.field = np.random.uniform(-1, 1, N)
        elif self.initCond == "delta_1":
            self.field = np.zeros(N)
            self.field[N // 2] = 1
        elif self.initCond == "gauss_1":
            self.field = np.random.normal(-1, 1, N)
        elif self.initCond.startswith("all"):
            self.initVal = float(self.initCond.split("_")[-1])
            self.field = self.initVal * np.ones(N)
        elif self.initCond.startswith("window"):
            s22 = self.system.side2 // 2
            wndwS = s22 - 1 if window_size > (s22 - 1) else window_size
            hS, hE = s22 - wndwS - 1, s22 + wndwS + 1
            initStatus = np.zeros(N)
            #
            if self.initCond.startswith("window_multiple"):
                self.nsquares = int(self.initCond.split("_")[-1])
                wndwSa = np.array([wndwS, -wndwS])
                sqTmp = np.random.randint(
                    wndwS, self.system.side1 - wndwS, size=(self.nsquares, 2)
                )
                result = np.column_stack((sqTmp, sqTmp + wndwS))
                result[:, [1, 2]] = result[:, [2, 1]]
                result = result.reshape(-1, 2, 2)
                sqIdx = np.concatenate(
                    [
                        np.concatenate(
                            [
                                [
                                    j + i * self.system.side1
                                    for j in range(*iSq[0])
                                ]
                                for i in range(*iSq[1])
                            ]
                        )
                        for iSq in result
                    ]
                )
                if type(win_val) is str:
                    if win_val == "gauss":
                        win_val = np.random.normal(0, 1, size=len(sqIdx))
                    elif win_val == "uniform":
                        win_val = np.random.uniform(-1, 1, len(sqIdx))
                initStatus[sqIdx] = np.ones(len(sqIdx)) * win_val
            else:
                shiftsX = [hS + window_shift_x, hE + window_shift_x]
                shiftsY = [hS + window_shift_y, hE + window_shift_y]
                sqIdx = np.concatenate(
                    [
                        [j + i * self.system.side1 for j in range(*shiftsX)]
                        for i in range(*shiftsY)
                    ]
                )
                if self.initCond == "window_1":
                    initStatus[sqIdx] = np.ones(len(sqIdx))
                elif self.initCond == "window_val":
                    initStatus[sqIdx] = np.ones(len(sqIdx)) * win_val
                elif self.initCond == "window_gauss_1":
                    initStatus[sqIdx] = np.random.normal(-1, 1, len(sqIdx))
                elif self.initCond.startswith("window_ground_state"):
                    self.eigenModeInit = int(self.initCond.split("_")[-1])
                    self.compute_k_eigvV(howmany=self.eigenModeInit + 1)
                    initStatus = self.eigV.T[self.eigenModeInit]
                    outIdx = np.setxor1d(np.arange(N), sqIdx)
                    initStatus[outIdx] = np.random.uniform(
                        initStatus.min(), initStatus.max(), len(outIdx)
                    )
                else:
                    print("Error, no mode for init laplacian dynamic chosen.")
            self.field = initStatus
        #
        if self.system.pbc is False:
            L = int(np.sqrt(N))
            self.fixed_border_idxs = np.array(
                sorted(
                    [i for i in range(L)]
                    + [(L - 1) * L + i for i in range(L)]
                    + [i * L for i in range(1, L - 1)]
                    + [(i + 1) * L - 1 for i in range(1, L - 1)]
                )
            )
            self.field[self.fixed_border_idxs] = self.system.fbc_val

    #
    def run_laplacian_dynamics(self, rescaled=False, saveFrames=False):
        self.frames_dynsys = []
        #
        x = self.field
        #
        def stop_conditions_lapdyn(self, x_tm1, xx):
            ERRTOL = self.ACCERR_LAPL_DYN * np.ones(self.system.N)
            C1 = (np.abs(x_tm1 / x_tm1.max() - xx / xx.max()) < ERRTOL).all()
            C2 = np.abs(np.log10(np.max(np.abs(xx)))) > self.MAXVAL_LAPL_DYN
            return C1, C2
        #
        if rescaled:
            if rescaled == 'dynamic':
                lap = lambda t: np.exp(-self.eigv[0] * t) * self.sLp
            else:
                self.system.rescaled_laplacian(rescaled)
                lap = lambda _: self.system.resLp
        else:
            lap = lambda _: self.sLp
        #
        if not self.system.pbc:
            def set_bc(self):
                x[self.fixed_border_idxs] = self.fbc_val
        else:
            def set_bc():
                pass
        #
        if saveFrames:
            def save_frames(self, x, t):
                if t in self.sampling:
                    self.frames_dynsys.append(x)
        else:
            def save_frames(*_):
                pass
        #
        print("Beginning Laplacian dynamics.")
        #
        for t in tqdm(range(1, self.simulationTime)):
            save_frames(self, x, t)
            x_tm1 = x
            x = x - self.Deltat * lap(t) @ x # + np.sqrt(Deltat)*noise 
            set_bc()
            C1, C2 = stop_conditions_lapdyn(self, x_tm1, x)
            if C2 or C1:
                if C1:
                    print("Convergence reached.")
                if C2:
                    print("Max val. reached.")
                break
        self.field = x
    #
    class IsingDynamics:
        magn = []
        ene = []
        magn_array_save = []
        Ising_clusters = []
        def __init__(self, system, T: float,  IsingIC: str,  nstepsIsing: int = 100, save_magnetization: bool = False,) -> None:
            self.system = system
            self.T = T
            self.IsingIC = IsingIC
            self.nstepsIsing = nstepsIsing
            self.save_magnetization=save_magnetization
            pass

        def boltzmann_factor(self, energy: float) -> float:
            return np.exp(-energy / self.T)
        #
        def neigh_weight_magn(self, node: int) -> list:
            node_dict = dict(self.system.H[node])
            return [w["weight"] * self.m[nn] for nn, w in node_dict.items()]
        #
        def neigh_ene(self, neigh: list) -> float:
            return np.sum(neigh) / len(neigh)
        #
        def flip_spin(self, node: int):
            self.m[node] = -self.m[node]

        def metropolis(self, node):
            neigh = self.neigh_weight_magn(node)
            neighene = self.neigh_ene(neigh)
            E_old = -self.m[node] * neighene
            E_new = +self.m[node] * neighene
            DeltaE = E_new - E_old
            if DeltaE < 0:
                self.flip_spin(node)
            elif np.random.uniform() < self.boltzmann_factor(DeltaE):
                self.flip_spin(node)

        #
        def calc_full_energy(self):
            return np.array(
                [
                    self.neigh_ene(self.neigh_weight_magn(node))
                    for node in range(self.system.N)
                ]
            ).sum()
        #
        def init_ising_dynamics(self):
            if self.IsingIC == "uniform":
                self.m = np.random.choice([-1, 1], size=self.system.N)
            elif self.IsingIC.startswith("ground_state"):
                self.m = np.array([1 if i in self.system.bin_eigV else -1 for i in range(self.system.N)])
                if len(self.IsingIC) > len("ground_state"):
                    IICnoise = float(self.IsingIC.split("_")[-1])
                    if IICnoise < 0 or IICnoise > 1:
                        print("error")
                    self.m = self.m*np.random.choice([-1, 1], size=self.system.N, p=[IICnoise, 1-IICnoise])
        #
        def run(self, MODE = 'py'):
            if MODE == 'C':
                self.id_string_isingdyn = randstring()
                output_file = open(f"src/LRGSG_package/tmp_stuff/m_{self.id_string_isingdyn}.bin", 'wb')
                self.m.tofile(output_file)
                call(["./spa", "args", "to", "spa"])
                # execute C program,  calling with proper arguments
                return
            metropolis_1step = np.vectorize(self.metropolis, excluded='self')
            if self.save_magnetization:
                def save_magn_array():
                    self.magn_array_save.append(self.m)
            else:
                def save_magn_array():
                    pass
            sample = rd.sample(list(self.system.H.nodes()), self.system.N)
            for _ in range(self.nstepsIsing):
                self.magn.append(np.sum(self.m))
                self.ene.append(self.calc_full_energy())
                # for i in range(self.system.N):
                #     self.metropolis(sample[i])
                metropolis_1step(sample)
                save_magn_array()
        #
        def find_ising_clusters(self):
            if self.Ising_clusters:
                print("exit function")
                return
            lnodes = list(self.system.H.nodes())
            lnodes_tmp = lnodes[:]
            #
            self.system.compute_k_eigvV()
            eigVbin = self.system.eigV.T[:, 0]
            eigVbin[eigVbin >= 0] = +1
            eigVbin[eigVbin < 0] = -1

            def recursive_search(seed, magn_i, clustertmp):
                neighs = get_kth_order_neighbours(self.system.H, seed, 1)
                neighs = np.array([e for e in neighs if e not in set(clustertmp)])
                if not neighs.size:
                    return
                samecluster = np.array(eigVbin[neighs] == magn_i)
                if not samecluster.any():
                    return
                neighs_samecluster = list(neighs[samecluster])
                clustertmp.extend(neighs_samecluster)
                for ss in neighs_samecluster:
                    recursive_search(ss, magn_i, clustertmp)

            #
            for i in lnodes:
                if i not in lnodes_tmp:
                    continue
                if not lnodes_tmp:
                    break
                #
                clustertmp = []
                clustertmp.extend([i])
                #
                recursive_search(i, eigVbin[i], clustertmp)
                lnodes_tmp = [e for e in lnodes_tmp if e not in set(clustertmp)]
                self.Ising_clusters.append(clustertmp)
            self.numIsing_cl = len(self.Ising_clusters)
            self.biggestIsing_clidx = np.argmax(
                [len(c) for c in self.Ising_clusters]
            )
            self.biggestIsing_cl = self.Ising_clusters[self.biggestIsing_clidx]
        #
        def mapping_nodes_to_clusters(self):
            if not self.Ising_clusters:
                self.find_ising_clusters()
            loc = [x for x in range(len(self.Ising_clusters))]
            self.loc = loc
            node_with_inherclust = [
                [[j, loc[i]] for j in clus]
                for i, clus in enumerate(self.Ising_clusters)
            ]
            self.node_with_inherclust = node_with_inherclust
            node_inherclust_flat = [i for j in node_with_inherclust for i in j]
            self.node_inherclust_flat = node_inherclust_flat
            sorted_list = sorted(node_inherclust_flat, key=lambda x: x[0])
            self.sorted_list = sorted_list
            result_array = np.empty(
                (self.system.side1, self.system.side2), dtype=object
            )
            self.result_array = result_array

            # Fill the result_array with tuples from sorted_list
            for i, sublist in enumerate(sorted_list):
                row, col = divmod(
                    i, self.system.side1
                )  # Calculate the row and column index
                result_array[row, col] = sublist[1]
            # self.distr = np.unique(result_array, return_counts=True)
            self.mapping = result_array

    #
    def rescaled_field_regularization(self):
        status = self.field.reshape(self.system.side1, self.system.side2)
        restatus = np.log10(np.max(status) - status)
        nnans = restatus[(restatus != np.inf) & (restatus != -np.inf)]
        self.restatus = np.nan_to_num(
            restatus, posinf=np.max(nnans), neginf=np.min(nnans)
        )

    #
    def make_animation_fromFrames(self, savename="output", fps=10, dpi=200):
        no_frames = len(self.frames_dynsys)
        print("# of frames: ", no_frames)
        #
        fig, ax = plt.subplots()
        #
        animate = self.system.make_animation(fig, ax, self.frames_dynsys)
        ani = animation.FuncAnimation(fig, animate, frames=no_frames)
        #
        fig.tight_layout()
        ani.save(
            f"{savename}{eMP4}", writer=animation.FFMpegWriter(fps=fps), dpi=dpi
        )
        plt.close(fig)        







#
def lsp_read_values(folder_path, fpattern="Sm1_avg_p", sort=True):
    file_pattern = rf"{fpattern}=(\d+\.\d+)"
    value_pattern = r"p=(\d+\.\d+)"
    # Get all files in the folder
    files = os.listdir(folder_path)
    # Filter files based on the pattern
    file_names = [
        file_name for file_name in files if re.match(file_pattern, file_name)
    ]
    # Extract values from file names
    values = []
    for file_name in file_names:
        match = re.search(value_pattern, file_name)
        if match:
            value = float(match.group(1))
            values.append(value)
    # Sort the values if needed
    if sort:
        values.sort()
    return np.array(values)


#                           ,((((((((((((((((((((((((.
#                      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                   @@@@@@                              @@@@@@
#                  @@@@@@                               *@@@@@/
#    /@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                #    TRASH CODE
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,
#    @@@@@@@@@@@@@@    @@@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    &@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@&
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    .@@@@@@@@@@@@@
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    (@@@@@@@@@@@@@
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    @@@@@@@@@@@@@@
#     @@@@@@@@@@@@@.    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    @@@@@@@@@@@@@@
#     @@@@@@@@@@@@@/    @@@@@@@@@@@@@@     @@@@@@@@@@@@@&    @@@@@@@@@@@@@@
#     *@@@@@@@@@@@@@    @@@@@@@@@@@@@@     @@@@@@@@@@@@@(    @@@@@@@@@@@@@.
#      @@@@@@@@@@@@@    &@@@@@@@@@@@@@     @@@@@@@@@@@@@.    @@@@@@@@@@@@@
#      @@@@@@@@@@@@@    *@@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      &@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@#
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#       *@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@@    @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&
#          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&



    # def run_ising_dynamics(
    #     self,
    #     T=0.1,
    #     nstepsIsing=100,
    #     IsingIC="uniform",
    #     save_magnetization=False,
    #     tqdm_on=False,
    #     MODE='C'
    # ):

    #     self.init_ising_dynamics(IsingIC)
    #     self.id_string_isingdyn = randstring()
    #     if MODE == 'C':
    #         output_file = open(f"src/LRGSG_package/tmp_stuff/m_{self.id_string_isingdyn}.bin", 'wb')
    #         self.m.tofile(output_file)
    #         return
        
    #     magn = []
    #     ene = []
    #     if save_magnetization:
    #         self.magn_array_save = []
    #     ising = self.IsingDynamics(self.system.N)
    #     #
    #     def boltzmann_factor(energy, temp):
    #         return np.exp(-energy / temp)
    #     #questo ci puo calcolare una sola volta
    #     def neigh_weight_magn(m: NDArray, node_dict: dict) -> list:
    #         return [w["weight"] * m[nn] for nn, w in node_dict.items()]

    #     def neigh_ene(neigh: list) -> float:
    #         return - np.sum(neigh) / neigh.__len__()

    #     def calc_full_energy(m: NDArray, graph):
    #         """Energy of a given configuration"""
    #         return np.array(
    #             [
    #                 m[i] * neigh_ene(neigh_weight_magn(m, dict(graph.H[i])))
    #                 for i in range(graph.N)
    #             ]
    #         ).sum()

    #     # import time

    #     # neigh_time1 = 0
    #     # neigh_time2 = 0
    #     # fliptime = 0

    #     def flip_spin(node, m, graph, T):
    #         # nonlocal neigh_time1
    #         # nonlocal neigh_time2
    #         # nonlocal fliptime
    #         # nt1 = time.time()
    #         m_flp = -m[node]
    #         neigh = neigh_weight_magn(m, dict(graph.H[node]))
    #         # ent1 = time.time()
    #         # nt2 = time.time()
    #         neighEn = neigh_ene(neigh)
    #         E_old = m[node] * neighEn
    #         E_new = m_flp * neighEn
    #         # ent2 = time.time()
    #         # et1 = time.time()
    #         DeltaE = E_new - E_old
    #         if DeltaE < 0:
    #             m[node] = m_flp
    #         elif np.random.uniform() < boltzmann_factor(DeltaE, T):
    #             m[node] = m_flp
    #         # net1 = time.time()
    #         # neigh_time1 += ent1 - nt1
    #         # neigh_time2 += ent2 - nt2
    #         # fliptime += net1 - et1

    #     flip_spin_all = np.vectorize(flip_spin)
    #     flip_spin_all.excluded.add(1)
    #     flip_spin_all.excluded.add(2)
    #     flip_spin_all.excluded.add(3)
    #     #
    #     # s1 = time.time()
    #     # # a = np.random.randint(0, self.system.N, nstepsIsing*self.system.N)
    #     # t2_tot = 0.0
    #     # t3_tot = 0.0

    #     iterator = tqdm(range(nstepsIsing)) if tqdm_on else range(nstepsIsing)
        
    #     for _ in range(nstepsIsing):
    #         magn.append(np.sum(m))
    #         ene.append(calc_full_energy(m, self.system))
    #         # s2 = time.time()
    #         sample = rd.sample(self.system.H.nodes(), self.system.N)
    #         # e2 = time.time()
    #         # t2_tot += e2 - s2
    #         # s3 = time.time()
    #         # for i in range(self.system.N):
    #         #     # node =   # a[t*nstepsIsing + i]#sample[i]
    #         #     flip_spin(sample[i], m, self.system, T)
    #         flip_spin_all(sample, m, self.system, T)
    #         # e3 = time.time()
    #         # t3_tot += e3 - s3
    #         if save_magnetization:
    #             self.magn_array_save.append(np.array(m))
    #     self.magn_array = m
    #     # e1 = time.time()
    #     # print("time for full:", e1 - s1)
    #     # print("time for sampling:", t2_tot)
    #     # print("time for flipping:", t3_tot)
    #     # print("time for flipping_neig1:", neigh_time1)
    #     # print("time for flipping_neig2:", neigh_time2)
    #     # print("time for flipping_flip:", fliptime)
    #     return magn, ene








#     def run_ising_dynamics(self):
#         import random as rd

# N = SLRG_obj.system.N
# nmax= 1000
# T = 0.001
# eigV = SLRG_obj.eigV
# H = SLRG_obj.system.H

# def calcEnergy(m, H):
#     '''Energy of a given configuration'''
#     energy = 0
#     for i in range(N):
#         neigh=[w['weight']*m[nn] for nn, w in dict(H[i]).items()]
#         energy += -m[i]*np.sum(neigh)
#     return energy/4.

# #Initialize magnetization
# labels=np.where((eigV[0] < 0))[0]
# lista=['red' if i in labels else 'blue' for i in range(H.number_of_nodes())]
# bigene = []
# for replica in range(1):
#     if replica != 0:
#         m=[1 if i in labels else -1 for i in range(H.number_of_nodes())]
#     else:
#         m = np.random.choice([-1, 1], size=N)

#     #Metropolis
#     magn = []
#     ene = []
#     sample = rd.sample(H.nodes, N)
#     for nsteps in tqdm(range(100)):
#         magn.append(np.sum(m))
#         ene.append(calcEnergy(m, H))
#         for i in range(N):
#             node=sample[i]
#             m_new = -m[node]

#             #Metropolis thing
#             neigh=[w['weight']*m[nn] for nn, w in dict(H[node]).items()]
#             E_old=-m[node]*np.sum(neigh)
#             E_new=-m_new*np.sum(neigh)

#             if E_new<E_old:
#                 m[node]=m_new
#             else:
#                 r=rd.uniform(0, 1)
#                 if (r<np.exp(-(E_new-E_old)/T)):
#                     m[node]=m_new
#     bigene.append(ene)


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
        L = D - A
        w = eigvals(L)
    else:
        L = nx.laplacian_matrix(G).todense()
        w = nx.laplacian_spectrum(G)
    return L, w


def get_graph_lspectrum_rw(G, is_signed=False):
    A = nx.adjacency_matrix(G).toarray()
    D = np.diag(np.abs(A).sum(axis=1))
    L = np.eye(D.shape[0]) - fractional_matrix_power(
        D, -0.5
    ) @ A @ fractional_matrix_power(D, -0.5)
    if is_signed:
        w = eigvals(L)
    else:
        w = nx.laplacian_spectrum(G)
    return L, w


def entropy(G, steps=600, is_signed=False, wTresh=1e-15, t1=-2, t2=5):
    N = G.number_of_nodes()

    L, w = get_graph_lspectrum_rw(G, is_signed=is_signed)
    wSig = w[w > wTresh]

    t = np.logspace(t1, t2, int(steps))
    S = np.zeros(len(t))
    VarL = np.zeros(len(t))

    for i, tau in enumerate(t):
        rhoTr = np.exp(-tau * w)
        Tr = np.nansum(rhoTr)
        rho = np.divide(rhoTr, Tr)
        S[i] = -np.nansum(rho * np.log(rho)) / np.log(N)
        avgrho = np.nansum(np.multiply(w, rhoTr)) / Tr
        av2rho = np.nansum(np.multiply(np.multiply(w, w), rhoTr)) / Tr
        VarL[i] = av2rho - avgrho**2

    dS = np.log(N) * np.diff(1 - S) / np.diff(np.log(t))
    return 1 - S, dS, VarL, t


def averaged_Sm1(t1Sm1Avg):
    lenLs = t1Sm1Avg[0][0].__len__()
    avgSm1 = np.zeros(lenLs)
    all_t1 = np.concatenate([i[0] for i in t1Sm1Avg])
    commonLs = np.logspace(np.log10(min(all_t1)), np.log10(max(all_t1)), lenLs)
    digitizedLs = []
    for tS in t1Sm1Avg:
        digiTmp = np.digitize(tS[0], bins=commonLs) - 1
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
    num = expm((-tau * L))
    rho = num / np.trace(num)
    Trho = np.copy(1.0 / rho)
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
    linkage_matrix = linkage(dists / tmax, method=method)
    return linkage_matrix, w


def Cspe_plot_ax(ax):
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\log(N)\langle{C}\rangle$")
    ax.set_xscale("log")


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
            correlation_node = (
                np.mean(
                    (flattened_field[central_index] - mean_field)
                    * flattened_field[neighbor_index]
                )
                / np.std(flattened_field) ** 2
            )
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
    m2 = np.sum((x_i - x_avg) ** 2) / N
    I_i = (x_i - x_avg) / m2 * np.sum(w_ij * (x_j - x_avg))
    # # Extract the values of the cell and its neighbors
    # cell_value = field_array[row, col]
    # cell_neighbors = field_array[max(0, row - wndwS):min(n_rows, row + wndwS + 1),
    #                  max(0, col - wndwS):min(n_cols, col + wndwS + 1)].flatten()
    # # Calculate the Moran's I numerator (sum of (x_i - x_mean) * (x_j - x_mean) for all i and j)
    # numerator = (cell_value - field_mean) * np.sum(cell_neighbors - field_mean)
    # # Calculate the Moran's I denominator (sum of (x_i - x_mean)^2 for all i)
    # denominator = np.sum((cell_value - field_mean) ** 2)/(n_rows * n_cols)
    # # Calculate Moran's I statistic
    # moran_i = (wndwS + 1) * numerator / denominator
    return I_i


def ising_spinglass_pmJ_2D_Tcrit(L):
    return L ** (-1.0 / 2)


def make_animation_fromFrames(frames, savename="output.mp4", fps=10, dpi=200):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # I like to position my colorbars this way, but you don't have to
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", "5%", "5%")
    cv0 = frames[0]
    im = ax.imshow(cv0)  # Here make an AxesImage rather than contour
    cb = fig.colorbar(im, cax=cax)

    # tx = ax.set_title('Frame 0')
    def animate(i):
        arr = frames[i]
        vmax = np.max(arr)
        vmin = np.min(arr)
        im.set_data(arr)
        im.set_clim(vmin, vmax)
        # tx.set_text('Frame {0}'.format(i))
        # In this version you don't have to do anything to the colorbar,
        # it updates itself when the mappable it watches (im) changes

    print("# of frames: ", len(frames))
    ani = animation.FuncAnimation(fig, animate, frames=len(frames))
    fig.tight_layout()
    writervideo = animation.FFMpegWriter(fps=fps)
    ani.save(savename, writer=writervideo, dpi=dpi)


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
    noflip = int(p * nedges)
    #
    if noflip < 1:
        print("p too small")
        exit()
    rndsmpl = random.sample(range(nedges), noflip)
    #
    # all_weights = {e: 1 for e in eset}
    neg_weights = {e: -1 for i, e in enumerate(eset) if i in rndsmpl}
    #
    nx.set_edge_attributes(G, values=1, name="weight")
    nx.set_edge_attributes(G, values=neg_weights, name="weight")


#
def flip_one_2dgraph(G, coord1, coord2):
    neg_weights = {(coord1, coord2): -1}
    #
    nx.set_edge_attributes(G, values=neg_weights, name="weight")


#
lambdaPath_l2d = lambda geometry: f"l2d_{geometry}/"
pltPath_l2d = lambda geometry: f"data/plot/{lambdaPath_l2d(geometry)}"
datPath_l2d = lambda geometry: f"data/{lambdaPath_l2d(geometry)}"
setPath_ERp = "conf/ERp/"
pltPath_Sm1C = "plot/Sm1_and_C/"
#
pflip_fmt = ".3g"
