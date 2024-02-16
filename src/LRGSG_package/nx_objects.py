import os
import pickle
import random
import scipy

import networkx as nx
import numpy as np
import scipy.sparse as scsp

from .config.LRGSG_const import *
from .config.LRGSG_errwar import *
from .config.LRGSG_plots import imshow_colorbar_caxdivider
from .config.LRGSG_utils import round_sigfig_n, sum_tuples
from matplotlib.colors import Colormap
from networkx.classes.graph import Graph
from .nx_patches import signed_spectral_layout, get_neighbors_within_distance
from scipy.sparse import csr_array
from typing import Union


class SignedGraph:
    p_c = None
    lsp = None
    slspectrum = None
    pflmin = DEFAULT_MIN_PFLIPVAL
    pflmax = DEFAULT_MAX_PFLIPVAL

    def __init__(
        self,
        G: Graph,
        import_on: bool = False,
        pflip: float = 0.0,
        lsp_mode: str = "intervals",
        expathc: str = "",
    ):
        self.__init_paths__()
        if not self.pflmin <= pflip <= self.pflmax:
            raise ValueError(f"pflip must be between {self.pflmin} and {self.pflmax}, inclusive. Received: {pflip}")
        else:
            self.pflip = pflip
        self.lsp_mode = lsp_mode
        self.import_on = import_on
        self.expath = (
            expathc if expathc else f"{self.DEFAULT_GRAPHDIR}{self.syshapePTH}"
        )
        self.isingpath = f"{self.DEFAULT_ISINGDIR}{self.syshapePTH}"
        self.voterpath = f"{self.DEFAULT_VOTERDIR}{self.syshapePTH}"
        self.lrgsgpath = f"{self.DEFAULT_LRGSGDIR}{self.syshapePTH}"
        self.phtrapath = f"{self.DEFAULT_PHTRADIR}{self.syshapePTH}"
        self.__make_dirs__()
        self.stdFname = self.stdFname + f"_p={self.pflip:.3g}"
        if import_on:
            self.graphfname = self.expath + self.stdFname
            self.G = self.__init_graph_fromfile__()
        else:
            self.G = G
            self.graphfname = self.expath + self.stdFname
        self.init_sgraph()

    #
    def __init_graph_fromfile__(self):
        return pickle.load(open(f"{self.graphfname}.pickle", "rb"))

    #
    def __init_paths__(self):
        self.datPath = f"{DEFAULT_DATA_OUTDIR}{self.sgpath}"
        self.pltPath = (
            f"{DEFAULT_DATA_OUTDIR}{DEFAULT_PLOT_OUTDIR}{self.sgpath}"
        )
        self.DEFAULT_GRAPHDIR = self.datPath + DEFAULT_GRAPH_OUTDIR
        self.DEFAULT_ISINGDIR = self.datPath + DEFAULT_ISING_OUTDIR
        self.DEFAULT_VOTERDIR = self.datPath + DEFAULT_VOTER_OUTDIR
        self.DEFAULT_LRGSGDIR = self.datPath + DEFAULT_LRGSG_OUTDIR
        self.DEFAULT_PHTRADIR = self.datPath + DEFAULT_PHTRA_OUTDIR

    #
    def __make_dirs__(self):
        for _ in [self.expath, self.isingpath, self.voterpath, self.lrgsgpath, 
                  self.phtrapath]:
            os.makedirs(_, exist_ok=True)

    #
    def init_weights(self):
        nx.set_edge_attributes(self.G, values=1, name="weight")

    #
    def number_of_negative_links(self):
        self.Ne_n = (
            np.array(list(nx.get_edge_attributes(self.H, "weight").values()))
            < 0
        ).sum()

    #
    def upd_G_graph(self):
        self.invnode_map = {v: k for k, v in self.node_map.items()}
        self.invedge_map = {v: k for k, v in self.edge_map.items()}
        self.G = nx.relabel_nodes(self.H, self.invnode_map)
        self.esetG = list(self.G.edges())
        self.number_of_negative_links()

    #
    def upd_H_graph(self):
        self.H = nx.convert_node_labels_to_integers(self.G)
        self.esetH = list(self.H.edges())
        self.node_map = dict(zip(self.G, self.H))
        self.edge_map = dict(zip(self.G.edges(), self.H.edges()))
        self.number_of_negative_links()

    def init_n_nodes_edges(self):
        self.N = self.G.number_of_nodes()
        self.Ne = self.G.number_of_edges()

    def init_sgraph(self):
        self.init_n_nodes_edges()
        if self.import_on:
            self.upd_H_graph()
            self.nflip = self.Ne_n
            self.pflip = self.nflip / self.Ne
            self.upd_graph_matrices()
            self.randsample = np.where(
                np.array(self.Adj.todense()).flatten() < 0
            )
        else:
            self.nflip = int(self.pflip * self.Ne)
            self.randsample = random.sample(range(self.Ne), self.nflip)
            self.init_weights()
            self.upd_H_graph()
        self.upd_G_graph()

    #
    def adjacency_matrix(self, weight: str = "weight"):
        return nx.to_scipy_sparse_array(self.H, weight=weight, format="csr")

    #
    def degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(scsp.spdiags(A.sum(axis=1), 0, *A.shape, format="csr"))

    #
    def absolute_degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(
            scsp.spdiags(abs(A).sum(axis=1), 0, *A.shape, format="csr")
        )

    #
    def laplacian_matrix(self) -> csr_array:
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
    def signed_laplacian(self) -> csr_array:
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
    def upd_graph_matrices(self, on_graph="H"):
        if on_graph == "G":
            motherNx = self.G
        elif on_graph == "H":
            motherNx = self.H
        self.Adj = self.adjacency_matrix()
        self.Deg = self.degree_matrix(self.Adj)
        self.sDeg = self.absolute_degree_matrix(self.Adj)
        self.Lap = self.laplacian_matrix()
        self.sLp = self.signed_laplacian()

    #
    def flip_sel_edges(self, neg_weights_dict=None, on_graph="H"):
        """Flips a specific edges of a graph G."""
        #
        if on_graph == "G":
            if neg_weights_dict is None:
                neg_weights_dict = self.neg_weights_dict.DEFAULT_NEG_WEIGHTS_DICT_G
            nx.set_edge_attributes(
                self.G, values=neg_weights_dict, name="weight"
            )
        elif on_graph == "H":
            if neg_weights_dict is None:
                neg_weights_dict = self.neg_weights_dict.DEFAULT_NEG_WEIGHTS_DICT_H
            nx.set_edge_attributes(
                self.H, values=neg_weights_dict, name="weight"
            )
        self.upd_graph(on_graph=on_graph)
        self.upd_graph_matrices()

    #
    def check_pflip(self):
        if self.nflip < 1:
            raise NflipError(
                """The probability of flipping an edge times the 
                             number of edges is < 1, then no edges would be
                             flipped. Skipping the analysis for this value."""
            )
    #
    def upd_graph(self, on_graph: str = "G"):
        if on_graph == "G":
            self.upd_H_graph()
        elif on_graph == "H":
            self.upd_G_graph()
    #
    def unflip_all(self, on_graph="H"):
        if on_graph == "G":
            eset = self.esetG
        elif on_graph == "H":
            eset = self.esetH
        self.WEIGHTS_DICT_H_PFLIP = {e: 1 for i, e in enumerate(eset)}
        self.flip_sel_edges(neg_weights_dict=self.WEIGHTS_DICT_H_PFLIP, on_graph=on_graph)

    #
    def flip_random_fract_edges(self, on_graph="H"):
        """Flips a fraction p of edges (+1 to -1) of a graph G."""

        #
        try:
            self.check_pflip()
        except NflipError:
            return None
        if on_graph == "G":
            eset = self.esetG
        elif on_graph == "H":
            eset = self.esetH
        self.NEG_WEIGHTS_DICT_H_PFLIP = {
            e: -1 for i, e in enumerate(eset) if i in self.randsample
        }
        self.flip_sel_edges(neg_weights_dict=self.NEG_WEIGHTS_DICT_H_PFLIP, on_graph=on_graph)

    #
    def compute_laplacian_spectrum(self, MODE_lapspec: str = "numpy") -> None:
        if MODE_lapspec == "networkx":
            self.slspectrum = nx.laplacian_spectrum(self.system.G)
        elif MODE_lapspec == "numpy":
            self.slspectrum = np.linalg.eigvalsh(self.sLp.toarray())
    #
    def compute_k_eigvV(self, MODE_dynspec: str = "scipy", howmany: int = 1, which: str = "SM"):
        if MODE_dynspec == "numpy" or howmany == self.N:
            self.eigv, self.eigV = np.linalg.eigh(self.sLp.astype(np.float64).todense())
            self.eigV = self.eigV.T
        if MODE_dynspec == "scipy":
            self.eigv, self.eigV = scsp.linalg.eigsh(
                self.sLp.astype(np.float64), k=howmany, which=which
            )
            self.eigV = self.eigV.T

    #
    def bin_eigV(self, which: int = 0):
        try:
            eigVbin = np.sign(np.where(self.eigV[which] == 0, +1, self.eigV[which]))
        except (AttributeError, IndexError):
            self.compute_k_eigvV(howmany = which + 1)
            eigVbin = np.sign(np.where(self.eigV[which] == 0, +1, self.eigV[which]))

        return eigVbin

    def bin_eigV_all(self):
        try:
            eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
        except (AttributeError, IndexError):
            self.compute_k_eigvV()
            eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
        return eigVbin
    
    def calc_fluct_Pinf(self, which: int = 0):
        eigV = self.bin_eigV(which)
        eV_p = np.count_nonzero(eigV >= 0)
        eV_n = self.N - eV_p

        self.eigV_fluct = abs(eV_p - eV_n) / self.N
        self.Pinf = np.min([eV_p, eV_n]) / self.N

        if hasattr(self, 'eigV_fluct_dict'):
            self.eigV_fluct_dict[which] = self.eigV_fluct
        else:
            self.eigV_fluct_dict = {which: self.Pinf}

        if hasattr(self, 'Pinf_dict'):
            self.Pinf_dict[which] = self.Pinf
        else:
            self.Pinf_dict = {which: self.Pinf}


    #
    def rescaled_signed_laplacian(self, MODE: str = "field"):
        if MODE == "field":
            self.resLp = self.sLp - self.eigv[0] * scsp.identity(self.N)
        elif MODE == "double":
            self.resLp = self.sLp - np.array([self.eigv[0]])
            new_eigv0 = scipy.linalg.eigvalsh(
                self.resLp.astype(np.float64), subset_by_index=[0, 0]
            )
            self.resLp = self.resLp - new_eigv0 * np.identity(self.N)

    #
    def lsp_selection(self, custom_list):
        if self.lsp_mode == "custom":
            self.lsp = np.array(custom_list)
        elif self.lsp_mode == "intervals":
            intervals = []
            tmp = max([vset["rsf"] for vset in custom_list])
            for vset in custom_list:
                if vset["kind"] == "log":
                    spacing_f = np.logspace
                    vset["start"] = np.log10(vset["start"])
                    vset["stop"] = np.log10(vset["stop"])
                elif vset["kind"] == "lin":
                    spacing_f = np.linspace
                intervals.append(  #
                    round_sigfig_n(  #
                        spacing_f(
                            vset["start"],
                            vset["stop"],
                            num=vset["num"],
                            endpoint=False,
                        ),
                        vset["rsf"],
                    )
                )
            self.lsp = (intervals := np.concatenate(intervals))
            while set(self.lsp).__len__() == intervals.__len__():
                tmp = tmp - 1
                self.lsp = np.round(self.lsp, tmp)
            tmp = tmp + 1
            self.lsp = np.round(intervals, tmp)

    def default_dict_lsp(self, num_low=3, num_at=6, num_high=3):
        d = (  #
            {
                "kind": "lin",
                "start": 0.001,
                "stop": self.p_c - self.p_c * num_at / 100,
                "num": num_low,
                "rsf": 1,
            },
            {
                "kind": "lin",
                "start": self.p_c - self.p_c * num_at / 100,
                "stop": self.p_c + self.p_c * num_at / 100,
                "num": num_at,
                "rsf": 3,
            },
            {
                "kind": "lin",
                "start": self.p_c + self.p_c * num_at / 100,
                "stop": 1,
                "num": num_high,
                "rsf": 1,
            },
        )
        return d
    #
    def dfs_list(self, node, visited, sign):
        if visited[node] or sign[node] <= 0:
            return 0
        visited[node] = True
        size = 1
        for neighbor in self.H[node]:
            if not visited[neighbor]:
                size += self.dfs_list(neighbor, visited, sign)
        return size
    #
    def cluster_distribution_list(self, sv = None):
        visited = [False] * len(self.H)
        distribution = {}
        if sv is None:
            try:
                sv = self.eigV[0]
            except AttributeError:
                self.compute_k_eigvV()
                sv = self.eigV[0]
        for node in range(len(self.H)):
            if not visited[node] and sv[node] > 0:
                size = self.dfs_list(node, visited, sv)
                distribution[size] = distribution.get(size, 0) + 1
        return distribution
    #
    def export_graph(self, MODE: str = "pickle"):
        if MODE == "pickle":
            pickle.dump(
                self.G,
                open(f"{self.graphfname}.pickle", "wb"),
                pickle.HIGHEST_PROTOCOL,
            )
        elif MODE == "gml":
            nx.write_gml(self.G, f"{self.graphfname}.gml")

    #
    def export_adj_bin(self, print_msg: bool = False) -> None:
        rowarr = [row[i:] for i, row in enumerate(self.Adj.todense())]
        exname = f"{self.expath}adj_{self.stdFname}.bin"
        if print_msg:
            print(f"exporting {exname}\n")
        with open(exname, "wb") as f:
            for i in range(len(rowarr)):
                rowarr[i].astype("float64").tofile(f)

    #
    def export_edgel(self):
        # TO BE FIXED
        a = list(self.H.edges(data="weight"))
        with open(r"src/LRGSG_package/tmp_stuff/prova.txt", "w") as fp:
            for item in a:
                # write each item on a new line
                fp.write("%s %s %s\n" % item)
            print("Done")


class FullyConnected(SignedGraph):
    def __init__(self, side1: int, anigemb: str = "sle"):
        self.side1 = side1
        self.G = nx.complete_graph(self.side1)
        super(FullyConnected, self).__init__(self.G)
        self.init_graph()
        self.pbc = True
        self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.esetG[len(self.esetG) // 2]: -1}
        self.DEFAULT_NEG_WEIGHTS_DICT_H = self.DEFAULT_NEG_WEIGHTS_DICT_G
        self.animation_graph_embedding = anigemb

    def init_graph(self):
        self.H = self.G
        self.upd_graph_matrices()

    def init_paths(self):
        self.sgpath = f"fc/"
        self.pltPath = f"data/plot/{self.sgpath}"
        self.datPath = f"data/{self.sgpath}"

    #
    def make_animation(
        self, fig, ax, frames, cmap: Union[str, Colormap] = "viridis"
    ):
        # I like to position my colorbars this way, but you don't have to
        #
        G_nodecol = frames[0]
        G_edgecol = [
            "b" if (e[2]["weight"] > 0) else "r"
            for e in self.G.edges(data=True)
        ]
        if self.animation_graph_embedding == "sle":
            pos = signed_spectral_layout(self.G)
        elif self.animation_graph_embedding == "circular":
            pos = nx.circular_layout(self.G)
        # nx.draw(G, ax=ax, pos=pos, edge_color=G_edgecol, node_color=G_nodecol, cmap='viridis')
        nodes = nx.draw_networkx_nodes(
            self.G, pos=pos, node_color=G_nodecol, cmap=cmap
        )
        nx.draw_networkx_edges(self.G, pos=pos, edge_color=G_edgecol)
        #
        cbar = fig.colorbar(nodes)
        # tx = ax.set_title('Frame 0')

        def animate(i):
            G_nodecol = frames[i]
            vmax = np.max(G_nodecol)
            vmin = np.min(G_nodecol)
            nx.draw_networkx_nodes(
                self.G, pos=pos, node_color=G_nodecol, cmap=cmap
            )
            cbar.mappable.set_clim(vmin, vmax)

        return animate
#
class Lattice2D(SignedGraph):
    #
    def __init__(
        self,
        side1: int = 32,
        geometry: str = DEFAULT_LATTICE2D_GEOMETRY,
        side2: int = 0,
        pbc: bool = True,
        fbc_val: float = 1.0,
        stdFnameSFFX: str = "",
        sgpath: str = "",
        **kwargs,
    ) -> None:
        try:
            self.geometry = geometry
            if geometry not in DEFLIST_LATTICE2D_GEOMETRIES:
                raise Lattice2DError(
                    """The selected geometry of the 2D lattice
                                     is not available. Setting it to 'squared' 
                                     for a 2d regular grid."""
                )
        except:
            self.geometry = DEFAULT_LATTICE2D_GEOMETRY
        self.side1 = side1
        if side2:
            self.side2 = side2
        else:
            if self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[0]:
                self.side2 = int(self.side1 * np.sqrt(3))
            elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[1]:
                self.side2 = self.side1
            elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[2]:
                self.side1 = int(self.side1 * np.sqrt(3))
                self.side2 = side1
        self.pbc = pbc
        self.fbc_val = fbc_val
        self.G = self.lattice_selection()
        self.sgpath = (
            f"{DEFAULT_LATTICE2D_PATHABBRV}{self.geometry}/"
            if not sgpath
            else sgpath
        )
        self.init_stdFname(stdFnameSFFX)
        super(Lattice2D, self).__init__(self.G, **kwargs)
        self.init_graph()
    #
    def init_graph(self):
        if self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[1]:
            self.posG = dict(zip(self.G, self.G))
            nx.set_node_attributes(self.G, values=self.posG, name="pos")
        self.upd_graph_matrices()
        self.neg_weights_dict = self.neg_weights_dicts_container(self)
    #
    class neg_weights_dicts_container:
        def __init__(self, lattice: SignedGraph):
            self.lattice = lattice
            self.midway_e = lattice.Ne // 2
            self.midway_H = lattice.N // 2 + lattice.side1 //2
            self.H_cent_edge = self.lattice.esetH[self.midway_e +1]
            #
            self.DEFAULT_NEG_WEIGHTS_DICT_H = {self.H_cent_edge: -1}
            try:
                self.DEFAULT_NEG_WEIGHTS_DICT_G = {lattice.invedge_map[self.H_cent_edge]: -1}
            except KeyError:
                print('KeyError')
            #
            # self.NEG_WEIGHTS_DICT_H_2ADJ = {lattice.esetH[self.midway_e]: -1, 
            #                                 lattice.esetH[self.midway_e+4]: -1}
            self.NEG_WEIGHTS_DICT_H_2CONT = {lattice.esetH[self.midway_e]: -1, 
                                            lattice.esetH[self.midway_e+2]: -1}
            #
            self.NEG_WEIGHTS_DICT_H_CROSS = self.get_neg_weights_dict_h_cross(self.midway_H)
            self.NEG_WEIGHTS_DICT_H_SQUARE = self.get_neg_weights_dict_h_square(self.midway_H)
            #
            self.flip_selection = np.random.choice(lattice.N, int(lattice.pflip * lattice.N))
            self.NEG_WEIGHTS_DICT_H_PCROSS = self.get_neg_weights_dict_h_pattern('cross')
            self.NEG_WEIGHTS_DICT_H_PSQUARE = self.get_neg_weights_dict_h_pattern('square')
            # _ = [graph[node][neighbor].update({'weight': -1}) for node in neighs_to_flip for neighbor in graph.neighbors(node)]
        #
        def get_neg_weights_dict_h_right(self, node: int):
            dictH = {
                (node, node+1): -1,
            }
            return dictH
        #
        def get_neg_weights_dict_h_square(self, node: int):
            dictH = {
                (node, node+1): -1,
                (node, node+self.lattice.side1): -1,
                (node + 1, node+self.lattice.side1 + 1): -1,
                (node + self.lattice.side1, node+self.lattice.side1 + 1): -1,
            }
            return dictH
        #
        def get_neg_weights_dict_h_rann(self, node: int):
            random_neighbor = random.choice(list(self.lattice.H.neighbors(node)))
            dictH = {
                (node, random_neighbor): -1
            }
            return dictH
        #
        def get_neg_weights_dict_h_cross(self, node: int):
            dictH = {
                (node-self.lattice.side1, node): -1,
                (node-1, node): -1,
                (node, node+1): -1,
                (node, node + self.lattice.side1): -1
            }
            return dictH
        #
        def get_neg_weights_dict_h_pattern(self, mode: str):
            if mode == "cross":
                merged_list = [item for i in self.flip_selection for item in self.get_neg_weights_dict_h_cross(i).items()]
            elif mode == "square":
                merged_list = [item for i in self.flip_selection for item in self.get_neg_weights_dict_h_square(i).items()]
            return dict(merged_list)
        #
        def get_neg_weights_dict_h_rball(self, R: int = 5):
            neighs_to_flip = get_neighbors_within_distance(self.lattice.H, self.midway_H, R)
            dictH = {(node, neighbor): -1 for node in neighs_to_flip for neighbor in self.lattice.H.neighbors(node)}
            return dictH
        #
        def get_neg_weights_dict_h_rball_EXT(self, R: int = 5):
            neighs_to_flip = get_neighbors_within_distance(self.lattice.H, self.midway_H, R)

            # Get the boundary nodes
            boundary_nodes = set()
            for node in neighs_to_flip:
                for neighbor in self.lattice.H.neighbors(node):
                    if neighbor not in neighs_to_flip:
                        boundary_nodes.add(node)
                        break

            # Create the dictionary with modified outer links
            dictH = {}
            for node in boundary_nodes:
                for neighbor in self.lattice.H.neighbors(node):
                    if neighbor not in neighs_to_flip:
                        dictH[(node, neighbor)] = -1
                    else:
                        dictH[(node, neighbor)] = 0

            return dictH
    #
    def init_stdFname(self, SFFX):
        if self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[0]:
            self.stdFname = DEFLIST_LATTICE2D_GEOABBRV[0]
        elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[1]:
            self.stdFname = DEFLIST_LATTICE2D_GEOABBRV[1]
        elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[2]:
            self.stdFname = DEFLIST_LATTICE2D_GEOABBRV[2]
        self.stdFname = self.stdFname + SFFX

    def lattice_selection(self, pbc = None) -> Graph:
        if pbc is None:
            pbc = self.pbc
        else:
            pbc = False
        if self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[0]:
            nxfunc = nx.triangular_lattice_graph
            self.p_c = 0.146
            self.r_c = np.sqrt(1.128/(np.pi*self.p_c))
            kwdict = {"with_positions": True}
        elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[1]:
            nxfunc = nx.grid_2d_graph
            self.p_c = 0.103
            self.r_c = np.sqrt(1.128/(np.pi*self.p_c))
            kwdict = {}
            self.syshape = (self.side1, self.side2)
        elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[2]:
            nxfunc = nx.hexagonal_lattice_graph
            self.p_c = 0.065
            self.r_c = np.sqrt(1.128/(np.pi*self.p_c))
            kwdict = {"with_positions": True}
        self.syshapePTH = (
            f"N={self.side1**2}/"
            if self.side1 == self.side2
            else f"L1={self.side1}_L2={self.side2}/"
        )
        return nxfunc(self.side1, self.side2, periodic=pbc, **kwdict)

    #
    def make_animation(self, fig, ax, frames):
        # I like to position my colorbars this way, but you don't have to

        cv0 = frames[0].reshape(self.syshape)
        im = ax.imshow(cv0)  # Here make an AxesImage rather than contour
        _, _, cbar = imshow_colorbar_caxdivider(im, ax)
        # tx = ax.set_title('Frame 0')

        def animate(i):
            arr = frames[i].reshape(self.syshape)
            vmax = np.max(arr)
            vmin = np.min(arr)
            im.set_data(arr)
            cbar.mappable.set_clim(vmin, vmax)
            # tx.set_text('Frame {0}'.format(i))
            # In this version you don't have to do anything to the colorbar,
            # it updates itself when the mappable it watches (im) changes

        return animate

class ErdosRenyi(SignedGraph):
    def __init__(self, n, p, sgpath: str = "", stdFname: str = "", **kwargs):
        self.G = nx.erdos_renyi_graph(n, p)
        self.sgpath = DEFAULT_ERGRAPH_SGPATH if not sgpath else sgpath
        self.syshapePTH = n
        self.stdFname = DEFAULT_ERGRAPH_ABBRV if not stdFname else stdFname
        super(ErdosRenyi, self).__init__(self.G, **kwargs)
        self.init_graph()
        self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.esetG[len(self.esetG) // 2]: -1}
        self.DEFAULT_NEG_WEIGHTS_DICT_H = self.DEFAULT_NEG_WEIGHTS_DICT_G

    def init_graph(self):
        self.H = self.G
        self.upd_graph_matrices()

