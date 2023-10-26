import pickle
import random
import scipy

import networkx as nx
import numpy as np
import scipy.sparse as scsp

from .LRGSG_const import *
from .LRGSG_errwar import *
from .LRGSG_plots import imshow_colorbar_caxdivider
from .LRGSG_utils import round_sigfig_n
from matplotlib.colors import Colormap
from networkx.classes.graph import Graph
from .nx_patches import signed_spectral_layout
from scipy.sparse import csr_array
from typing import Union


class SignedGraph:
    p_c = None
    lsp = None
    DEFAULT_OUTDIR = "data/l2d_sq_ising/graphs/"

    def __init__(
        self,
        G: Graph,
        lsp_mode: str = "intervals",
        stdFname: str = "graph",
        import_on: bool = False,
        pflip: float = 0.0,
        expathc: str = ""
    ):
        self.lsp_mode = lsp_mode
        self.stdFname = stdFname
        self.import_on = import_on
        self.expath = self.DEFAULT_OUTDIR + expathc
        if import_on:
            self.G = self.__init_graph_fromfile__()
        else:
            self.G = G
            self.pflip = pflip
        self.init_sgraph()

    #
    def __init_graph_fromfile__(self):
        return pickle.load(open(f"{self.expath}{self.stdFname}.pickle", "rb"))

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
        self.esetG = list(self.G.edges())
        self.init_n_nodes_edges()
        if self.import_on:
            self.upd_H_graph()
            self.nflip = self.Ne_n
            self.pflip = self.nflip / self.Ne
            self.upd_graph_matrices() 
            self.randsample = np.where(np.array(self.Adj.todense()).flatten() < 0)
        else:
            self.nflip = int(self.pflip * self.Ne)
            self.randsample = random.sample(range(self.Ne), self.nflip) 
            self.init_weights()
            self.upd_H_graph()

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
                neg_weights_dict = self.DEFAULT_NEG_WEIGHTS_DICT_G
            nx.set_edge_attributes(
                self.G, values=neg_weights_dict, name="weight"
            )
            self.upd_H_graph()
        elif on_graph == "H":
            if neg_weights_dict is None:
                neg_weights_dict = self.DEFAULT_NEG_WEIGHTS_DICT_H
            nx.set_edge_attributes(
                self.H, values=neg_weights_dict, name="weight"
            )
            self.upd_G_graph()
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
        neg_weights = {
            e: -1 for i, e in enumerate(eset) if i in self.randsample
        }
        self.flip_sel_edges(neg_weights_dict=neg_weights, on_graph=on_graph)

    #
    def compute_laplacian_spectrum(self, MODE_lapspec: str = "numpy") -> None:
        if MODE_lapspec == "networkx":
            self.slspectrum = nx.laplacian_spectrum(self.system.G)
        elif MODE_lapspec == "numpy":
            self.slspectrum = np.linalg.eigvalsh(self.sLp.toarray())

    #
    def compute_k_eigvV(self, MODE_dynspec: str = "scipy", howmany: int = 1):
        if MODE_dynspec == "scipy":
            self.eigv, self.eigV = scsp.linalg.eigsh(
                self.sLp.astype(np.float64), k=howmany, which="SM"
            )
            self.eigV = self.eigV.T

    #
    def bin_eigV(self, which=0):
        eigVbin = np.sign(self.eigV[which])
        eigVbin[eigVbin == 0] = +1
        return eigVbin
    def bin_eigV_all(self):
        eigVbin = np.sign(self.eigV)
        eigVbin[eigVbin == 0] = +1
        return eigVbin
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
    def export_graph(self, MODE: str = "pickle", expath: str = DEFAULT_OUTDIR):
        fname = f"{expath}N={self.N:d}/{self.stdFname}"
        if MODE == "pickle":
            pickle.dump(
                self.G, open(f"{fname}.pickle", "wb"), pickle.HIGHEST_PROTOCOL
            )
        elif MODE == "gml":
            nx.write_gml(self.G, f"{fname}.gml")
    #
    def export_adj_bin(self, expath: str = DEFAULT_OUTDIR):
        rowarr = [row[i:] for i, row in enumerate(self.Adj.todense())]
        with open(f"{expath}N={self.N:d}/adj_{self.stdFname}.bin", "wb") as f:
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
        self.lambdaPath = f"fc/"
        self.pltPath = f"data/plot/{self.lambdaPath}"
        self.datPath = f"data/{self.lambdaPath}"

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
        side1: int,
        geometry: str = DEFAULT_LATTICE2D_GEOMETRY,
        side2: int = 0,
        pbc: bool = True,
        fbc_val: float = 1.0,
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
            self.geometry = self.DEFAULT_GEOMETRY
        self.side1 = side1
        if side2:
            self.side2 = side2
        else:
            if self.geometry == "triangular":
                self.side2 = int(self.side1 * np.sqrt(3))
            elif self.geometry == "hexagonal":
                self.side1 = int(self.side1 * np.sqrt(3))
                self.side2 = side1
            elif self.geometry == "squared":
                self.side2 = self.side1
        self.pbc = pbc
        self.fbc_val = fbc_val
        self.G = self.lattice_selection()
        super(Lattice2D, self).__init__(self.G, **kwargs)
        self.init_graph()
        self.init_paths()
        self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.esetG[len(self.esetG) // 2]: -1}
        self.DEFAULT_NEG_WEIGHTS_DICT_H = {self.esetH[len(self.esetH) // 2]: -1}

    #
    def init_graph(self):
        if self.geometry == "squared":
            self.posG = dict(zip(self.G, self.G))
            nx.set_node_attributes(self.G, values=self.posG, name="pos")
        self.upd_graph_matrices()

    #
    def init_paths(self):
        if self.geometry == "triangular":
            self.stdFname = "trLattice" + f"_p={self.pflip:.3g}"
        elif self.geometry == "squared":
            self.stdFname = "sqLattice" + f"_p={self.pflip:.3g}"
        elif self.geometry == "hexagonal":
            self.stdFname = "hxLattice" + f"_p={self.pflip:.3g}"
        self.lambdaPath = f"l2d_{self.geometry}/"
        self.pltPath = f"data/plot/{self.lambdaPath}"
        self.datPath = f"data/{self.lambdaPath}"

    def lattice_selection(self, pbc=None) -> Graph:
        if pbc is None:
            pbc = self.pbc
        else:
            pbc = False
        if self.geometry == "triangular":
            nxfunc = nx.triangular_lattice_graph
            self.p_c = 0.146
            kwdict = {"with_positions": True}
        elif self.geometry == "squared":
            nxfunc = nx.grid_2d_graph
            self.p_c = 0.103
            kwdict = {}
            self.syshape = (self.side1, self.side2)
        elif self.geometry == "hexagonal":
            nxfunc = nx.hexagonal_lattice_graph
            self.p_c = 0.065
            kwdict = {"with_positions": True}
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
