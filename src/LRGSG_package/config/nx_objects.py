import os
import pickle
import random
import scipy

import networkx as nx
import numpy as np
import scipy.sparse as scsp

from .LRGSG_const import *
from .LRGSG_errwar import *
from .LRGSG_plots import imshow_colorbar_caxdivider
from .LRGSG_utils import round_sigfig_n, sum_tuples
from matplotlib.colors import Colormap
from networkx.classes.graph import Graph
from ..nx_patches import signed_spectral_layout, get_neighbors_within_distance
from scipy.sparse import csr_array
from typing import Union
from.SignedGraph import SignedGraph



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
        elif self.geometry == DEFLIST_LATTICE2D_GEOMETRIES[0]:
            self.posG = dict((node, (node[0] + 0.5 * (node[1] % 2), -node[1] * np.sqrt(3) / 2)) for node in self.G.nodes())
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
            self.NEG_WEIGHTS_DICT_H_PTRIA = self.get_neg_weights_dict_h_pattern('triangle')
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
        def get_neg_weights_dict_h_triangle(self, node: int):
            node1 = node
            node2 = random.sample(list(self.lattice.H.neighbors(node1)), 1)[0]
            common_neighbors = list(nx.common_neighbors(self.lattice.H, node1, node2))
            node3 = random.sample(common_neighbors, 1)[0]

            dictH = {
                (node1, node2): -1,
                (node2, node3): -1,
                (node1, node3): -1
            }
            return dictH
        #
        def get_neg_weights_dict_h_pattern(self, mode: str):
            if mode == "cross":
                merged_list = [item for i in self.flip_selection for item in self.get_neg_weights_dict_h_cross(i).items()]
            elif mode == "square":
                merged_list = [item for i in self.flip_selection for item in self.get_neg_weights_dict_h_square(i).items()]
            elif mode == "triangle":
                merged_list = [item for i in self.flip_selection for item in self.get_neg_weights_dict_h_triangle(i).items()]
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

