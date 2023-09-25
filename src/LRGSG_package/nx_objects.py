import networkx as nx
import numpy as np

from .LRGSG_const import *
from .LRGSG_errwar import *
from .LRGSG_plots import imshow_colorbar_caxdivider
from .LRGSG_utils import round_sigfig_n
from matplotlib.colors import Colormap
from networkx.classes.graph import Graph
from .nx_patches import signed_spectral_layout
from typing import Union


class SignedGraph():
    p_c = None
    lsp = None
    def __init__(self, G: Graph, lsp_mode: str = 'intervals'):
        self.G = G
        self.lsp_mode = lsp_mode
    #
    def init_sgraph(self):
        self.N = self.G.number_of_nodes()
        self.Ne = self.G.number_of_edges()
        self.esetG = list(self.G.edges())
        self.init_H_graph()
    #
    def init_H_graph(self):
        self.H = nx.convert_node_labels_to_integers(self.G)
        self.upd_H_graph()
    #
    def upd_H_graph(self):
        self.esetH = list(self.H.edges())
        self.posH = nx.get_node_attributes(self.H, 'pos')
        self.node_map = dict(zip(self.G, self.H))
        self.edge_map = dict(zip(self.G.edges(), self.H.edges()))
        self.number_of_negative_links()
    #
    def upd_G_graph(self):
        self.invnode_map = {v: k for k, v in self.node_map.items()}
        self.invedge_map = {v: k for k, v in self.edge_map.items()}
        self.G = nx.relabel_nodes(self.H, self.invnode_map)
        self.posG = nx.get_node_attributes(self.G, 'pos')
        self.esetG = list(self.G.edges())
        self.number_of_negative_links()
    #
    def number_of_negative_links(self):
        self.Ne_n = (np.array(list(
            nx.get_edge_attributes(self.H, 'weight').values())) < 0).sum()
    #
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

class FullyConnected(SignedGraph):
    def __init__(self, side1: int, anigemb: str = 'sle'):
        self.side1 = side1
        self.G = nx.complete_graph(self.side1)
        super().__init__(self.G)
        self.init_graph()
        self.pbc = True
        self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.esetG[len(self.esetG)//2]: -1}
        self.DEFAULT_NEG_WEIGHTS_DICT_H = self.DEFAULT_NEG_WEIGHTS_DICT_G
        self.animation_graph_embedding = anigemb
    
    def init_graph(self):
        self.H = self.G
        self.init_sgraph()
    
    def init_paths(self):
        self.lambdaPath = f"fc/"
        self.pltPath = f"data/plot/{self.lambdaPath}"
        self.datPath = f"data/{self.lambdaPath}"
    #
    def make_animation(self, fig, ax, frames, cmap: Union[str, Colormap] = 'viridis'):
        # I like to position my colorbars this way, but you don't have to
        #
        G_nodecol = frames[0]
        G_edgecol = ['b' if (e[2]['weight'] > 0) else 'r' 
                     for e in self.G.edges(data=True)]
        if self.animation_graph_embedding == 'sle':
            pos = signed_spectral_layout(self.G)
        elif self.animation_graph_embedding == 'circular':
            pos = nx.circular_layout(self.G)
        # nx.draw(G, ax=ax, pos=pos, edge_color=G_edgecol, node_color=G_nodecol, cmap='viridis')
        nodes = nx.draw_networkx_nodes(self.G, pos=pos, 
                                       node_color=G_nodecol, cmap=cmap)
        nx.draw_networkx_edges(self.G, pos=pos, 
                                       edge_color=G_edgecol)
        #
        cbar = fig.colorbar(nodes)
        # tx = ax.set_title('Frame 0')

        def animate(i):
            G_nodecol = frames[i]
            vmax = np.max(G_nodecol)
            vmin = np.min(G_nodecol)
            nx.draw_networkx_nodes(self.G, pos=pos, 
                                       node_color=G_nodecol, cmap=cmap)
            cbar.mappable.set_clim(vmin, vmax)
        return animate


#
class Lattice2D(SignedGraph):
    #
    def __init__(self, side1: int, geometry: str = DEFAULT_LATTICE2D_GEOMETRY, 
                 side2: int = 0, pbc: bool = True, fbc_val: float = 1.) -> None:
        try: 
            self.geometry = geometry
            if geometry not in DEFLIST_LATTICE2D_GEOMETRIES:
                raise Lattice2DError("""The selected geometry of the 2D lattice
                                     is not available. Setting it to 'squared' 
                                     for a 2d regular grid.""")
        except:
            self.geometry = self.DEFAULT_GEOMETRY
        self.side1 = side1
        if side2:
            self.side2 = side2
        else:
            if self.geometry == 'triangular':
                self.side2 = int(self.side1 * np.sqrt(3))
            elif self.geometry == 'hexagonal':
                self.side1 = int(self.side1 * np.sqrt(3))
                self.side2 = side1
            elif self.geometry == 'squared':
                self.side2 = self.side1
        self.pbc = pbc
        self.fbc_val = fbc_val
        self.G = self.lattice_selection()
        super().__init__(self.G)
        self.init_graph()
        self.init_paths()
        self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.esetG[len(self.esetG)//2]: -1}
        self.DEFAULT_NEG_WEIGHTS_DICT_H = {self.esetH[len(self.esetH)//2]: -1}
    #
    def init_graph(self):
        if self.geometry == 'squared':
            self.posG = dict(zip(self.G, self.G))
            nx.set_node_attributes(self.G, values=self.posG, name='pos')
        self.init_sgraph()
    #
    def init_paths(self):
        self.lambdaPath = f"l2d_{self.geometry}/"
        self.pltPath = f"data/plot/{self.lambdaPath}"
        self.datPath = f"data/{self.lambdaPath}"

    
    def lattice_selection(self, pbc=None) -> Graph:
        if pbc is None:
            pbc = self.pbc
        else:
            pbc = False
        if self.geometry == 'triangular':
            nxfunc = nx.triangular_lattice_graph
            self.p_c = 0.146
            kwdict = {'with_positions': True}
        elif self.geometry == 'squared':
            nxfunc = nx.grid_2d_graph
            self.p_c = 0.103
            kwdict = {}
            self.syshape = (self.side1, self.side2)
        elif self.geometry == 'hexagonal':
            nxfunc = nx.hexagonal_lattice_graph
            self.p_c = 0.065
            kwdict = {'with_positions': True}
        return nxfunc(self.side1, self.side2, periodic=pbc, **kwdict)
    
    #
    def make_animation(self, fig, ax, frames):
        # I like to position my colorbars this way, but you don't have to

        cv0 = frames[0].reshape(self.syshape)
        im = ax.imshow(cv0) # Here make an AxesImage rather than contour
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