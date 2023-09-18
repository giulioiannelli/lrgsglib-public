import networkx as nx
import numpy as np

class FullyConnected():
    p_c = None
    lsp = None
    def __init__(self, side1: int):
        self.side1 = side1
        self.init_graph()
        self.pbc = True
    
    def init_graph(self):
        self.G = nx.complete_graph(self.side1)
        self.N = self.G.number_of_nodes()
        self.esetG = list(self.G.edges())
        self.Ne = self.G.number_of_edges()
    
    def init_paths(self):
        self.lambdaPath = f"fc_{self.geometry}/"
        self.pltPath = f"data/plot/{self.lambdaPath}"
        self.datPath = f"data/{self.lambdaPath}"
    #
    def init_H_graph(self):
        self.upd_H_graph()
    #
    def upd_H_graph(self):
        self.H = nx.convert_node_labels_to_integers(self.G)
        self.node_map = dict(zip(self.G, self.H))
        self.edge_map = dict(zip(self.G.edges(), self.H.edges()))
        self.esetH = list(self.H.edges())
    #
    def upd_G_graph(self):
        self.invnode_map = {v: k for k, v in self.node_map.items()}
        self.invedge_map = {v: k for k, v in self.edge_map.items()}
        self.G = nx.relabel_nodes(self.H, self.invnode_map)
        self.esetG = list(self.G.edges())
    #
    def number_of_negative_links(self):
        self.Ne_n = (np.array(list(
            nx.get_edge_attributes(self.H, 'weight').values())) < 0).sum()
        return self.Ne_n
    #
    def make_animation(self, fig, ax, frames):
        # I like to position my colorbars this way, but you don't have to
        #
        G_nodecol = frames[0]
        G_edgecol = ['b' if (e[2]['weight'] > 0) else 'r' 
                     for e in self.G.edges(data=True)]
        pos = nx.circular_layout(self.G)
        # nx.draw(G, ax=ax, pos=pos, edge_color=G_edgecol, node_color=G_nodecol, cmap='viridis')
        nodes = nx.draw_networkx_nodes(self.G, pos=pos, 
                                       node_color=G_nodecol, cmap='viridis')
        nx.draw_networkx_edges(self.G, pos=pos, 
                                       edge_color=G_edgecol)
        #
        fig.colorbar(nodes)
        # tx = ax.set_title('Frame 0')

        def animate(i):
            G_nodecol = frames[i]
            vmax = np.max(G_nodecol)
            vmin = np.min(G_nodecol)
            nodes = nx.draw_networkx_nodes(self.G, pos=pos, 
                                       node_color=G_nodecol, cmap='viridis')
            # nodes.set_clim(vmin, vmax)
            fig.colorbar(nodes)
        return animate