from .common import *
from .funcs import *
from .SignedGraph import SignedGraph
from ..config.plotlib import Colormap
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
        self.upd_graphMtr()

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