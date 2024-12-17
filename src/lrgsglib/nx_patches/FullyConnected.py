from .common import *
from .funcs import *
from .SignedGraph import SignedGraph
from ..config.plotlib import Colormap
class FullyConnected(SignedGraph):
    #
    def __init__(
        self, 
        N: int = FC_N,
        with_positions: bool = False,
        mode_positions: Union[str, Callable] = "circular",
        anigemb: str = "sle",
        **kwargs
    ) -> None:
        self.N = N
        self.with_positions = with_positions
        self.mode_positions = mode_positions
        self.__init_fullyconnected__()
        super(FullyConnected, self).__init__(self.G, **kwargs)
        self.animation_graph_embedding = anigemb
    #
    def __init_fullyconnected__(self) -> None:
        self.G = nx.complete_graph(self.N)
        if self.with_positions:
            match self.mode_positions:
                case "circular":
                    pos = nx.circular_layout(self.G)
                case "sspectral":
                    pos = signed_spectral_layout(self.G)
                case callable(mode) if callable(mode):
                    pos = mode(self.G)
                case _:
                    raise ValueError(f"Unsupported mode_positions: {self.mode_positions}")
            nx.set_node_attributes(self.G, pos, "pos")
    #
    def __repr__(self) -> str:
        return f"FullyConnected(N={self.N})"
    #
    def __str__(self) -> str:
        return f"FullyConnected(N={self.N})"
    #
    def compute_hopfield_edges(self, mode: str) -> None:
        if mode == "random":
            self.hopfield_edges = {
                (u, v): np.random.choice([-1, 1])
                for u, v in self.G.edges()
            }
        elif mode == "all+":
            self.hopfield_edges = {
                (u, v): 1
                for u, v in self.G.edges()
            }
        elif mode == "all-":
            self.hopfield_edges = {
                (u, v): -1
                for u, v in self.G.edges()
            }
        else:
            raise ValueError(f"Unsupported mode: {mode}")
    #
    def set_hopfield_edges(self, on_g: str = SG_GRAPH_REPR, **kwcompute_hopf_edges) -> None:
        if not hasattr(self, "hopfield_edges"):
            self.compute_hopfield_edges(**kwcompute_hopf_edges)
        nx.set_edge_attributes(self.gr[on_g], self.hopfield_edges, "weight")
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