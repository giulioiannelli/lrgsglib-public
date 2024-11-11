from .common import *
from .SignedGraph import SignedGraph

class WattStrogatz(SignedGraph):
    def __init__(
            self, 
            n: int,
            k: int, 
            p: float, 
            sgpathn: str = WS_SGPATH, 
            stdFnameSFFX: str = WS_STDFN,
            **kwargs):
        
        self.sgpathn = WS_PHTABB if not sgpathn else sgpathn
        self.__init_stdFname__(stdFnameSFFX)
        self.__init_network__(n, k, p)
        super(WattStrogatz, self).__init__(self.G, **kwargs)
    
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = WS_PHTABB + SFFX
    #
    def __init_network__(self, n, k, p):
        self.G = nx.connected_watts_strogatz_graph(n, k, p)
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={n}_p={p:.3g}"
    