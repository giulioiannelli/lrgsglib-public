from .common import *
from .SignedGraph import SignedGraph

class ErdosRenyi(SignedGraph):
    def __init__(
            self, 
            n: int, 
            p: float, 
            sgpath: str = ER_SGPATH, 
            stdFnameSFFX: str = ER_STDFN,
            **kwargs):
        
        self.sgpath = ER_PHTABB if not sgpath else sgpath
        self.__init_stdFname__(stdFnameSFFX)
        self.__init_network__(n, p)
        super(ErdosRenyi, self).__init__(self.G, **kwargs)
    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = "er" + SFFX
    #
    def __init_network__(self, n, p):
        G = nx.erdos_renyi_graph(n, p)
        # Get the connected components (returns a list of sets of nodes)
        CC = nx.connected_components(G)
        # Find the largest connected component
        GC = max(CC, key=len)
        # Create a subgraph from the nodes in the largest component
        self.G = G.subgraph(GC).copy()
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={n}_p={p:.3g}"
    class nwContainer(dict):
        def __init__(self, er: SignedGraph, iterable=[], constant=None, 
                     **kwargs):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.er = er
            self.rd = self.er.GraphReprs
        #     #
        #     self['single'] = {g: [self.centedge[g]] for g in self.rd}
        #     self['singleZERR'] = {g: self.get_links_ZERR(
        #         self.centedge[g][0], g, self.l.geo) for g in self.rd}
        #     self['singleXERR'] = {g: self.get_links_XERR(
        #         self.centedge[g][0], g) for g in self.rd}
        #     self['rand'] = {g: [e for e in self.l.fleset[g]] 
        #                     for g in self.rd}
        #     self['randZERR'] = {g: self.get_rand_pattern('ZERR', on_g=g) 
        #                            for g in self.rd}
        #     self['randXERR'] = {g: self.get_rand_pattern('XERR', on_g=g) 
        #                            for g in self.rd}
        # #
        # def get_links_XERR(self, node: Any, on_g: str = L2D_ONREP):
        #     return [(node, nn) for nn in self.l.graph_neighbors(node, on_g)]