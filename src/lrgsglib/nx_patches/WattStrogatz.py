from .common import *
from .SignedGraph import SignedGraph

class WattStrogatz(SignedGraph):
    def __init__(
            self, 
            n: int, 
            p: float, 
            sgpath: str = WS_SGPATH, 
            stdFnameSFFX: str = WS_STDFN,
            **kwargs):
        
        self.sgpath = WS_PHTABB if not sgpath else sgpath
        self.__init_stdFname__(stdFnameSFFX)
        self.__init_network__(n, p)
        super(WattStrogatz, self).__init__(self.G, **kwargs)
    
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = WS_PHTABB + SFFX
    #
    def __init_network__(self, n, p):
        self.G = nx.connected_watts_strogatz_graph(n, p)
        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={n}_p={p:.3g}"
    # class nwContainer(dict):
    #     def __init__(self, er: SignedGraph, iterable=[], constant=None, 
    #                  **kwargs):
    #         super().__init__(**kwargs)
    #         self.update((key, constant) for key in iterable)
    #         self.er = er
    #         self.rd = self.er.GraphReprs
    #         self.rNodeFlip = {g: random.sample(
    #                                 list(self.er.nodesIn[g]), 
    #                                 self.er.nflip
    #                             ) for g in self.rd}
    #         self['randXERR'] = {g: self.get_rand_pattern('XERR', on_g=g) 
    #                                for g in self.rd}
    #     #
    #     def get_links_XERR(self, node: Any, on_g: str = ER_ONREP):
    #         return [(node, nn) for nn in self.er.get_graph_neighbors(node, on_g)]
    #     #
    #     def get_rand_pattern(self, mode: str, on_g: str = ER_ONREP):
    #         match mode:
    #             case 'XERR':
    #                 if COUNT_XERR_PATTERNS:
    #                     patternList = [k for i in self.rNodeFlip[on_g] 
    #                                 for k in self.get_links_XERR(i, on_g)]
    #                 else:
    #                     tmplst = self.rNodeFlip[on_g]
    #                     grph = self.er.gr[on_g]
    #                     _ = 0
    #                     patternList = []
    #                     while _ < len(tmplst):
    #                         leval = [all([nnn['weight'] == -1 
    #                                     for nnn in grph[nn].values()])
    #                                     for nn in grph.neighbors(tmplst[_])]
    #                         if any(leval):
    #                             tmplst.pop(_)  # Removing the element
    #                         else:
    #                             glXERR = self.get_links_XERR(tmplst[_], 
    #                                                          on_g)
    #                             patternList.extend([k for k in glXERR])
    #                             _ += 1
    #         return patternList