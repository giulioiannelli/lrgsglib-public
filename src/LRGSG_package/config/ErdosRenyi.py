from .nx_objects import *

class ErdosRenyi(SignedGraph):
    def __init__(self, n, p, sgpath: str = "", stdFname: str = "", **kwargs):
        self.G = nx.erdos_renyi_graph(n, p)
        self.sgpath = DEFAULT_ERGRAPH_SGPATH if not sgpath else sgpath
        self.syshapePth = n
        self.stdFname = DEFAULT_ERGRAPH_ABBRV if not stdFname else stdFname
        super(ErdosRenyi, self).__init__(self.G, **kwargs)
        self.init_graph()
        self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.esetG[len(self.esetG) // 2]: -1}
        self.DEFAULT_NEG_WEIGHTS_DICT_H = self.DEFAULT_NEG_WEIGHTS_DICT_G

    def init_graph(self):
        self.H = self.G
        self.upd_graph_matrices()