from .objects import *
from scipy.sparse import spdiags
from scipy.sparse import identity as scsp_identity
from scipy.sparse.linalg import eigsh as scsp_eigsh
class SignedGraph:
    sgpath = "custom_graph"
    syshapePth = ""
    stdFname = ""

    def __init__(
        self,
        G: Graph, pflip: float = 0.0,
        import_on: bool = False,
        init_nw_dict: bool = False,
        init_weights_val: float = 1,
        make_dir_tree: bool = True,
        on_graph: str = "",
        expOutdir: str = "",
        dataOutdir: str = "",
        plotOutdir: str = "",
        seed: int = None,
    ):
        if seed:
            random.seed(seed)
            np.random.seed(seed)
        self.GraphReprDict = {}
        self.nodeMap = {}
        self.edgeMap = {}
        self.eset = {}
        self.fleset = {}
        self.nodesIn = {}
        #
        self.__init_paths__(
            dataOutdir=dataOutdir, 
            plotOutdir=plotOutdir, 
            expOutdir=expOutdir
        )
        self.__make_dirs__(make_dir_tree)
        # pflip check
        if not is_in_range(pflip, LB_PFLIP, UB_PFLIP):
            raise ValueError(SG_ERRMSG_PFLIP)
        else:
            self.pflip = pflip
        self.onGraph = on_graph or SG_GRAPH_REPR
        self.init_weights_val = init_weights_val
        self.slspectrum = None
        self.import_on = import_on
        self.stdFname = self.stdFname + f"_p={self.pflip:.3g}"
        self.graphPath = pth_join(self.expOutdir, self.stdFname)
        #
        if import_on:
            self.G = self.__init_graph_fromfile__()
        else:
            self.G = G
        self.__init_reprdict__()
        self.__init_sgraph__()
        if init_nw_dict:
            if hasattr(self, 'nwContainer'):
                self.nwDict = self.nwContainer(self)
            else:
                raise AttributeError(SG_ERRMSG_NW_DICT)
    #
    @property
    def gr(self):
        return self.GraphReprDict
    @property
    def adj(self):
        return self.Adj
    @property
    def lap(self):
        return self.Lap
    @property
    def slp(self):
        return self.sLp
    #
    def __init_paths__(
        self, dataOutdir: str = "", plotOutdir: str = "", expOutdir: str = ""
    ):
        self.dataOutdir = dataOutdir or DIR_DAT
        self.plotOutdir = plotOutdir or DIR_PLT
        #
        self.datPath = f"{self.dataOutdir}{self.sgpath}"
        self.pltPath = f"{self.dataOutdir}{self.plotOutdir}{self.sgpath}"
        #
        self.DEFAULT_GRAPHDIR = pth_join(self.datPath, DIR_GRAPH)
        self.DEFAULT_ISINGDIR = pth_join(self.datPath, DIR_ISING)
        self.DEFAULT_VOTERDIR = pth_join(self.datPath, DIR_VOTER)
        self.DEFAULT_LRGSGDIR = pth_join(self.datPath, DIR_LRGSG)
        self.DEFAULT_PHTRADIR = pth_join(self.datPath, DIR_PHTRA)
        self.DEFAULT_SPECTDIR = pth_join(self.datPath, DIR_SPECT)
        #
        self.expOutdir = expOutdir or pth_join(self.DEFAULT_GRAPHDIR, self.syshapePth)
        #
        self.isingpath = pth_join(self.DEFAULT_ISINGDIR, self.syshapePth)
        self.voterpath = pth_join(self.DEFAULT_VOTERDIR, self.syshapePth)
        self.lrgsgpath = pth_join(self.DEFAULT_LRGSGDIR, self.syshapePth)
        self.phtrapath = pth_join(self.DEFAULT_PHTRADIR, self.syshapePth)
        self.spectpath = pth_join(self.DEFAULT_SPECTDIR, self.syshapePth)
        self.dirMakeList = [self.expOutdir, self.isingpath, self.voterpath,
                self.phtrapath, self.lrgsgpath, self.spectpath]
    #
    def __make_dirs__(self, exist_ok: bool = True):
        for _ in self.dirMakeList: os.makedirs(_, exist_ok=exist_ok)
    #
    def __init_graph_fromfile__(self):
        return pk.load(open(self.graphPath + PKL, "rb"))
    #
    def __init_reprdict__(self):
        self.GraphReprDict[self.onGraph] = self.G
        for grd in ['H', 'I', 'J', 'K', 'L', 'M', 'N']:
            try:
                self.GraphReprDict[grd] = getattr(self, grd)
            except AttributeError:
                pass
        self.GraphReprs = list(self.GraphReprDict.keys())
    #
    def __init_weights__(self, values: Union[float, List] = 1) -> None:
        nx.set_edge_attributes(self.gr[self.onGraph], values, 'weight')
        self.upd_GraphRepr_All(self.onGraph)
    #
    def __init_nNodesEdges__(self) -> None:
        self.N = self.gr[self.onGraph].number_of_nodes()
        self.Ne = self.gr[self.onGraph].number_of_edges()
        self.nodesIn[self.onGraph] = list(self.gr[self.onGraph].nodes())
    #
    def __init_sgraph__(self, init_weights_val: Union[float, List] = 1) -> None:
        self.__init_nNodesEdges__()
        if self.import_on:
            self.neflip = self.Ne_n
            self.pflip = self.neflip / self.Ne
            edgesWithData = self.gr[self.onGraph].edges(data=True)
            self.fleset[self.onGraph] = [(u, v) for u, v, _ in edgesWithData 
                                         if _.get('weight', 1) < 0]
        else:
            self.neflip = int(self.pflip * self.Ne)
            self.nflip = int(self.pflip * self.N)
            self.__init_weights__(init_weights_val)
            edges = self.gr[self.onGraph].edges
            self.eset[self.onGraph] = list(edges())
            self.fleset[self.onGraph] = random.sample(self.eset[self.onGraph], 
                                                      self.neflip)
        self.upd_GraphRepr_All(self.onGraph)
        self.upd_graphMtr()
    #
    def upd_Degree(self, on_graph: str = SG_GRAPH_REPR):
        self.degrees = list(dict(self.gr[on_graph].degree).values())
    #
    def upd_Nen(self, on_graph: str = SG_GRAPH_REPR):
        self.Ne_n = len(self.fleset[on_graph])
    #
    def get_adjMtrx(self, on_graph: str = SG_GRAPH_REPR, weight: str = 'weight',
                    format: str = 'csr') -> csr_array:
        return nx.to_scipy_sparse_array(self.GraphReprDict[on_graph], 
                                        weight=weight, format=format)
    #
    def zip_reprNodes(self, x, on_graph: str = SG_GRAPH_REPR):
        graph = self.gr[on_graph]
        return dict(zip(graph, self.gr[x]))
    #
    def zip_reprEdges(self, x, on_graph: str = SG_GRAPH_REPR):
        graph = self.gr[on_graph]
        return dict(zip(graph.edges(), self.gr[x].edges()))
    #
    def upd_NodeMap(self, on_graph: str = SG_GRAPH_REPR):
        self.nodeMap[on_graph] = {x: 
            {v: k for k, v in self.zip_reprNodes(x, on_graph).items()} 
            for x in self.GraphReprs if x != on_graph}
    #
    def upd_EdgeMap(self, on_graph: str = SG_GRAPH_REPR):
        graph = self.gr[on_graph]
        egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
        self.edgeMap[on_graph] = {x: {v: k for k, v in egraph_to(x).items()} 
                                   for x in self.GraphReprs if x != on_graph}
    #
    def upd_GraphRelabel(self, on_graph: str = SG_GRAPH_REPR, to_graph: str = SG_GRAPHINT_REPR):
        self.gr[to_graph] = nx.relabel_nodes(self.gr[on_graph], self.nodeMap[to_graph][on_graph])
    #
    def upd_ReprMaps(self, on_graph: str = SG_GRAPH_REPR):
        graph = self.gr[on_graph]
        ngraph_to = lambda x: dict(zip(graph, self.gr[x]))
        self.nodeMap[on_graph] = {x: {v: k for k, v in ngraph_to(x).items()} 
                                   for x in self.GraphReprs if x != on_graph}
        egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
        self.edgeMap[on_graph] = {x: {v: k for k, v in egraph_to(x).items()} 
                                   for x in self.GraphReprs if x != on_graph}
    #
    def upd_GraphRepr(self, on_graph: str = SG_GRAPH_REPR):
        self.upd_ReprMaps(on_graph)
    #
    def upd_GraphRepr_All(self, on_graph: str = SG_GRAPH_REPR, also_itself: bool = True):
        if also_itself:
            self.upd_GraphRepr(on_graph)
        for i in self.GraphReprs:
            if i != on_graph:
                self.upd_GraphRepr(i)
                self.upd_GraphRelabel(on_graph, i)
                edges = self.gr[i].edges
                self.nodesIn[i] = list(self.gr[i].nodes())
                self.eset[i] = list(edges())
                self.fleset[i] = [(u, v) for u, v, _ in edges(data=True) 
                                         if _.get('weight', 1) < 0]
                self.upd_Nen(i)
    #
    def get_degMtrx(self, A: csr_array, fmt: str = 'csr') -> csr_array:
        return csr_array(spdiags(A.sum(axis=1), 0, *A.shape, format="csr"))

    #
    def absolute_degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(spdiags(abs(A).sum(axis=1), 0, *A.shape, format="csr"))

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
    def upd_graphMtr(self, on_graph: str = SG_GRAPH_REPR):
        self.Adj = self.get_adjMtrx(on_graph=on_graph)
        self.Deg = self.get_degMtrx(self.Adj)
        self.sDeg = self.absolute_degree_matrix(self.Adj)
        self.Lap = self.laplacian_matrix()
        self.sLp = self.signed_laplacian()

    def get_edge_weight(self, u: Any, v: Any, on_graph: str = SG_GRAPH_REPR):
        return self.gr[on_graph].get_edge_data(u, v)['weight']
    #
    def flip_sel_edges(self, links: List[Any] = [], on_graph: str = SG_GRAPH_REPR):
        """Flips specific edges of a graph G."""
        #
        neg_weights_dict = {}
        if links:
            neg_weights_dict = {
                (u, v): -1 * self.get_edge_weight(u, v)
                for u, v in links
            }
        nx.set_edge_attributes(
            self.gr[on_graph], values=neg_weights_dict, name='weight'
        )
        self.upd_GraphRepr_All(on_graph)
        self.upd_graphMtr(on_graph)

    #
    def check_pflip(self):
        if self.neflip < 1: raise NflipError(SG_ERRMSG_NFLIP)
    #
    def flip_random_fract_edges(self, on_graph: str = SG_GRAPH_REPR):
        """Flips a fraction p of edges (+1 to -1) of a graph G."""
        #
        try:
            self.check_pflip()
            self.flip_sel_edges(self.fleset[on_graph], on_graph)
        except NflipError:
            pass
    #
    def unflip_all(self, on_graph: str = SG_GRAPH_REPR):
        self.flip_sel_edges(1, on_graph)

    #
    def compute_k_eigvV(
        self,
        howmany: int = 1,
        MODE_dynspec: str = "scipy",
        which: str = "SM",
        typf: type = np.float64,
    ):
        if MODE_dynspec == "numpy" or howmany == self.N:
            self.eigv, self.eigV = np.linalg.eigh(
                self.sLp.astype(typf).todense()
            )
            self.eigV = self.eigV.T
        if MODE_dynspec.startswith("scipy"):
            # scsp_eigshMode = MODE_dynspec.split('_')
            # if len(scsp_eigshMode) == 1:
            #     scsp_eigshMode = 'normal'
            #     sigma = None
            # else:
            #     scsp_eigshMode = scsp_eigshMode[-1]
            #     sigma = 0
            # self.eigv, self.eigV = scsp_eigsh(
            #     self.sLp.astype(typf), k=howmany, which=which, sigma=sigma, mode=scsp_eigshMode
            # )
            self.eigv, self.eigV = scsp_eigsh(
                self.sLp.astype(typf), k=howmany, which=which, mode="caley"
            )
            self.eigV = self.eigV.T

    #
    def bin_eigV(self, which: int = 0):
        try:
            eigVbin = np.sign(
                np.where(self.eigV[which] == 0, +1, self.eigV[which])
            )
        except (AttributeError, IndexError):
            self.compute_k_eigvV(howmany=which + 1)
            eigVbin = flip_to_positive_majority(np.sign(
                    np.where(self.eigV[which] == 0, +1, self.eigV[which])
                )
            )

        return eigVbin
    
    def compute_rbim_energy_eigV(self, which: int = 0, on_graph: str = SG_GRAPH_REPR):
        spins = self.bin_eigV(which)
        edges = self.GraphReprDict[on_graph].edges(data=True)
        Aij = np.array([(u, v) for u, v, data in edges])
        Jij = np.array([data['weight'] for u, v, data in edges])
        # Compute the product of weights and spins for each edge
        Eij = Jij * spins[Aij[:, 0]] * spins[Aij[:, 1]]
        # Sum the product terms
        E = -np.sum(Eij)
        return E/self.N

    def bin_eigV_all(self):
        try:
            eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
        except (AttributeError, IndexError):
            self.compute_k_eigvV()
            eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
        return eigVbin

    def calc_Pinf(self, which: int = 0, on_graph: str = SG_GRAPH_REPR):
        cd = self.cluster_sizes(which, on_graph)
        size = cd[0]
        # eV_p = np.count_nonzero(eigV >= 0)
        # eV_n = self.N - eV_p
        # self.eigV_fluct = abs(eV_p - eV_n) / self.N
        self.Pinf = size / self.N
        # if hasattr(self, "eigV_fluct_dict"):
        #     self.eigV_fluct_dict[which] = self.eigV_fluct
        # else:
        #     self.eigV_fluct_dict = {which: self.Pinf}
        if hasattr(self, "Pinf_dict"):
            self.Pinf_dict[which] = self.Pinf
        else:
            self.Pinf_dict = {which: self.Pinf}

    #
    def rescaled_signed_laplacian(self, MODE: str = "field"):
        if MODE == "field":
            self.resLp = self.sLp - self.eigv[0] * scsp_identity(self.N)
        elif MODE == "double":
            from scipy.linalg import eigvalsh

            self.resLp = self.sLp - np.array([self.eigv[0]])
            new_eigv0 = eigvalsh(
                self.resLp.astype(np.float64), subset_by_index=[0, 0]
            )
            self.resLp = self.resLp - new_eigv0 * np.identity(self.N)


    def graph_neighbors(self, node, on_graph: str = SG_GRAPH_REPR):
        graph = self.GraphReprDict[on_graph]
        return list(graph.neighbors(node))

    def get_node_attr(self, attr, on_graph: str = SG_GRAPH_REPR):
        node_attr = nx.get_node_attributes(self.GraphReprDict[on_graph], attr)
        node_attrs = [
            node_attr[node] for node in self.GraphReprDict[on_graph].nodes()
        ]
        return node_attrs

    def load_vec_on_nodes(self, vec: np.ndarray, attr: str, on_graph: str = SG_GRAPH_REPR):
        vecNodeAttr = {
            nd: v for v, nd in zip(vec, self.GraphReprDict[on_graph].nodes)
        }
        nx.set_node_attributes(self.GraphReprDict[on_graph], vecNodeAttr, attr)
    def load_eigV_on_graph(self, which: int = 0, on_graph: str = SG_GRAPH_REPR, 
                           binarize: bool = False):
        try:
            eigV = self.eigV[which]
        except (IndexError, AttributeError):
            self.compute_k_eigvV(howmany=which + 1)
            eigV = self.eigV[which]
        if binarize:
            eigV = self.bin_eigV(which=which)
        eigVNodeAttr = {
            nd: v for v, nd in zip(eigV, self.GraphReprDict[on_graph].nodes)
        }
        nx.set_node_attributes(
            self.GraphReprDict[on_graph], eigVNodeAttr, f"eigV{which}"
        )

    def group_nodes_by_kv(self, k, val, on_graph: str = SG_GRAPH_REPR):
        G = self.GraphReprDict[on_graph]
        G_yes, G_no = G.copy(), G.copy()
        for node, v in G.nodes(data=k):
            if v == val:
                G_yes.remove_node(node)
            else:
                G_no.remove_node(node)
        return G_yes, G_no

    def cluster_sizes(self, which: int = 0, on_graph: str = SG_GRAPH_REPR, binarize: bool = True):
        self.load_eigV_on_graph(which, on_graph, binarize)
        G_pos, G_neg = self.group_nodes_by_kv(f"eigV{which}", -1, on_graph)
        clusters = list(nx.connected_components(G_neg))
        clusterLen = sorted(list(map(lambda x: len(x), clusters)), reverse=True)
        return clusterLen
    
    def cluster_distribution(self, which: int = 0, 
                                  on_graph: str = SG_GRAPH_REPR,
                                  binarize: bool = True):
        clusterLen = self.cluster_sizes(which, on_graph, binarize)
        distNeg = {
            size: clusterLen.count(size) for size in set(clusterLen)
        }
        # self.load_eigV_on_graph(which, on_graph, binarize)
        # _, G_neg = self.group_nodes_by_kv(f"eigV{which}", -1, on_graph)
        # clusters = sorted(list(nx.connected_components(G_neg)), 
        #               key=lambda x: len(x), reverse=True)
        # clNeg = list(map(lambda x: len(x), clusters))
        # distNeg = {
        #     size: clNeg.count(size) for size in set(clNeg)
        # }
        return distNeg
    def classify_ferroAntiferro_regions(self, attr_str: str = 's', on_graph: str = SG_GRAPH_REPR):
        antiGroup = []
        ferroGroup = []
        graph = self.GraphReprDict[on_graph]
        for node, att in graph.nodes(data=attr_str):
            neighbors = graph.neighbors(node)
            if all(graph.nodes[n][attr_str] == -att for n in neighbors):
                antiGroup.append(node)
            elif all(graph.nodes[n][attr_str] == att for n in neighbors):
                ferroGroup.append(node)
        
        return ferroGroup, antiGroup
    
    # def classify_disorder_regions(self, attr_str: str = 's', on_graph: str = SG_GRAPH_REPR):
    #     n_negnei = {}
    #     graph = self.GraphReprDict[on_graph]
    #     for node, att in graph.nodes(data=attr_str):
    #         neighbors = graph.neighbors(node)
    #         condition = (graph.nodes[n][attr_str] == -att for n in neighbors)
    #         if all(condition):
    #             n_negnei['all'].append(node)
    #         else any(condition):
    #             ferroGroup.append(node)
        
    #     return antiGroup, ferroGroup
    #
    def export_graph(self, MODE: str = "pickle"):
        if MODE == "pickle":
            pk.dump(
                self.G,
                open(f"{self.graphPath}.pkl", "wb"),
                pk.HIGHEST_PROTOCOL,
            )
        elif MODE == "gml":
            nx.write_gml(self.G, f"{self.graphPath}.gml")

    #
    def export_adj_bin(self, print_msg: bool = False) -> None:
        rowarr = [row[i:] for i, row in enumerate(self.Adj.todense())]
        exname = f"{self.expOutdir}adj_{self.stdFname}.bin"
        if print_msg:
            print(f"exporting {exname}\n")
        with open(exname, "wb") as f:
            for i in range(len(rowarr)):
                rowarr[i].astype("float64").tofile(f)
    #
    def export_edgel_bin(self, on_graph: str = SG_GRAPH_REPR, expoName: str = "",
                         print_msg: bool = False, mode: str = 'numpy') -> None:
        edges = self.GraphReprDict[on_graph].edges(data='weight')
        exname = expoName or self.stdFname
        match mode:
            case 'numpy':
                edge_array = np.array(list(edges), 
                                 dtype=[('i', np.uint64), #these has to be changed based on the repr
                                        ('j', np.uint64), #these has to be changed based on the repr
                                        ('w_ij', np.float64)])
                edge_array.tofile(f"{self.expOutdir}edgelist_{exname}.bin")
            case 'struct':
                with open(f"{self.expOutdir}edgelist_{exname}.bin", "wb") as f:
                    for edge in edges:
                        f.write(struct.pack("QQd", edge[0], edge[1], edge[2]))

    #
    def get_edge_color(
        self,
        on_graph: str = "G",
        pec: ColorType = "blue",
        nec: ColorType = "red",
    ):
        def map_values(value):
            if value == -1:
                return nec
            elif value == 1:
                return pec
            else:
                return value

        arr = nx.get_edge_attributes(self.GraphReprDict[on_graph], 'weight')
        return list(map(map_values, arr.values()))
        # return list(map(lambda x: pec if x == 1 else nec, nx.get_edge_attributes(self.GraphReprDict[on_graph], 'weight').values()))

    #
    def get_node_pos(self, on_graph: str = "G"):
        return nx.get_node_attributes(self.GraphReprDict[on_graph], "pos")