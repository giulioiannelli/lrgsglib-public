from .common import *
class SignedGraph:
    sgpathn = "signed_graph"
    syshapePth = ""
    stdFname = ""
    #
    def __init__(self,
        G: Graph, 
        pflip: float = SG_PFLIP,
        import_on: bool = SG_IMPORT_ON,
        init_nw_dict: bool = SG_INIT_NW_DICT,
        init_weights_val: float = SG_INIT_WVAL,
        xprt_modeL: str = SG_XPRT_MODE,
        mprt_modeL: str = SG_MPRT_MODE,
        make_dir_tree: bool = True,
        path_data: Path = None,
        on_g: str = "",
        expOut: str = "",
        dataOut: str = "",
        plotOut: str = "",
        seed: int = None,
    ):
        self.GraphReprDict = {}
        self.nodeMap = {}
        self.edgeMap = {}
        self.eset = {}
        self.fleset = {}
        self.lfeset = {}
        self.nodes_in = {}
        #
        self.seed = seed or int(time.time() * 1000) + os.getpid()
        random.seed(seed)
        np.random.seed(seed)
        #
        self.__init_paths__(
            path_data=path_data,
            dataOut=dataOut, 
            plotOut=plotOut, 
            expOut=expOut
        )
        if make_dir_tree:
            self.__make_dirs__()
        # pflip check
        if not is_in_range(pflip, LB_PFLIP, UB_PFLIP):
            raise ValueError(SG_ERRMSG_PFLIP)
        else:
            self.pflip = pflip
            self.pEqStr = f"p={self.pflip:.3g}"
        self.onGraph = on_g or SG_GRAPH_REPR
        self.init_weights_val = init_weights_val
        self.slspectrum = None
        self.import_on = import_on
        self.stdFname = '_'.join([self.stdFname, self.pEqStr])
        self.graphPath = pth_join(self.expOut, self.stdFname)
        #
        if import_on:
            self.G = self.__init_graph_fromfile__(mprt_modeL)
        else:
            self.G = G
        self.__init_reprdict__()
        self.__init_sgraph__()
        if init_nw_dict:
            if hasattr(self, 'nwContainer'):
                self.nwDict = self.nwContainer(self)
            else:
                raise AttributeError(SG_ERRMSG_NW_DICT)
        self.__make_graphCl_util__()
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
    @property
    def degm(self):
        return self.Deg
    @property
    def sdeg(self):
        return self.sDeg
    @property
    def gcl(self):
        return self.gclutil
    #
    def __init_paths__(
        self, path_data: Path = None, path_plot: Path = None, path_export: Path = None,
        dataOut: str = "", plotOut: str = "", expOut: str = ""
    ):
        #
        self.path_data = path_data or PATHDATA
        self.path_plot = path_plot or PATHPLOT
        self.path_sgdata = self.path_data / Path(self.sgpathn)
        #
        self.path_export = self.path_sgdata / Path(PATHNGRPH, self.syshapePth)
        self.path_ising = self.path_sgdata / Path(PATHNISNG, self.syshapePth)
        self.path_voter = self.path_sgdata / Path(PATHNVM, self.syshapePth)
        self.path_lrgsg = self.path_sgdata / Path(PATHNLRGS, self.syshapePth)
        self.path_phtra = self.path_sgdata / Path(PATHNPHTR, self.syshapePth)
        self.path_spect = self.path_sgdata / Path(PATHNSPEC, self.syshapePth)
        #
        self.dirMakeList = [self.path_export,
            self.path_ising,
            self.path_voter,
            self.path_lrgsg,
            self.path_phtra,
            self.path_spect]
        #
        self.dataOut = dataOut or str(self.path_data)
        self.plotOut = plotOut or str(self.path_plot)
        self.expOut = expOut or str(self.path_export)
        #
        self.sgdatpath = str(self.path_sgdata)
        self.isingpath = str(self.path_ising)
        self.voterpath = str(self.path_voter)
        self.lrgsgpath = str(self.path_lrgsg)
        self.phtrapath = str(self.path_phtra)
        self.spectpath = str(self.path_spect)
    #
    def __make_dirs__(self, exist_ok: bool = True):
        for _ in self.dirMakeList: os.makedirs(_, exist_ok=exist_ok)
    #
    def __make_graphCl_util__(self):
        self.gclutil = NestedDict()
    #
    def __init_graph_fromfile__(self, import_mode: str = SG_MPRT_MODE):
        match import_mode:
            case 'pkl'|'pk'|'pickle':
                return pk.load(open(self.graphPath+PKL, "rb"))
            case 'gml':
                return nx.read_gml(self.graphPath+GML)
    #
    def __init_reprdict__(self):
        self.GraphReprDict[self.onGraph] = self.G
        for grd in SG_LIST_REPR:
            try:
                self.GraphReprDict[grd] = getattr(self, grd)
            except AttributeError:
                pass
        self.GraphReprs = list(self.GraphReprDict.keys())
    #
    def __init_weights__(self, values: Union[float, List] = 1) -> None:
        nx.set_edge_attributes(self.gr[self.onGraph], values, 'weight')
        if type(values) is list:
            self.upd_GraphRepr_All(self.onGraph)
    #
    def __init_nNodesEdges__(self) -> None:
        self.N = self.gr[self.onGraph].number_of_nodes()
        self.Ne = self.gr[self.onGraph].number_of_edges()
        self.nodes_in[self.onGraph] = list(self.gr[self.onGraph].nodes())
    #
    def __init_sgraph__(self, init_weights_val: Union[float, List] = 1) -> None:
        self.__init_nNodesEdges__()
        on_g = self.onGraph
        if self.import_on:
            self.neflip = self.Ne_n
            self.pflip = self.neflip / self.Ne
            edgesWithData = self.gr[on_g].edges(data=True)
            self.fleset[on_g] = set([(u, v) for u, v, _ in edgesWithData 
                                         if _.get('weight', 1) < 0])
            self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
        else:
            self.neflip = int(self.pflip * self.Ne)
            self.nflip = int(self.pflip * self.N)
            edges = self.gr[on_g].edges
            self.eset[on_g] = set(list(edges()))
            self.fleset[on_g] = set(self.get_random_links(n=self.neflip, 
                                                          on_g=on_g))
            self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
            self.__init_weights__(init_weights_val)
        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices()
    #
    # export graph tools
    #
    def export_adj_bin(self, print_msg: bool = False) -> None:
        rowarr = [row[i:] for i, row in enumerate(self.Adj.todense())]
        exname = f"{self.expOut}adj_{self.stdFname}.bin"
        if print_msg:
            print(f"exporting {exname}\n")
        with open(exname, "wb") as f:
            for i in range(len(rowarr)):
                rowarr[i].astype("float64").tofile(f)
    #
    def export_edgel_bin(self, on_g: str = SG_GRAPH_REPR, exName: str = "",
                         print_msg: bool = False, mode: str = 'numpy') -> None:
        edges = self.GraphReprDict[on_g].edges(data='weight')
        exname = '_'.join([self.pEqStr, exName]) if exName else self.pEqStr
        fname = '_'.join(["edgelist", exname])+BIN
        edglFname = os.path.join(self.expOut, fname)
        match mode:
            case 'numpy':
                edge_array = np.array(list(edges), 
                                 dtype=[('i', np.uint64),
                                        ('j', np.uint64),
                                        ('w_ij', np.float64)])
                self.edglFile = open(edglFname, "wb")
                edge_array.tofile(self.edglFile)
            case 'struct':
                with open(edglFname, "wb") as f:
                    for edge in edges:
                        f.write(struct.pack("QQd", edge[0], edge[1], edge[2]))
    #
    def remove_edgl_file(self):
        os.remove(self.edglFile.name)
    
    def export_graphPKL(self):
        pk.dump(self.G, open(self.graphPath+PKL, "wb"), pk.HIGHEST_PROTOCOL)
    #
    def export_graphGML(self):
        nx.write_gml(self.G, self.graphPath+GML)
    #
    def __export_graph__(self, xprt_modeL: str = SG_XPRT_MODE):
        match xprt_modeL:
            case 'pk'|'pickle'|'pkl':
                self.export_graphPKL()
            case 'gml':
                self.export_graphGML()
    #
    # graph get attributes
    #
    def get_node_attributes(self, attr: str = 'pos', on_g: str = SG_GRAPH_REPR):
        return nx.get_node_attributes(self.gr[on_g], attr)
    #
    def get_edge_data(self, u: Any, v: Any, thedata: str = 'weight',
                      on_g: str = SG_GRAPH_REPR):
        return self.gr[on_g].get_edge_data(u, v)[thedata]
    #
    def get_edge_mapping_or_reverse(self, edge, target_g, on_g: str = SG_GRAPH_REPR):
        try:
            return self.edgeMap[target_g][on_g][edge]
        except KeyError:
            rev_edge = edge[::-1]
            return self.edgeMap[target_g][on_g].get(rev_edge, None)
    #
    def get_edge_color(self, pec: ColorType = "blue", nec: ColorType = "red",
                    thedata: str = 'weight', on_g: str = SG_GRAPH_REPR,
                    continuous: bool = False, cmap: str = 'coolwarm'):
        def map_values(value):
            if continuous:
                norm_value = (value - min_val) / (max_val - min_val)  # Normalize value to [0, 1]
                return plt.get_cmap(cmap)(norm_value)  # Get color from colormap
            return nec if value == -1 else pec if value == 1 else value
        
        arr = nx.get_edge_attributes(self.gr[on_g], thedata)

        if continuous:
            values = np.array(list(arr.values()))
            min_val, max_val = values.min(), values.max()

        return list(map(map_values, arr.values()))
    #  
    def get_graph_neighbors(self, node: Any, on_g: str = SG_GRAPH_REPR):
        return list(self.gr[on_g].neighbors(node))
    #
    def get_adjacency_matrix(self, on_g: str = SG_GRAPH_REPR, 
                             weight: str = 'weight', format: str = 'csr'):
        return nx.to_scipy_sparse_array(self.gr[on_g], 
                                        weight=weight, format=format)
    #
    def get_degree_matrix(self, format: str = 'csr'):
        return spdiags(self.adj.sum(axis=1), 0, *self.adj.shape, format=format)
    #
    def get_abs_degree_matrix(self, format: str = 'csr'):
        return spdiags(abs(self.adj).sum(axis=1), 0, *self.adj.shape, 
                       format=format)
    #
    def get_laplacian(self):
        return self.degm - self.adj
    #
    def get_signed_laplacian(self):
        return self.sdeg - self.adj
    #
    def get_random_links(self, only_in: str = '', n: int = 1,
                         on_g: str = SG_GRAPH_REPR):
        match only_in:
            case ''|'all':
                return random.sample(tuple(self.eset[on_g]), n)
            case '+'|'positive'|'+1'|'plus':
                return random.sample(tuple(self.lfeset[on_g]), n)
            case '-'|'negative'|'-1'|'minus':
                return random.sample(tuple(self.fleset[on_g]), n)
    #
    def get_eigV_binarized(self, which: int = 0):
        eigV = self.eigV[which].squeeze()
        return flip_to_positive_majority(bin_sign(eigV)).squeeze()
    #
    def get_eigV_bin_check(self, which: int = 0):
        if not hasattr(self, f"eigV") or which >= len(self.eigV):
            self.compute_k_eigvV(k=which+1)
        return self.get_eigV_binarized(which)
    #
    def get_bineigV_all(self):
        try:
            eigVbin = bin_sign(self.eigV)
        except (AttributeError, IndexError):
            self.compute_k_eigvV()
            eigVbin = bin_sign(self.eigV)
        for i in range(len(eigVbin)):
            eigVbin[i] = flip_to_positive_majority(eigVbin[i])
        return eigVbin
    #
    def get_signed_laplacian_embedding(self, k: int = 2):
        return self.eigV[:k]
    #
    def get_subgraph_from_nodes(self, list_of_nodes, on_g: str = SG_GRAPH_REPR):
        return self.gr[on_g].subgraph(list_of_nodes)
    #
    def get_nodes_subgraph_by_kv(self, k, val, on_g: str = SG_GRAPH_REPR):
        G = self.gr[on_g]
        G_yes, G_no = G.copy(), G.copy()
        for node, v in G.nodes(data=k):
            if v != val:
                G_yes.remove_node(node)
            else:
                G_no.remove_node(node)
        return G_yes, G_no
    #
    def get_bineigV_cluster_sizes(self, which: int = 0, 
                      binarize: bool = True, on_g: str = SG_GRAPH_REPR):
        if not all(f"eigV{which}" in self.gr[on_g].nodes[node] 
                   for node in self.gr[on_g].nodes):
            self.load_eigV_on_graph(which, on_g, binarize)
        if not hasattr(self, "clustersY"):
            self.make_clustersYN(f"eigV{which}", +1, on_g)
        cl_len = sorted(map(len, self.clustersY), reverse=True)
        return cl_len
    #
    def get_cluster_distribution(self, which: int = 0, 
                                  on_g: str = SG_GRAPH_REPR,
                                  binarize: bool = True):
        cl_len = self.get_bineigV_cluster_sizes(which, on_g, binarize)
        dictdist_cluster_sizes = {
            size: cl_len.count(size) for size in set(cl_len)
        }
        return dictdist_cluster_sizes
    #
    def get_ferroAntiferro_regions(self, attr_str: str = 's', 
                                        on_g: str = SG_GRAPH_REPR):
        antiGroup = []
        ferroGroup = []
        graph = self.gr[on_g]
        for node, att in graph.nodes(data=attr_str):
            neighbors = graph.neighbors(node)
            if all(graph.nodes[n][attr_str] == -att for n in neighbors):
                antiGroup.append(node)
            elif all(graph.nodes[n][attr_str] == att for n in neighbors):
                ferroGroup.append(node)
        return ferroGroup, antiGroup
    #
    # update graph methods
    #
    def upd_graph_matrices(self, format: str = 'csr', 
                           on_g: str = SG_GRAPH_REPR):
        self.Adj = self.get_adjacency_matrix(on_g=on_g)
        self.Deg = self.get_degree_matrix()
        self.sDeg = self.get_abs_degree_matrix()
        self.Lap = self.get_laplacian()
        self.sLp = self.get_signed_laplacian()
        self.upd_Nen(on_g)
        self.upd_Degree(on_g)
    #
    def upd_Degree(self, on_g: str = SG_GRAPH_REPR):
        self.degrees = list(dict(self.gr[on_g].degree).values())
    #
    def upd_Nen(self, on_g: str = SG_GRAPH_REPR):
        self.Ne_n = len(self.fleset[on_g])
    #
    def zip_reprNodes(self, x, on_g: str = SG_GRAPH_REPR):
        return dict(zip(self.gr[on_g], self.gr[x]))
    #
    def zip_reprEdges(self, x, on_g: str = SG_GRAPH_REPR):
        return dict(zip(self.gr[on_g].edges(), self.gr[x].edges()))
    #
    def upd_NodeMap(self, on_g: str = SG_GRAPH_REPR):
        self.nodeMap[on_g] = {x: 
            {v: k for k, v in self.zip_reprNodes(x, on_g).items()} 
            for x in self.GraphReprs if x != on_g}
    #
    def upd_EdgeMap(self, on_g: str = SG_GRAPH_REPR):
        graph = self.gr[on_g]
        egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
        self.edgeMap[on_g] = {x: {v: k for k, v in egraph_to(x).items()} 
                                   for x in self.GraphReprs if x != on_g}
    #
    def upd_GraphRelabel(self, on_g: str = SG_GRAPH_REPR, 
                         to_graph: str = SG_GRAPHINT_REPR):
        self.gr[to_graph] = nx.relabel_nodes(self.gr[on_g], 
                                             self.nodeMap[to_graph][on_g],
                                             copy=True)
    #
    def upd_ReprMaps(self, on_g: str = SG_GRAPH_REPR):
        graph = self.gr[on_g]
        ngraph_to = lambda x: dict(zip(graph, self.gr[x]))
        self.nodeMap[on_g] = {x: {v: k for k, v in ngraph_to(x).items()} 
                                   for x in self.GraphReprs if x != on_g}
        egraph_to =  lambda x: dict(zip(graph.edges(), self.gr[x].edges()))
        self.edgeMap[on_g] = {x: {v: k for k, v in egraph_to(x).items()} 
                                   for x in self.GraphReprs if x != on_g}
    #
    def upd_GraphRepr(self, on_g: str = SG_GRAPH_REPR):
        self.upd_ReprMaps(on_g)
    #
    def upd_GraphRepr_All(self, on_g: str = SG_GRAPH_REPR, 
                          also_itself: bool = True):
        if also_itself:
            self.upd_GraphRepr(on_g)
            self.upd_Nen(on_g)
        for i in self.GraphReprs:
            if i != on_g:
                self.upd_GraphRepr(i)
                self.upd_GraphRelabel(on_g, i)
                self.nodes_in[i] = list(self.gr[i].nodes())
                self.eset[i] = {x for e in self.eset[on_g] 
                                if (x := self.get_edge_mapping_or_reverse(e, i, on_g))}
                self.fleset[i] = {x for e in self.fleset[on_g]
                                  if (x := self.get_edge_mapping_or_reverse(e, i, on_g))}
                self.lfeset[i] = {x for e in self.lfeset[on_g]
                                  if (x := self.get_edge_mapping_or_reverse(e, i, on_g))}
                self.upd_Nen(i)
    #
    # graph operations
    #
    def check_pflip(self):
        if self.neflip < 1: raise NflipError(SG_ERRMSG_NFLIP)
    #
    def flip_sel_edges(self, links, on_g: str = SG_GRAPH_REPR):
        neg_weights_dict = {
            (u, v): -1 * self.get_edge_data(u, v, on_g=on_g) 
            for u, v in links
        }
        nx.set_edge_attributes(self.gr[on_g], values=neg_weights_dict, 
                               name='weight')
        self.fleset[on_g].update(links)
        self.lfeset[on_g].difference_update(links)
        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices(on_g)
    #
    def flip_random_fract_edges(self, pflip: float = None, 
                                on_g: str = SG_GRAPH_REPR):
        try:
            if pflip:
                self.pflip = pflip
                self.neflip = int(self.pflip * (self.Ne-self.Ne_n))
                self.check_pflip()
                links = self.get_random_links(self.neflip, on_g)
                self.flip_sel_edges(links, on_g)
            else:
                self.check_pflip()
                self.flip_sel_edges(self.fleset[on_g], on_g)
        except NflipError:
            pass
    #
    def unflip_all(self, on_g: str = SG_GRAPH_REPR):
        self.flip_sel_edges(1, on_g)
    #
    def make_edges_random_normal(self, mu: float = 1.0, sigma: float = 1.0, 
                                 on_g: str = SG_GRAPH_REPR):
        weights = {
            edge: random.normalvariate(mu, sigma) 
            for edge in self.gr[on_g].edges()
        }
        nx.set_edge_attributes(self.gr[on_g], values=weights, name='weight')
        self.fleset[on_g] = set([(u, v) 
                                 for u, v, _ in self.gr[on_g].edges(data=True)
                                 if _.get('weight', 1) < 0])
        self.lfeset[on_g] = self.eset[on_g].difference(self.fleset[on_g])
        self.upd_GraphRepr_All(on_g)
        self.upd_graph_matrices(on_g)
    #
    def load_vec_on_nodes(self, vec: NDArray, attr: str,
                          on_g: str = SG_GRAPH_REPR):
        vecNodeAttr = {nd: v for v, nd in zip(vec, self.gr[on_g].nodes)}
        nx.set_node_attributes(self.gr[on_g], values=vecNodeAttr, name=attr)
    #
    def load_eigV_on_graph(self, which: int = 0, on_g: str = SG_GRAPH_REPR, 
                           binarize: bool = False):
        if binarize: 
            eigV = self.get_eigV_bin_check(which=which)
        else:
            try:
                eigV = self.eigV[which]
            except (IndexError, AttributeError):
                self.compute_k_eigvV(k=which + 1)
                eigV = self.eigV[which]
        # print(eigV, self.gr[on_g].nodes)
        eigV_val_nd = {nd: v for v, nd in zip(eigV, self.gr[on_g].nodes)}
        nx.set_node_attributes(self.gr[on_g], eigV_val_nd, f"eigV{which}")
    #
    # computations
    #
    def compute_laplacian_spectrum(self, typf: type = np.float64):
        # NOT WORKING, NEEDS FIX
        self.eigv = np.linalg.eigvalsh(self.Lap.astype(typf).todense())
    #
    def compute_full_laplacian_spectrum(self, typf: type = np.float64):
        self.eigv, self.eigV = np.linalg.eigh(self.sLp.astype(typf).todense())
        self.make_eigV_transposed()
    #
    def compute_k_eigvV(self, k: int = 1, with_routine: str = "scipy",
        which: str = "SM", typf: type = np.float64):
        if with_routine == "numpy" or k > self.N//2:
            self.compute_full_laplacian_spectrum(typf)
        elif with_routine.startswith("scipy"):
            mode = with_routine.split("_")
            mode = mode[-1] if len(mode) > 1 else "caley"
            self.eigv, self.eigV = scsp_eigsh(
                self.sLp.astype(typf), k=k, which=which, mode=mode
            )
            self.make_eigV_transposed()
    #
    def compute_pinf(self, which: int = 0, on_g: str = SG_GRAPH_REPR):
        clustd = np.array(self.get_bineigV_cluster_sizes(which, on_g))
        mclust = clustd[0]
        self.Pinf = mclust / self.N
        self.Pinf_var = np.sum(clustd@clustd-mclust**2)/(np.sum(clustd)-mclust)
        if hasattr(self, "Pinf_dict"):
            self.Pinf_dict[which] = (self.Pinf, self.Pinf_var)
        else:
            self.Pinf_dict = {which: (self.Pinf, self.Pinf_var)}
    #
    def compute_rbim_energy_eigV(self, which: int = 0):
        spins = self.get_eigV_bin_check(which)
        edges = self.gr[SG_GRAPH_REPR].edges(data=True)
        ui, vi, w = zip(*[(u, v, data['weight']) for u, v, data in edges])
        ui, vi = np.array(ui, dtype=int), np.array(vi, dtype=int)
        w = np.array(w)
        if not hasattr(self, "energy_eigV_RBIM"):
            self.energy_eigV_RBIM = {}
        self.energy_eigV_RBIM[which] = -np.sum(w * spins[ui] * spins[vi])
    #
    def get_rbim_energy_eigV(self, which: int = 0):
        if not hasattr(self, "energy_eigV_RBIM"):
            self.compute_rbim_energy_eigV(which)
        elif which not in self.energy_eigV_RBIM.keys():
                self.compute_rbim_energy_eigV(which)
        return self.energy_eigV_RBIM[which]
    #
    # make methods
    #
    def make_eigV_transposed(self):
        self.eigV = self.eigV.T
    #
    def make_rescaled_signed_laplacian(self, MODE: str = "field"):
        if MODE == "field":
            self.resLp = self.sLp - self.eigv[0] * scsp_identity(self.N)
        elif MODE == "double":
            from scipy.linalg import eigvalsh

            self.resLp = self.sLp - np.array([self.eigv[0]])
            new_eigv0 = eigvalsh(
                self.resLp.astype(np.float64), subset_by_index=[0, 0]
            )
            self.resLp = self.resLp - new_eigv0 * np.identity(self.N)
    #
    def make_graphYN(self, k, val, on_g: str = SG_GRAPH_REPR):
        self.gclutil[k][val][on_g] = self.get_nodes_subgraph_by_kv(k, val, on_g)
    #
    def make_clustersYN(self, k, val, on_g: str = SG_GRAPH_REPR):
        try:
            graphY, graphN = self.gclutil[k][val][on_g]
        except:
            self.make_graphYN(k, val, on_g)
            graphY, graphN = self.gclutil[k][val][on_g]
        #
        self.clustersY = list(nx.connected_components(graphY))
        self.clustersN = list(nx.connected_components(graphN))
        #
        self.numClustersY = len(self.clustersY)
        self.numClustersN = len(self.clustersN)
        #
        lgcY = max(self.clustersY, key=len) if self.clustersY else []
        lgcN = max(self.clustersN, key=len) if self.clustersN else []
        self.biggestClSet = self.clustersY if len(lgcY) >= len(lgcN) \
            else self.clustersN
        self.biggestClSet.sort(key=len, reverse=True)
        self.numClustersBig = len(self.biggestClSet)
        self.gc = max(self.biggestClSet, key=len)
    #
    def make_connected_component_by_edge(self, edge_attr='weight', value=1, 
                                         on_g: str = SG_GRAPH_REPR):
        connected_components = nx.connected_components
        value_edges = [
            (u, v) for u, v, attrs in self.gr[on_g].edges(data=True)
            if attrs.get(edge_attr) == value
        ]
        if not value_edges:
            raise ValueError("No positive edges found in the graph.")
        G_pos = self.gr[on_g].edge_subgraph(value_edges).copy()
        all_connected_components = list(connected_components(G_pos))
        if not all_connected_components:
            raise ValueError("No connected components found with positive edges.")
        self.largest_cc = max(all_connected_components, key=len)
        self.largest_cc_subgraph = G_pos.subgraph(self.largest_cc).copy()