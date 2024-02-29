from .nx_objects import *

class SignedGraph:
    p_c = None
    lsp = None
    slspectrum = None
    pflmin = DEFAULT_MIN_PFLIPVAL
    pflmax = DEFAULT_MAX_PFLIPVAL

    def __init__(
        self,
        G: Graph,
        import_on: bool = False,
        pflip: float = 0.0,
        lsp_mode: str = "intervals",
        expathc: str = "",
        init_weight_dict: bool = True
    ):
        self.__init_paths__()
        if not self.pflmin <= pflip <= self.pflmax:
            raise ValueError(f"pflip must be between {self.pflmin} and {self.pflmax}, inclusive. Received: {pflip}")
        else:
            self.pflip = pflip
        self.lsp_mode = lsp_mode
        self.import_on = import_on
        self.expath = (
            expathc if expathc else f"{self.DEFAULT_GRAPHDIR}{self.syshapePTH}"
        )
        self.isingpath = f"{self.DEFAULT_ISINGDIR}{self.syshapePTH}"
        self.voterpath = f"{self.DEFAULT_VOTERDIR}{self.syshapePTH}"
        self.lrgsgpath = f"{self.DEFAULT_LRGSGDIR}{self.syshapePTH}"
        self.phtrapath = f"{self.DEFAULT_PHTRADIR}{self.syshapePTH}"
        self.__make_dirs__()
        self.stdFname = self.stdFname + f"_p={self.pflip:.3g}"
        if import_on:
            self.graphfname = self.expath + self.stdFname
            self.G = self.__init_graph_fromfile__()
        else:
            self.G = G
            self.graphfname = self.expath + self.stdFname
        self.init_sgraph()
        if init_weight_dict:
            self.neg_weights_dict = self.neg_weights_dicts_container(self)

    #
    def __init_graph_fromfile__(self):
        return pickle.load(open(f"{self.graphfname}.pickle", "rb"))

    #
    def __init_paths__(self):
        self.datPath = f"{DEFAULT_DATA_OUTDIR}{self.sgpath}"
        self.pltPath = (
            f"{DEFAULT_DATA_OUTDIR}{DEFAULT_PLOT_OUTDIR}{self.sgpath}"
        )
        self.DEFAULT_GRAPHDIR = self.datPath + DEFAULT_GRAPH_OUTDIR
        self.DEFAULT_ISINGDIR = self.datPath + DEFAULT_ISING_OUTDIR
        self.DEFAULT_VOTERDIR = self.datPath + DEFAULT_VOTER_OUTDIR
        self.DEFAULT_LRGSGDIR = self.datPath + DEFAULT_LRGSG_OUTDIR
        self.DEFAULT_PHTRADIR = self.datPath + DEFAULT_PHTRA_OUTDIR

    #
    def __make_dirs__(self):
        for _ in [self.expath, self.isingpath, self.voterpath, self.lrgsgpath, 
                  self.phtrapath]:
            os.makedirs(_, exist_ok=True)

    #
    def init_weights(self):
        nx.set_edge_attributes(self.G, values=1, name="weight")

    #
    def number_of_negative_links(self):
        self.Ne_n = (
            np.array(list(nx.get_edge_attributes(self.H, "weight").values()))
            < 0
        ).sum()

    #
    def upd_G_graph(self):
        self.invnode_map = {v: k for k, v in self.node_map.items()}
        self.invedge_map = {v: k for k, v in self.edge_map.items()}
        self.G = nx.relabel_nodes(self.H, self.invnode_map)
        self.esetG = list(self.G.edges())
        self.number_of_negative_links()

    #
    def upd_H_graph(self):
        self.H = nx.convert_node_labels_to_integers(self.G)
        self.esetH = list(self.H.edges())
        self.node_map = dict(zip(self.G, self.H))
        self.edge_map = dict(zip(self.G.edges(), self.H.edges()))
        self.number_of_negative_links()

    def init_n_nodes_edges(self):
        self.N = self.G.number_of_nodes()
        self.Ne = self.G.number_of_edges()

    def init_sgraph(self):
        self.init_n_nodes_edges()
        if self.import_on:
            self.upd_H_graph()
            self.nflip = self.Ne_n
            self.pflip = self.nflip / self.Ne
            self.upd_graph_matrices()
            self.randsample = np.where(
                np.array(self.Adj.todense()).flatten() < 0
            )
        else:
            self.nflip = int(self.pflip * self.Ne)
            self.randsample = random.sample(range(self.Ne), self.nflip)
            self.init_weights()
            self.upd_H_graph()
        self.upd_G_graph()

    #
    def adjacency_matrix(self, weight: str = "weight"):
        return nx.to_scipy_sparse_array(self.H, weight=weight, format="csr")

    #
    def degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(scsp.spdiags(A.sum(axis=1), 0, *A.shape, format="csr"))

    #
    def absolute_degree_matrix(self, A: csr_array) -> csr_array:
        return csr_array(
            scsp.spdiags(abs(A).sum(axis=1), 0, *A.shape, format="csr")
        )

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
    def upd_graph_matrices(self, on_graph="H"):
        if on_graph == "G":
            motherNx = self.G
        elif on_graph == "H":
            motherNx = self.H
        self.Adj = self.adjacency_matrix()
        self.Deg = self.degree_matrix(self.Adj)
        self.sDeg = self.absolute_degree_matrix(self.Adj)
        self.Lap = self.laplacian_matrix()
        self.sLp = self.signed_laplacian()

    #
    def flip_sel_edges(self, neg_weights_dict=None, on_graph="H"):
        """Flips a specific edges of a graph G."""
        #
        if on_graph == "G":
            if neg_weights_dict is None:
                neg_weights_dict = self.neg_weights_dict.DEFAULT_NEG_WEIGHTS_DICT_G
            nx.set_edge_attributes(
                self.G, values=neg_weights_dict, name="weight"
            )
        elif on_graph == "H":
            if neg_weights_dict is None:
                neg_weights_dict = self.neg_weights_dict.DEFAULT_NEG_WEIGHTS_DICT_H
            nx.set_edge_attributes(
                self.H, values=neg_weights_dict, name="weight"
            )
        self.upd_graph(on_graph=on_graph)
        self.upd_graph_matrices()

    #
    def check_pflip(self):
        if self.nflip < 1:
            raise NflipError(
                """The probability of flipping an edge times the 
                             number of edges is < 1, then no edges would be
                             flipped. Skipping the analysis for this value."""
            )
    #
    def upd_graph(self, on_graph: str = "G"):
        if on_graph == "G":
            self.upd_H_graph()
        elif on_graph == "H":
            self.upd_G_graph()
    #
    def unflip_all(self, on_graph="H"):
        if on_graph == "G":
            eset = self.esetG
        elif on_graph == "H":
            eset = self.esetH
        self.flip_sel_edges(neg_weights_dict={e: 1 for e in eset}, on_graph=on_graph)

    #
    def flip_random_fract_edges(self, on_graph="H"):
        """Flips a fraction p of edges (+1 to -1) of a graph G."""
        #
        try:
            self.check_pflip()
        except NflipError:
            return None
        if on_graph == "G":
            eset = self.esetG
        elif on_graph == "H":
            eset = self.esetH
        self.flip_sel_edges(neg_weights_dict={e: -1 for e in random.sample(eset, self.nflip)}, on_graph=on_graph)

    #
    def compute_laplacian_spectrum(self, MODE_lapspec: str = "numpy") -> None:
        if MODE_lapspec == "networkx":
            self.slspectrum = nx.laplacian_spectrum(self.system.G)
        elif MODE_lapspec == "numpy":
            self.slspectrum = np.linalg.eigvalsh(self.sLp.toarray())
    #
    def compute_k_eigvV(self, MODE_dynspec: str = "scipy", howmany: int = 1, which: str = "SM"):
        if MODE_dynspec == "numpy" or howmany == self.N:
            self.eigv, self.eigV = np.linalg.eigh(self.sLp.astype(np.float64).todense())
            self.eigV = self.eigV.T
        if MODE_dynspec == "scipy":
            self.eigv, self.eigV = scsp.linalg.eigsh(
                self.sLp.astype(np.float64), k=howmany, which=which
            )
            self.eigV = self.eigV.T

    #
    def bin_eigV(self, which: int = 0):
        try:
            eigVbin = np.sign(np.where(self.eigV[which] == 0, +1, self.eigV[which]))
        except (AttributeError, IndexError):
            self.compute_k_eigvV(howmany = which + 1)
            eigVbin = np.sign(np.where(self.eigV[which] == 0, +1, self.eigV[which]))

        return eigVbin

    def bin_eigV_all(self):
        try:
            eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
        except (AttributeError, IndexError):
            self.compute_k_eigvV()
            eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
        return eigVbin
    
    def calc_fluct_Pinf(self, which: int = 0):
        eigV = self.bin_eigV(which)
        eV_p = np.count_nonzero(eigV >= 0)
        eV_n = self.N - eV_p

        self.eigV_fluct = abs(eV_p - eV_n) / self.N
        self.Pinf = np.min([eV_p, eV_n]) / self.N

        if hasattr(self, 'eigV_fluct_dict'):
            self.eigV_fluct_dict[which] = self.eigV_fluct
        else:
            self.eigV_fluct_dict = {which: self.Pinf}

        if hasattr(self, 'Pinf_dict'):
            self.Pinf_dict[which] = self.Pinf
        else:
            self.Pinf_dict = {which: self.Pinf}


    #
    def rescaled_signed_laplacian(self, MODE: str = "field"):
        if MODE == "field":
            self.resLp = self.sLp - self.eigv[0] * scsp.identity(self.N)
        elif MODE == "double":
            self.resLp = self.sLp - np.array([self.eigv[0]])
            new_eigv0 = scipy.linalg.eigvalsh(
                self.resLp.astype(np.float64), subset_by_index=[0, 0]
            )
            self.resLp = self.resLp - new_eigv0 * np.identity(self.N)

    #
    def lsp_selection(self, custom_list):
        if self.lsp_mode == "custom":
            self.lsp = np.array(custom_list)
        elif self.lsp_mode == "intervals":
            intervals = []
            tmp = max([vset["rsf"] for vset in custom_list])
            for vset in custom_list:
                if vset["kind"] == "log":
                    spacing_f = np.logspace
                    vset["start"] = np.log10(vset["start"])
                    vset["stop"] = np.log10(vset["stop"])
                elif vset["kind"] == "lin":
                    spacing_f = np.linspace
                intervals.append(  #
                    round_sigfig_n(  #
                        spacing_f(
                            vset["start"],
                            vset["stop"],
                            num=vset["num"],
                            endpoint=False,
                        ),
                        vset["rsf"],
                    )
                )
            self.lsp = (intervals := np.concatenate(intervals))
            while set(self.lsp).__len__() == intervals.__len__():
                tmp = tmp - 1
                self.lsp = np.round(self.lsp, tmp)
            tmp = tmp + 1
            self.lsp = np.round(intervals, tmp)

    def default_dict_lsp(self, num_low=3, num_at=6, num_high=3):
        d = (  #
            {
                "kind": "lin",
                "start": 0.001,
                "stop": self.p_c - self.p_c * num_at / 100,
                "num": num_low,
                "rsf": 1,
            },
            {
                "kind": "lin",
                "start": self.p_c - self.p_c * num_at / 100,
                "stop": self.p_c + self.p_c * num_at / 100,
                "num": num_at,
                "rsf": 3,
            },
            {
                "kind": "lin",
                "start": self.p_c + self.p_c * num_at / 100,
                "stop": 1,
                "num": num_high,
                "rsf": 1,
            },
        )
        return d
    #
    def dfs_list(self, node, visited, sign):
        if visited[node] or sign[node] <= 0:
            return 0
        visited[node] = True
        size = 1
        for neighbor in self.H[node]:
            if not visited[neighbor]:
                size += self.dfs_list(neighbor, visited, sign)
        return size
    #
    def cluster_distribution_list(self, sv = None):
        visited = [False] * len(self.H)
        distribution = {}
        if sv is None:
            try:
                sv = self.eigV[0]
            except AttributeError:
                self.compute_k_eigvV()
                sv = self.eigV[0]
        for node in range(len(self.H)):
            if not visited[node] and sv[node] > 0:
                size = self.dfs_list(node, visited, sv)
                distribution[size] = distribution.get(size, 0) + 1
        return distribution
    #
    def export_graph(self, MODE: str = "pickle"):
        if MODE == "pickle":
            pickle.dump(
                self.G,
                open(f"{self.graphfname}.pickle", "wb"),
                pickle.HIGHEST_PROTOCOL,
            )
        elif MODE == "gml":
            nx.write_gml(self.G, f"{self.graphfname}.gml")

    #
    def export_adj_bin(self, print_msg: bool = False) -> None:
        rowarr = [row[i:] for i, row in enumerate(self.Adj.todense())]
        exname = f"{self.expath}adj_{self.stdFname}.bin"
        if print_msg:
            print(f"exporting {exname}\n")
        with open(exname, "wb") as f:
            for i in range(len(rowarr)):
                rowarr[i].astype("float64").tofile(f)

    #
    def export_edgel(self):
        # TO BE FIXED
        a = list(self.H.edges(data="weight"))
        with open(r"src/LRGSG_package/tmp_stuff/prova.txt", "w") as fp:
            for item in a:
                # write each item on a new line
                fp.write("%s %s %s\n" % item)
            print("Done")


"""





























"""


class SignedGraph_DEV(Graph):
    p_c = None
    lsp = None
    sgpath = ""
    syshapepth = ""
    stdFname = ""
    slspectrum = None
    pflmin = DEFAULT_MIN_PFLIPVAL
    pflmax = DEFAULT_MAX_PFLIPVAL

    def __init__(
        self,
        import_on: bool = False,
        pflip: float = 0.0,
        expathc: str = "",
    ):
        super().__init__()
        self.import_on = import_on
        self.__init_paths__(expathc)
        self.__make_dirs__()
        if not self.pflmin <= pflip <= self.pflmax:
            raise ValueError(f"pflip must be between {self.pflmin} and {self.pflmax}, inclusive. Received: {pflip}")
        else:
            self.pflip = pflip
        self.stdFname = self.stdFname + f"_p={self.pflip:.3g}"
        self.graphfname = self.graphpath + self.stdFname
        if import_on:
            self.__init_graph_fromfile__()
        self.__init_weights__()
    #
    def __init_paths__(self, expathc: str = ""):
        self.datPath = f"{DEFAULT_DATA_OUTDIR}{self.sgpath}"
        self.DEFAULT_GRAPHDIR = self.datPath + DEFAULT_GRAPH_OUTDIR
        self.DEFAULT_ISINGDIR = self.datPath + DEFAULT_ISING_OUTDIR
        self.DEFAULT_VOTERDIR = self.datPath + DEFAULT_VOTER_OUTDIR
        self.DEFAULT_LRGSGDIR = self.datPath + DEFAULT_LRGSG_OUTDIR
        self.DEFAULT_PHTRADIR = self.datPath + DEFAULT_PHTRA_OUTDIR
        self.graphpath = f"{self.DEFAULT_GRAPHDIR}{self.syshapepth}{expathc}"
        self.isingpath = f"{self.DEFAULT_ISINGDIR}{self.syshapepth}{expathc}"
        self.voterpath = f"{self.DEFAULT_VOTERDIR}{self.syshapepth}{expathc}"
        self.lrgsgpath = f"{self.DEFAULT_LRGSGDIR}{self.syshapepth}{expathc}"
        self.phtrapath = f"{self.DEFAULT_PHTRADIR}{self.syshapepth}{expathc}"
    #
    def __make_dirs__(self):
        for _ in [self.graphpath, self.isingpath, self.voterpath, self.lrgsgpath, 
                  self.phtrapath]:
            os.makedirs(_, exist_ok=True)
    #
    def __init_graph_fromfile__(self):
        try:
            with open(f"{self.graphfname}.pickle", "rb") as file:
                loaded_obj = pickle.load(file)
                for attr_name in loaded_obj.__dict__:
                    setattr(self, attr_name, getattr(loaded_obj, attr_name))
        except FileNotFoundError:
            print(f"Error: Pickle file '{self.graphfname}.pickle' not found.")


    #
    def __init_weights__(self):
        nx.set_edge_attributes(self, values=1, name="weight")

    #
    def nlinks_count(self):
        weights = list(nx.get_edge_attributes(self, "weight").values())
        return (np.array(weights) < 0).sum()
    #
    def assign_attributes_from_nx_graph(self, nx_graph):
        for attr_name in nx_graph.graph:
            setattr(self, attr_name, nx_graph.graph[attr_name])
    #
    def relabel_nodes_to_integers(self):
        self.mapping_to_integers = {node: idx for idx, node in enumerate(self.nodes())}
        self.graph = nx.relabel_nodes(self.graph, self.mapping_to_integers)

    # #
    # def upd_H_graph(self):
    #     self.H = nx.convert_node_labels_to_integers(self.G)
    #     self.esetH = list(self.H.edges())
    #     self.node_map = dict(zip(self.G, self.H))
    #     self.edge_map = dict(zip(self.G.edges(), self.H.edges()))
    #     self.number_of_negative_links()

    # def init_n_nodes_edges(self):
    #     self.N = self.G.number_of_nodes()
    #     self.Ne = self.G.number_of_edges()

    # def __init_sgraph__(self):
    #     self.init_n_nodes_edges()
    #     if self.import_on:
    #         self.upd_H_graph()
    #         self.nflip = self.Ne_n
    #         self.pflip = self.nflip / self.Ne
    #         self.upd_graph_matrices()
    #         self.randsample = np.where(
    #             np.array(self.Adj.todense()).flatten() < 0
    #         )
    #     else:
    #         self.nflip = int(self.pflip * self.Ne)
    #         self.randsample = random.sample(range(self.Ne), self.nflip)
    #         self.init_weights()
    #         self.upd_H_graph()
    #     self.upd_G_graph()

    # #
    # def adjacency_matrix(self, weight: str = "weight"):
    #     return nx.to_scipy_sparse_array(self.H, weight=weight, format="csr")

    # #
    # def degree_matrix(self, A: csr_array) -> csr_array:
    #     return csr_array(scsp.spdiags(A.sum(axis=1), 0, *A.shape, format="csr"))

    # #
    # def absolute_degree_matrix(self, A: csr_array) -> csr_array:
    #     return csr_array(
    #         scsp.spdiags(abs(A).sum(axis=1), 0, *A.shape, format="csr")
    #     )

    # #
    # def laplacian_matrix(self) -> csr_array:
    #     """Returns the signed Laplacian matrix of G.
    #     The graph Laplacian is the matrix L = D - A, where
    #     A is the adjacency matrix and D is the diagonal matrix of node degrees

    #     Returns
    #     -------
    #     L : SciPy sparse array
    #     The Laplacian matrix of G.
    #     """
    #     return self.Deg - self.Adj

    # #
    # def signed_laplacian(self) -> csr_array:
    #     """Returns the signed Laplacian matrix of G.
    #     The graph Laplacian is the matrix L = |D| - A, where
    #     A is the adjacency matrix and |D| is the diagonal matrix of absolute
    #     values of node degrees

    #     Returns
    #     -------
    #     L : SciPy sparse array
    #     The Laplacian matrix of G.
    #     """
    #     return self.sDeg - self.Adj

    # #
    # def upd_graph_matrices(self, on_graph="H"):
    #     if on_graph == "G":
    #         motherNx = self.G
    #     elif on_graph == "H":
    #         motherNx = self.H
    #     self.Adj = self.adjacency_matrix()
    #     self.Deg = self.degree_matrix(self.Adj)
    #     self.sDeg = self.absolute_degree_matrix(self.Adj)
    #     self.Lap = self.laplacian_matrix()
    #     self.sLp = self.signed_laplacian()

    # #
    # def flip_sel_edges(self, neg_weights_dict=None, on_graph="H"):
    #     """Flips a specific edges of a graph G."""
    #     #
    #     if on_graph == "G":
    #         if neg_weights_dict is None:
    #             neg_weights_dict = self.neg_weights_dict.DEFAULT_NEG_WEIGHTS_DICT_G
    #         nx.set_edge_attributes(
    #             self.G, values=neg_weights_dict, name="weight"
    #         )
    #     elif on_graph == "H":
    #         if neg_weights_dict is None:
    #             neg_weights_dict = self.neg_weights_dict.DEFAULT_NEG_WEIGHTS_DICT_H
    #         nx.set_edge_attributes(
    #             self.H, values=neg_weights_dict, name="weight"
    #         )
    #     self.upd_graph(on_graph=on_graph)
    #     self.upd_graph_matrices()

    # #
    # def check_pflip(self):
    #     if self.nflip < 1:
    #         raise NflipError(
    #             """The probability of flipping an edge times the 
    #                          number of edges is < 1, then no edges would be
    #                          flipped. Skipping the analysis for this value."""
    #         )
    # #
    # def upd_graph(self, on_graph: str = "G"):
    #     if on_graph == "G":
    #         self.upd_H_graph()
    #     elif on_graph == "H":
    #         self.upd_G_graph()
    # #
    # def unflip_all(self, on_graph="H"):
    #     if on_graph == "G":
    #         eset = self.esetG
    #     elif on_graph == "H":
    #         eset = self.esetH
    #     self.WEIGHTS_DICT_H_PFLIP = {e: 1 for i, e in enumerate(eset)}
    #     self.flip_sel_edges(neg_weights_dict=self.WEIGHTS_DICT_H_PFLIP, on_graph=on_graph)

    # #
    # def flip_random_fract_edges(self, on_graph="H"):
    #     """Flips a fraction p of edges (+1 to -1) of a graph G."""

    #     #
    #     try:
    #         self.check_pflip()
    #     except NflipError:
    #         return None
    #     if on_graph == "G":
    #         eset = self.esetG
    #     elif on_graph == "H":
    #         eset = self.esetH
    #     self.NEG_WEIGHTS_DICT_H_PFLIP = {
    #         e: -1 for i, e in enumerate(eset) if i in self.randsample
    #     }
    #     self.flip_sel_edges(neg_weights_dict=self.NEG_WEIGHTS_DICT_H_PFLIP, on_graph=on_graph)

    # #
    # def compute_laplacian_spectrum(self, MODE_lapspec: str = "numpy") -> None:
    #     if MODE_lapspec == "networkx":
    #         self.slspectrum = nx.laplacian_spectrum(self.system.G)
    #     elif MODE_lapspec == "numpy":
    #         self.slspectrum = np.linalg.eigvalsh(self.sLp.toarray())
    # #
    # def compute_k_eigvV(self, MODE_dynspec: str = "scipy", howmany: int = 1, which: str = "SM"):
    #     if MODE_dynspec == "numpy" or howmany == self.N:
    #         self.eigv, self.eigV = np.linalg.eigh(self.sLp.astype(np.float64).todense())
    #         self.eigV = self.eigV.T
    #     if MODE_dynspec == "scipy":
    #         self.eigv, self.eigV = scsp.linalg.eigsh(
    #             self.sLp.astype(np.float64), k=howmany, which=which
    #         )
    #         self.eigV = self.eigV.T

    # #
    # def bin_eigV(self, which: int = 0):
    #     try:
    #         eigVbin = np.sign(np.where(self.eigV[which] == 0, +1, self.eigV[which]))
    #     except (AttributeError, IndexError):
    #         self.compute_k_eigvV(howmany = which + 1)
    #         eigVbin = np.sign(np.where(self.eigV[which] == 0, +1, self.eigV[which]))

    #     return eigVbin

    # def bin_eigV_all(self):
    #     try:
    #         eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
    #     except (AttributeError, IndexError):
    #         self.compute_k_eigvV()
    #         eigVbin = np.sign(np.where(self.eigV == 0, +1, self.eigV))
    #     return eigVbin
    
    # def calc_fluct_Pinf(self, which: int = 0):
    #     eigV = self.bin_eigV(which)
    #     eV_p = np.count_nonzero(eigV >= 0)
    #     eV_n = self.N - eV_p

    #     self.eigV_fluct = abs(eV_p - eV_n) / self.N
    #     self.Pinf = np.min([eV_p, eV_n]) / self.N

    #     if hasattr(self, 'eigV_fluct_dict'):
    #         self.eigV_fluct_dict[which] = self.eigV_fluct
    #     else:
    #         self.eigV_fluct_dict = {which: self.Pinf}

    #     if hasattr(self, 'Pinf_dict'):
    #         self.Pinf_dict[which] = self.Pinf
    #     else:
    #         self.Pinf_dict = {which: self.Pinf}


    # #
    # def rescaled_signed_laplacian(self, MODE: str = "field"):
    #     if MODE == "field":
    #         self.resLp = self.sLp - self.eigv[0] * scsp.identity(self.N)
    #     elif MODE == "double":
    #         self.resLp = self.sLp - np.array([self.eigv[0]])
    #         new_eigv0 = scipy.linalg.eigvalsh(
    #             self.resLp.astype(np.float64), subset_by_index=[0, 0]
    #         )
    #         self.resLp = self.resLp - new_eigv0 * np.identity(self.N)

    # #
    # def lsp_selection(self, custom_list):
    #     if self.lsp_mode == "custom":
    #         self.lsp = np.array(custom_list)
    #     elif self.lsp_mode == "intervals":
    #         intervals = []
    #         tmp = max([vset["rsf"] for vset in custom_list])
    #         for vset in custom_list:
    #             if vset["kind"] == "log":
    #                 spacing_f = np.logspace
    #                 vset["start"] = np.log10(vset["start"])
    #                 vset["stop"] = np.log10(vset["stop"])
    #             elif vset["kind"] == "lin":
    #                 spacing_f = np.linspace
    #             intervals.append(  #
    #                 round_sigfig_n(  #
    #                     spacing_f(
    #                         vset["start"],
    #                         vset["stop"],
    #                         num=vset["num"],
    #                         endpoint=False,
    #                     ),
    #                     vset["rsf"],
    #                 )
    #             )
    #         self.lsp = (intervals := np.concatenate(intervals))
    #         while set(self.lsp).__len__() == intervals.__len__():
    #             tmp = tmp - 1
    #             self.lsp = np.round(self.lsp, tmp)
    #         tmp = tmp + 1
    #         self.lsp = np.round(intervals, tmp)

    # def default_dict_lsp(self, num_low=3, num_at=6, num_high=3):
    #     d = (  #
    #         {
    #             "kind": "lin",
    #             "start": 0.001,
    #             "stop": self.p_c - self.p_c * num_at / 100,
    #             "num": num_low,
    #             "rsf": 1,
    #         },
    #         {
    #             "kind": "lin",
    #             "start": self.p_c - self.p_c * num_at / 100,
    #             "stop": self.p_c + self.p_c * num_at / 100,
    #             "num": num_at,
    #             "rsf": 3,
    #         },
    #         {
    #             "kind": "lin",
    #             "start": self.p_c + self.p_c * num_at / 100,
    #             "stop": 1,
    #             "num": num_high,
    #             "rsf": 1,
    #         },
    #     )
    #     return d
    # #
    # def dfs_list(self, node, visited, sign):
    #     if visited[node] or sign[node] <= 0:
    #         return 0
    #     visited[node] = True
    #     size = 1
    #     for neighbor in self.H[node]:
    #         if not visited[neighbor]:
    #             size += self.dfs_list(neighbor, visited, sign)
    #     return size
    # #
    # def cluster_distribution_list(self, sv = None):
    #     visited = [False] * len(self.H)
    #     distribution = {}
    #     if sv is None:
    #         try:
    #             sv = self.eigV[0]
    #         except AttributeError:
    #             self.compute_k_eigvV()
    #             sv = self.eigV[0]
    #     for node in range(len(self.H)):
    #         if not visited[node] and sv[node] > 0:
    #             size = self.dfs_list(node, visited, sv)
    #             distribution[size] = distribution.get(size, 0) + 1
    #     return distribution
    # #
    # def export_graph(self, MODE: str = "pickle"):
    #     if MODE == "pickle":
    #         pickle.dump(
    #             self.G,
    #             open(f"{self.graphfname}.pickle", "wb"),
    #             pickle.HIGHEST_PROTOCOL,
    #         )
    #     elif MODE == "gml":
    #         nx.write_gml(self.G, f"{self.graphfname}.gml")

    # #
    # def export_adj_bin(self, print_msg: bool = False) -> None:
    #     rowarr = [row[i:] for i, row in enumerate(self.Adj.todense())]
    #     exname = f"{self.expath}adj_{self.stdFname}.bin"
    #     if print_msg:
    #         print(f"exporting {exname}\n")
    #     with open(exname, "wb") as f:
    #         for i in range(len(rowarr)):
    #             rowarr[i].astype("float64").tofile(f)

    # #
    # def export_edgel(self):
    #     # TO BE FIXED
    #     a = list(self.H.edges(data="weight"))
    #     with open(r"src/LRGSG_package/tmp_stuff/prova.txt", "w") as fp:
    #         for item in a:
    #             # write each item on a new line
    #             fp.write("%s %s %s\n" % item)
    #         print("Done")
