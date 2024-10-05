from .common import *
from .BinDynSys import BinDynSys
class IsingDynamics:
    magn = []
    ene = []
    magnc1 = []
    magn_array_save = []
    Ising_clusters = []
    k_B = 1

    def __init__(
        self,
        sg: SignedGraph,
        T: float = 0.,
        ic: str = "uniform",
        runlang: str = "py",
        NoClust: int = 1,
        nstepsIsing: int = 100,
        save_magnetization: bool = False,
        upd_mode: str = "asynchronous",
        out_suffix: str = "",
        rndStr: bool = False,
        id_string: str = "",
    ) -> None:
        self.runlang = runlang
        self.sg = sg
        self.N = self.sg.N
        self.T = T
        self.ic = ic
        self.nstepsIsing = nstepsIsing
        self.save_magnetization = save_magnetization
        self.NoClust = NoClust
        self.upd_mode = upd_mode
        self.randstring_OPT = rndStr
        randstr_tmp = randstring() if rndStr else ""
        self.id_string_isingdyn = join_non_empty('_', id_string, randstr_tmp)
        self.pflip_id = f"p={self.sg.pflip:.3g}"
        self.in_suffix = join_non_empty('_', self.pflip_id, 
                                        self.id_string_isingdyn)
        self.out_suffix = out_suffix
    #
    def bltzmnn_fact(self, E: float) -> float:
        return np.exp(-E / (self.k_B * self.T))
    #
    def flip_spin(self, nd: int) -> None:
        self.s[nd] *= -1
    #
    def neigh_ene(self, neigh: list) -> float:
        return np.sum(neigh) / len(neigh)
    #
    def neigh_wghtmagn(self, node: int) -> list:
        nd = dict(self.sg.G[node])
        return [w["weight"] * self.s[nn] for nn, w in nd.items()]
    #
    def metropolis(self, node):
        neigh = self.neigh_wghtmagn(node)
        neighene = self.neigh_ene(neigh)
        DeltaE = 2 * self.s[node] * neighene
        if DeltaE < 0:
            self.flip_spin(node)
        elif np.random.uniform() < self.bltzmnn_fact(DeltaE):
            self.flip_spin(node)
    #
    def calc_full_energy(self):
        return np.array([
                -self.s[node] * self.neigh_ene(self.neigh_wghtmagn(node))
                for node in range(self.sg.N)]).sum()
    #
    def init_ising_dynamics(self, custom: Any = None):
        match self.ic:
            case "uniform":
                self.s = np.random.choice([-1, 1], size=self.sg.N)
            case _ if self.ic.startswith("ground_state"):
                number = int(self.ic.split("_")[-1])
                bineigv = self.sg.get_eigV_bin_check(which=number)
                self.s = bineigv
            case "custom":
                self.s = custom
            case "delta":
                self.s = np.zeros(self.sg.N)
                self.s[np.random.randint(self.sg.N)] = 1
            case _:
                raise ValueError("Invalid initial condition.")
        if self.runlang.startswith("C") and not self.sg.import_on:
            self.export_s_init()
            self.CbaseName = f"IsingSimulator{self.runlang[-1]}"
    #
    def run(  #
        self,
        tqdm_on: bool = True,
        thrmSTEP: int = 2,
        eqSTEP: int = 10,
        freq: int = 10,
        T_ising: float = None,
        verbose: bool = False,
    ):
        # name C0 and C1 to be modified, in C0 -> CS that is a single eigenstate is studied with all of his subclusters
        # name C1 -> CM where one studies the largest component of all the clusters
        try:
            getattr(self, f"CbaseName")
        except AttributeError:
            self.init_ising_dynamics()
        if T_ising:
            self.T = T_ising
        if self.runlang.startswith("C"):
            run_id = self.id_string_isingdyn
            self.run_id = f"_{run_id}" if run_id else ''
            arglist = [f"{self.N}",
                       f"{self.T:.3g}",
                       f"{self.sg.pflip:.3g}",
                       f"{self.NoClust}",
                       f"{thrmSTEP:.3g}",
                       f"{eqSTEP}",
                       self.sg.datPath,
                       self.sg.syshapePth,
                       self.run_id,
                       self.out_suffix,
                       self.upd_mode,
                       f"{freq}"
                       ]
            self.cprogram = [pth_join(LRGSG_LIB_CBIN, self.CbaseName)] + arglist
            print(self.cprogram)
            stderrFname = join_non_empty('_', f"err{self.CbaseName}", f"{self.N}", 
                                    f"{self.id_string_isingdyn}", 
                                    f"{self.out_suffix}")+LOG
            stderrFname = os.path.join(LRGSG_LOG, stderrFname)
            self.stderr = open(stderrFname, 'w')
            if verbose:
                print('\rExecuting: ', ' '.join(self.cprogram), end='', 
                      flush=True)
            result = subprocess.run(self.cprogram, stderr=self.stderr, 
                                    stdout=subprocess.PIPE)
            output = result.stdout
            self.s = np.frombuffer(output, dtype=np.int8)
        else:
            metropolis_1step = np.vectorize(self.metropolis, excluded="self")
            if self.save_magnetization:
                def save_magn_array():
                    self.magn_array_save.append(self.s)
            else:
                def save_magn_array():
                    pass
            sample = list(
                range(self.sg.N)
            )  # rd.sample(list(self.system.H.nodes()), self.system.N)
            iterator = (
                tqdm(range(self.nstepsIsing))
                if tqdm_on
                else range(self.nstepsIsing)
            )
            self.ene = []
            for _ in iterator:
                self.magn.append(np.sum(self.s))
                self.ene.append(self.calc_full_energy())
                # for i in range(self.system.N):
                #     self.metropolis(sample[i])
                metropolis_1step(sample)
                save_magn_array()

    def find_ising_clusters(self, import_cl: bool = False):
        #can be easily reworked
        if import_cl:
            for i in range(self.NoClust):
                self.Ising_clusters.append(
                    np.fromfile(
                        f"{self.sg.isingpath}cl{i}_{self.sg.stdFname}.bin",
                        dtype=int,
                    )
                )
            self.numIsing_cl = len(self.Ising_clusters)
        if self.Ising_clusters:
            print("Ising Clusters already computed.")
            return
        #
        self.sg.compute_k_eigvV(howmany=self.NoClust)
        eigVbin = self.sg.get_bineigV_all()
        #
        self.Ising_clusters = []
        for j in range(self.NoClust):
            lnodes = list(self.sg.H.nodes())
            lnodes_tmp = lnodes[:]

            def recursive_search(seed, magn_i, clustertmp):
                neighs = get_kth_order_neighbours(self.sg.H, seed, 1)
                neighs = np.array(
                    [e for e in neighs if e not in set(clustertmp)]
                )
                if not neighs.size:
                    return
                samecluster = np.array(eigVbin[j][neighs] == magn_i)
                if not samecluster.any():
                    return
                neighs_samecluster = list(neighs[samecluster])
                clustertmp.extend(neighs_samecluster)
                for ss in neighs_samecluster:
                    recursive_search(ss, magn_i, clustertmp)

            allclusters = []
            for i in lnodes:
                if i not in lnodes_tmp:
                    continue
                if not lnodes_tmp:
                    break
                #
                clustertmp = []
                clustertmp.extend([i])
                #
                recursive_search(i, eigVbin[j][i], clustertmp)
                lnodes_tmp = [e for e in lnodes_tmp if e not in set(clustertmp)]
                allclusters.append(clustertmp)
            allclusters.sort(key=len, reverse=True)
            self.Ising_clusters.append(allclusters[0])
        self.numIsing_cl = len(self.Ising_clusters)
        if self.runlang.startswith("C"):
            self.export_ising_clust()

    #
    def mapping_nodes_to_clusters(self):
        if not self.Ising_clusters:
            self.find_ising_clusters()
        loc = [x for x in range(len(self.Ising_clusters))]
        self.loc = loc
        node_with_inherclust = [
            [[j, loc[i]] for j in clus]
            for i, clus in enumerate(self.Ising_clusters)
        ]
        self.node_with_inherclust = node_with_inherclust
        node_inherclust_flat = [i for j in node_with_inherclust for i in j]
        self.node_inherclust_flat = node_inherclust_flat
        sorted_list = sorted(node_inherclust_flat, key=lambda x: x[0])
        self.sorted_list = sorted_list
        result_array = np.empty((self.sg.side1, self.sg.side2), dtype=object)
        self.result_array = result_array

        # Fill the result_array with tuples from sorted_list
        for i, sublist in enumerate(sorted_list):
            row, col = divmod(
                i, self.sg.side1
            )  # Calculate the row and column index
            result_array[row, col] = sublist[1]
        self.mapping = result_array

    #
    def export_s_init(self):
        fname = os.path.join(self.sg.isingpath,
                     '_'.join(['s', self.in_suffix])+BIN)
        self.sfout = open(fname, "wb")
        self.s.astype("int8").tofile(self.sfout)
        self.s_0 = self.s.copy()
    #
    def export_ising_clust(self, val: Any = +1, which: int = 0, 
                           on_g: str = SG_GRAPH_REPR, binarize: bool = True): 
        self.sg.load_eigV_on_graph(which, on_g, binarize)
        self.sg.make_clustersYN(f"eigV{which}", val, on_g=on_g)
        try:
            if self.NoClust > self.sg.numClustersBig:
                raise NoClustError(
                    "Requested number of Cluster files is bigger than the one"
                    " in selected topology."
                )
        except NoClustError as excpt:
            print(excpt)
        for i in range(self.NoClust):
            fname = os.path.join(self.sg.isingpath,
                                 join_non_empty('_', f"cl{i}", self.in_suffix)+BIN)
            self.clfout = open(fname, "wb")
            np.array(list(self.sg.biggestClSet[i])).astype(int).tofile(self.clfout)

    def remove_run_c_files(self, remove_stderr: bool = True):
        try:
            os.remove(self.sfout.name)
        except FileNotFoundError:
            pass
        if self.runlang == "C1":
            try:
                os.remove(self.clfout.name)
            except FileNotFoundError:
                pass
        if remove_stderr:
            try:
                os.remove(self.stderr.name)
            except FileNotFoundError:
                pass



























class IsingDynamics_DEV(BinDynSys):
    dyn_UVclass = "ising_model"
    id_string_isingdyn = ""

    magn = []
    ene = []
    magnc1 = []
    magn_array_save = []
    Ising_clusters = []
    k_B = 1

    def __init__(
        self,
        sg: SignedGraph = Lattice2D,
        T: float = 1.0,
        runlang: str = "py",
        NoClust: int = 1,
        nstepsIsing: int = 100,
        save_magnetization: bool = False,
        **kwargs
    ) -> None:
        self.T = T
        self.sg = sg
        self.dynpath = self.isingpath
        super(IsingDynamics_DEV, self).__init__(self.sg, **kwargs)
        self.nstepsIsing = nstepsIsing
        self.save_magnetization = save_magnetization
        self.NoClust = NoClust

    def bltzmnn_fact(self, enrgy: float) -> float:
        """
        Calculate the Boltzmann factor for a given energy.

        The Boltzmann factor is given by:
        exp(-enrgy / (k_B * T))

        Parameters:
        -----------
        - enrgy (float): Energy value for which the Boltzmann factor is calculated.

        Returns:
        -----------
        - float: Boltzmann factor.
        """
        return np.exp(-enrgy / (self.k_B * self.T))

    #
    def flip_spin(self, nd: int) -> None:
        """
        Flips the spin of a particle at a specified index within the system.

        Parameters:
        -----------
        - nd (int): The index of the particle/spin to be flipped.

        Returns:
        -----------
        - None

        Modifies the state of the system by flipping the sign of the particle
        at the given index `nd`.
        """
        self.s[nd] *= -1
    #
    def neigh_ene(self, neigh: list) -> float:
        """
        Calculate the average energy of a list of neighboring energy values.

        Parameters:
        -----------
        - neigh (list): A list of neighboring energy values.

        Returns:
        -----------
        - float: The average energy.
        """
        return np.sum(neigh) / len(neigh)
    #
    def neigh_wghtmagn(self, node: int) -> list:
        """
        Calculate the weighted magnitudes of neighboring nodes' values for a given node.

        Parameters:
        -----------
        - node (int): The node for which to calculate the weighted magnitudes of neighbors.

        Returns:
        -----------
        - list: A list of weighted magnitudes.
        """
        node_dict = dict(self.sg.H[node])
        return [w["weight"] * self.s[nn] for nn, w in node_dict.items()]

    # #
    # def metropolis(self, node):
    #     neigh = self.neigh_wghtmagn(node)
    #     neighene = self.neigh_ene(neigh)
    #     E_old = -self.s[node] * neighene
    #     E_new = +self.s[node] * neighene
    #     DeltaE = E_new - E_old
    #     if DeltaE < 0:
    #         self.flip_spin(node)
    #     elif np.random.uniform() < self.bltzmnn_fact(DeltaE):
    #         self.flip_spin(node)

    # #
    # def calc_full_energy(self):
    #     return np.array(
    #         [
    #             -self.s[node] * self.neigh_ene(self.neigh_wghtmagn(node))
    #             for node in range(self.sg.N)
    #         ]
    #     ).sum()

    # #
    # def init_ising_dynamics(self, randstring_OPT: bool = True):
    #     self.id_string_isingdyn = randstring() if randstring_OPT else ""
    #     if self.ic == "uniform":
    #         self.s = np.random.choice([-1, 1], size=self.sg.N)
    #     elif self.ic.startswith("ground_state"):
    #         number = int(self.ic.split("_")[-1])
    #         bineigv = self.sg.bin_eigV(which=number)
    #         self.s = bineigv
    #     if self.runlang.startswith("C") and not self.sg.import_on:
    #         self.export_s_init()

    # #
    # def run(  #
    #     self,
    #     adjfname: str = "",
    #     out_suffix: str = "",
    #     tqdm_on: bool = True,
    #     thrmSTEP: int = 2,
    #     eqSTEP: int = 10,
    # ):
    #     # name C0 and C1 to be modified, in C0 -> CS that is a single eigenstate is studied with all of his subclusters
    #     # name C1 -> CM where one studies the largest component of all the clusters
    #     if self.runlang.startswith("C"):
    #         if adjfname == "":
    #             adjfname = self.sg.stdFname
    #         out_suffix = out_suffix + self.id_string_isingdyn
    #         if out_suffix == "":
    #             out_suffix = '""'
    #         self.cprogram = [
    #             f"src/lrgsglib/IsingSimulator{self.runlang[-1]}",
    #             f"{self.sg.N}",
    #             f"{self.T}",
    #             f"{self.sg.pflip}",
    #             adjfname,
    #             f"{self.NoClust}",
    #             f"{thrmSTEP}",
    #             f"{eqSTEP}",
    #             self.sg.datPath,
    #             out_suffix,
    #         ]
    #         call(self.cprogram)
    #     else:
    #         metropolis_1step = np.vectorize(self.metropolis, excluded="self")
    #         if self.save_magnetization:

    #             def save_magn_array():
    #                 self.magn_array_save.append(self.s)

    #         else:

    #             def save_magn_array():
    #                 pass

    #         sample = list(
    #             range(self.sg.N)
    #         )  # rd.sample(list(self.system.H.nodes()), self.system.N)
    #         iterator = (
    #             tqdm(range(self.nstepsIsing))
    #             if tqdm_on
    #             else range(self.nstepsIsing)
    #         )
    #         self.ene = []
    #         for _ in iterator:
    #             self.magn.append(np.sum(self.s))
    #             self.ene.append(self.calc_full_energy())
    #             # for i in range(self.system.N):
    #             #     self.metropolis(sample[i])
    #             metropolis_1step(sample)
    #             save_magn_array()

    # def find_ising_clusters(self, import_cl: bool = False):
    #     if import_cl:
    #         for i in range(self.NoClust):
    #             self.Ising_clusters.append(
    #                 np.fromfile(
    #                     f"{self.sg.isingpath}cl{i}_{self.sg.stdFname}.bin",
    #                     dtype=int,
    #                 )
    #             )
    #         self.numIsing_cl = len(self.Ising_clusters)
    #     if self.Ising_clusters:
    #         print("Ising Clusters already computed.")
    #         return
    #     #
    #     self.sg.compute_k_eigvV(howmany=self.NoClust)
    #     eigVbin = self.sg.bin_eigV_all()
    #     #
    #     self.Ising_clusters = []
    #     for j in range(self.NoClust):
    #         lnodes = list(self.sg.H.nodes())
    #         lnodes_tmp = lnodes[:]

    #         def recursive_search(seed, magn_i, clustertmp):
    #             neighs = get_kth_order_neighbours(self.sg.H, seed, 1)
    #             neighs = np.array(
    #                 [e for e in neighs if e not in set(clustertmp)]
    #             )
    #             if not neighs.size:
    #                 return
    #             samecluster = np.array(eigVbin[j][neighs] == magn_i)
    #             if not samecluster.any():
    #                 return
    #             neighs_samecluster = list(neighs[samecluster])
    #             clustertmp.extend(neighs_samecluster)
    #             for ss in neighs_samecluster:
    #                 recursive_search(ss, magn_i, clustertmp)

    #         allclusters = []
    #         for i in lnodes:
    #             if i not in lnodes_tmp:
    #                 continue
    #             if not lnodes_tmp:
    #                 break
    #             #
    #             clustertmp = []
    #             clustertmp.extend([i])
    #             #
    #             recursive_search(i, eigVbin[j][i], clustertmp)
    #             lnodes_tmp = [e for e in lnodes_tmp if e not in set(clustertmp)]
    #             allclusters.append(clustertmp)
    #         allclusters.sort(key=len, reverse=True)
    #         self.Ising_clusters.append(allclusters[0])
    #     self.numIsing_cl = len(self.Ising_clusters)
    #     if self.runlang.startswith("C"):
    #         self.export_ising_clust()

    # #
    # def mapping_nodes_to_clusters(self):
    #     if not self.Ising_clusters:
    #         self.find_ising_clusters()
    #     loc = [x for x in range(len(self.Ising_clusters))]
    #     self.loc = loc
    #     node_with_inherclust = [
    #         [[j, loc[i]] for j in clus]
    #         for i, clus in enumerate(self.Ising_clusters)
    #     ]
    #     self.node_with_inherclust = node_with_inherclust
    #     node_inherclust_flat = [i for j in node_with_inherclust for i in j]
    #     self.node_inherclust_flat = node_inherclust_flat
    #     sorted_list = sorted(node_inherclust_flat, key=lambda x: x[0])
    #     self.sorted_list = sorted_list
    #     result_array = np.empty((self.sg.side1, self.sg.side2), dtype=object)
    #     self.result_array = result_array

    #     # Fill the result_array with tuples from sorted_list
    #     for i, sublist in enumerate(sorted_list):
    #         row, col = divmod(
    #             i, self.sg.side1
    #         )  # Calculate the row and column index
    #         result_array[row, col] = sublist[1]
    #     self.mapping = result_array

    # #
    # def export_s_init(self):
    #     output_file = open(
    #         f"{self.sg.isingpath}s_{self.sg.stdFname}.bin",
    #         "wb",
    #     )
    #     self.s.astype("int8").tofile(output_file)

    # #
    # def export_ising_clust(self):
    #     try:
    #         if self.NoClust > self.numIsing_cl:
    #             raise NoClustError(
    #                 "Requested number of Cluster files is bigger than the one"
    #                 " in selected topology."
    #             )
    #     except NoClustError as excpt:
    #         print(excpt)

    #     for i in range(self.NoClust):
    #         output_file = open(
    #             f"{self.sg.isingpath}cl{i}_{self.sg.stdFname}.bin",
    #             "wb",
    #         )
    #         np.array(self.Ising_clusters[i]).astype(int).tofile(output_file)











