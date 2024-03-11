from .nx_objects import *

class Lattice2D(SignedGraph):
    #
    def __init__(
        self,
        side1: int = DEFLattice2D_side1,
        side2: int = DEFLattice2D_side2,
        geo: str = DEFLattice2D_geo,
        pbc: bool = DEFLattice2D_pbc,
        fbc_val: float = DEFLattice2D_fbcv,
        stdFnameSFFX: str = DEFLattice2D_stdFn,
        sgpath: str = DEFLattice2D_sgpath,
        with_positions: bool = False,
        **kwargs,
    ) -> None:
        #
        self.side1 = side1
        self.side2 = side2 if side2 else side1
        #
        self.__init_geo__(geo)
        #
        self.pbc = pbc
        self.fbc_val = fbc_val
        #
        self.__init_stdFname__(stdFnameSFFX)
        #
        _ = DEFLattice2D_pthdict[self.geo]
        self.sgpath = sgpath + _ if sgpath else _
        self.with_positions = with_positions
        #
        self.__init_lattice__()
        super(Lattice2D, self).__init__(self.G, **kwargs)
    #
    def __init_geo__(self, geo: str):
        self.geo = geo
        if geo not in DEFLattice2D_geolist:
            if geo not in DEFLattice2D_geoabblist:
                warnings.warn(DEFLattice2D_geowarnmsg, Lattice2DWarning)
                self.geo = DEFLattice2D_geo
            else:
                self.geo = DEFLattice2D_geodictabb[geo]
        if geo == 'hexagonal':
            if self.side1 % 2 or self.side2 % 2:
                raise ValueError(DEFLattice2D_geoerrmsg)
    #
    def __init_stdFname__(self, SFFX: str = ""):
        self.stdFname = DEFLattice2D_geodict[self.geo] + SFFX
    #
    def __init_lattice__(self, pbc = None) -> None:
        if pbc is None:
            pbc = self.pbc
        else:
            pbc = False
        #
        if self.geo == DEFLattice2D_geodictabb['tri']:
            nxfunc = triangular_lattice_graph_FastPatch
            self.syshape = (self.side1, self.side2)
        elif self.geo == DEFLattice2D_geodictabb['sqr']:
            nxfunc = squared_lattice_graph_FastPatch
            self.syshape = (self.side1, self.side2)
        elif self.geo == DEFLattice2D_geodictabb['hex']:
            nxfunc = hexagonal_lattice_graph_FastPatch
        #
        if self.side1 == self.side2:
            self.syshapeStr = f"N={self.side1**2}"
        elif self.side1 > self.side2:
            self.syshapeStr = f"L1={self.side1}_L2={self.side2}"
        elif self.side2 > self.side1:
            self.syshapeStr = f"L1={self.side2}_L2={self.side1}"
        self.syshapePth = f"{self.syshapeStr}/"
        #
        self.p_c = DEFLattice2D_p_cdict[self.geo]
        self.r_c = np.sqrt(1.128/(np.pi*self.p_c))
        #
        self.G = nxfunc(self.side1, self.side2, periodic=pbc, with_positions=self.with_positions)
    #
    def degree_check(self, degree):
        return np.where(np.array(list(map(lambda x: x[1], list(self.G.degree())))) != degree)
    #
    def get_central_edge(self, on_graph: str = 'G'):
        G = self.G
        if self.geo == 'triangular':
            cnode = (self.side1//2, self.side2//2)
            cnode_t = (self.side1//2+1, self.side2//2)
        else:
            cnode = (self.side1//2-1, self.side2//2)
            cnode_t = (self.side1//2, self.side2//2)
        edge_t = (cnode, cnode_t)
        if self.geo == 'squared':
            return edge_t
        if not G.has_edge(cnode, cnode_t):
            if self.geo =='hexagonal':
                cnode = cnode_t
                cnode_t = (self.side1//2+1, self.side2//2)
            edge_t = (cnode, cnode_t)
        if on_graph == 'G':
            return edge_t
        elif on_graph == 'H':
            return self.edge_map[edge_t]
    #
    class neg_weights_dicts_container(dict):
        NEG_WEIGHTS_DICT_G_PFLIP = {}
        def __init__(self, l: SignedGraph, iterable=[], constant=None, **kwargs):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.l = l
            self.rNodeFlip = {g: random.sample(
                list(self.l.GraphReprDict[g].nodes()), int(self.l.pflip*self.l.N))
                for g in ['G', 'H']}
            for n in self.rNodeFlip['G']:
                if n not in set(self.l.G.nodes()):
                    print(n, ' node error')
            # #
            self.centedge = {g: self.l.get_central_edge(g) for g in ['G', 'H']}
            self['rand'] = {g: [e for e in self.l.rEdgeFlip[g]] for g in ['G', 'H']}
            self['randXERR'] = {g: self.get_nwd_pattern('cross', on_graph=g) 
                                   for g in ['G', 'H']}
            self['randZERR'] = {g: self.get_nwd_pattern('unit_cell', on_graph=g) 
                                   for g in ['G', 'H']}
            self['central'] = {g: [self.centedge[g]]
                                   for g in ['G', 'H']}
        #
        def get_links_cross(self, node: Any, on_graph: str = 'G'):
            return [(node, nn) for nn in self.l.graph_neighbors(node, on_graph)]
        #
        def get_links_triangle(self, node: Any, on_graph = 'G'):
            node2 = list(self.l.graph_neighbors(node, on_graph))[0]
            common_neighbors = list(nx.common_neighbors(self.l.GraphReprDict[on_graph], node, node2))
            node3 = common_neighbors[0]
            links = [(node, node2), (node2, node3), (node, node3)]
            return links
        #
        def get_links_square(self, node: Any, on_graph = 'G'):
            graph = self.l.GraphReprDict[on_graph]
            neighbors = list(graph.neighbors(node))    
            for i in range(1, len(neighbors)):
                first_neighbor = neighbors[0]  # Always the first neighbor
                second_neighbor = neighbors[i]  # Iterating through the rest
                # Find common neighbors excluding the start_node
                common_neighbors = set(graph.neighbors(first_neighbor)) & set(graph.neighbors(second_neighbor))
                common_neighbors.discard(node)
                # If there is a common neighbor, we found a square
                if common_neighbors:
                    common_neighbor = common_neighbors.pop()  # Get one common neighbor to form a square
                    # Now, extract the links forming the square
                    links = [(node, first_neighbor),
                            (node, second_neighbor),
                            (first_neighbor, common_neighbor),
                            (second_neighbor, common_neighbor)]
                    
                    return links
            return None
        #
        def get_links_hexagon(self, node: int, on_graph: str = 'G'):
            graph = self.l.GraphReprDict[on_graph]
            nodes_in_cycle = [node]
            node_nn = list(self.l.GraphReprDict[on_graph].neighbors(node))

            samp_node_nn_1 = random.sample(node_nn, 1)[0]
            node_nn.remove(samp_node_nn_1)
            # #
            node_nn_1 = list(graph.neighbors(samp_node_nn_1))
            node_nn_1.remove(node)
            samp_node_nn_2 = random.sample(node_nn_1, 1)[0]
            node_nn_1.remove(samp_node_nn_2)
            #
            flag = True
            for nn in node_nn:
                if flag:
                    common_nn = None
                    node_nn_1b = list(graph.neighbors(nn))
                    node_nn_1b.remove(node)
                    for i in node_nn_1b:
                        common_neighs = list(nx.common_neighbors(graph, samp_node_nn_2, i))
                        if common_neighs != []:
                            common_nn = i
                            nodes_in_cycle.append(nn)
                            flag = False
                            break
            nodes_in_cycle.extend([samp_node_nn_1, samp_node_nn_2, common_nn, common_neighs[0]])
            subH = graph.subgraph(nodes_in_cycle)
            links = [tuple(sorted(edge)) for edge in subH.edges()]
            return links
        #
        def get_nwd_pattern(self, mode: str, on_graph: str = 'G'):
            if mode == "unit_cell":
                if self.l.geo == 'squared':
                    mode = "square"
                elif self.l.geo == 'triangular':
                    mode = "triangle"
                elif self.l.geo == 'hexagonal':
                    mode = "hexagon"
            if mode == "cross":
                patternList = [k for i in self.rNodeFlip[on_graph] 
                               for k in self.get_links_cross(i, on_graph)]
            elif mode == "hexagon":
                patternList = [k for i in self.rNodeFlip[on_graph]
                                for k in self.get_links_hexagon(i, on_graph)]
            elif mode == "square":
                patternList = [k for i in self.rNodeFlip[on_graph] 
                               for k in self.get_links_square(i, on_graph)]
            elif mode == "triangle":
                patternList = [k for i in self.rNodeFlip[on_graph] 
                               for k in self.get_links_triangle(i, on_graph)]
            return patternList
        #
        def get_links_rball(self, R: int = 1, center: Any = None, on_graph: str = 'G'):
            if not center:
                center = self.centedge[on_graph][0]
            neighs_to_flip = get_neighbors_within_distance(self.l.GraphReprDict[on_graph], center, R)
            links = {(node, neighbor) for node in neighs_to_flip for neighbor in self.l.GraphReprDict[on_graph].neighbors(node)}
            return links
        # #
    def make_animation(self, fig, ax, frames):

        cv0 = frames[0].reshape(self.syshape)
        im = ax.imshow(cv0)  # Here make an AxesImage rather than contour
        _, _, cbar = imshow_colorbar_caxdivider(im, ax)
        # tx = ax.set_title('Frame 0')

        def animate(i):
            arr = frames[i].reshape(self.syshape)
            vmax = np.max(arr)
            vmin = np.min(arr)
            im.set_data(arr)
            cbar.mappable.set_clim(vmin, vmax)
            # tx.set_text('Frame {0}'.format(i))
            # In this version you don't have to do anything to the colorbar,
            # it updates itself when the mappable it watches (im) changes

        return animate