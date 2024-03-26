from .nx_objects import *

class Lattice2D(SignedGraph):
    #
    def __init__(
        self,
        side1: int = L2D_SIDE1,
        side2: int = L2D_SIDE2,
        geo: str = L2D_GEO,
        pbc: bool = L2D_PBC,
        fbc_val: float = L2D_FBCV,
        stdFnameSFFX: str = L2D_STDFN,
        sgpath: str = L2D_SGPATH,
        with_positions: bool = False,
        **kwargs,
    ) -> None:
        self.__init_side__(side1, side2)
        self.pbc = pbc
        self.__init_geo__(geo)
        self.fbc_val = fbc_val
        #
        #
        self.__init_stdFname__(stdFnameSFFX)
        #
        _ = L2D_PATH_DICT[self.geo]
        self.sgpath = sgpath + _ if sgpath else _
        self.with_positions = with_positions
        #
        self.__init_lattice__()
        super(Lattice2D, self).__init__(self.G, **kwargs)
    #
    def __init_side__(self, side1: int, side2: int) -> None:
        # to do add check for variables to be int and so on...
        if side2:
            if side2 > side1:
                self.side2 = side1
                self.side1 = side2
            else:
                self.side1 = side1
                self.side2 = side2
        else:
            self.side1 = side1
        #
    #
    def __init_geo__(self, geo: str) -> None:
        self.geo = geo
        if geo not in L2D_GEO_LIST:
            if geo not in L2D_GEO_SHRT_LIST:
                warnings.warn(L2D_WARNMSG_GEO, Lattice2DWarning)
                self.geo = L2D_GEO
            else:
                self.geo = L2D_SHRT_GEO_DICT[geo]
        if not hasattr(self, 'side2'):
            if self.geo == 'hexagonal':
                self.side2 = self.side1
                self.side1 = adjust_to_even(self.side1/np.sqrt(3))
                if (self.side1 % 2 or self.side2 % 2) and self.pbc:
                    raise ValueError(L2D_ERRMSG_GEO)
            else:
                self.side2 = self.side1
            # if (self.side1 % 2 or self.side2 % 2) and self.pbc:
            #     raise ValueError(DEFLattice2D_geoerrmsg)
        
    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = L2D_GEO_SHRT_DICT[self.geo] + SFFX
    #
    def __init_lattice__(self) -> None:
        #
        if self.geo == L2D_SHRT_GEO_DICT['tri']:
            nxfunc = triangular_lattice_graph_FastPatch
            self.syshape = (self.side1, self.side2)
        elif self.geo == L2D_SHRT_GEO_DICT['sqr']:
            nxfunc = squared_lattice_graph_FastPatch
            self.syshape = (self.side1, self.side2)
        elif self.geo == L2D_SHRT_GEO_DICT['hex']:
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
        self.p_c = L2D_P_C_DICT[self.geo]
        self.r_c = np.sqrt(1.128/(np.pi*self.p_c))
        #
        self.G = nxfunc(self.side1, self.side2, periodic=self.pbc, 
                        with_positions=self.with_positions)
    #
    def degree_check(self, degree):
        return np.where(np.array(list(map(lambda x: x[1], 
                                          list(self.G.degree())))) != degree)
    #
    def get_central_edge(self, on_graph: str = L2D_ONREP):
        G = self.G
        cnode = (self.side1//2-1, self.side2//2)
        cnode_t = (self.side1//2, self.side2//2)
        if self.geo == 'triangular':
            cnode = (self.side2//2, self.side1//2-1)
            cnode_t = (self.side2//2, self.side1//2)
        edge_t = (cnode, cnode_t)
        if not G.has_edge(*edge_t):
            if self.geo =='hexagonal':
                cnode = cnode_t
                cnode_t = (self.side1//2+1, self.side2//2)
                edge_t = (cnode, cnode_t)
        if on_graph == 'G':
            return edge_t
        elif on_graph == 'H':
            return self.edge_map[edge_t]
    #
    class nwContainer(dict):
        def __init__(self, l: SignedGraph, iterable=[], constant=None, 
                     **kwargs):
            super().__init__(**kwargs)
            self.update((key, constant) for key in iterable)
            self.l = l
            self.reprDict = list(self.l.GraphReprDict.keys())
            self.rNodeFlip = {g: random.sample(
                                    list(self.l.GraphReprDict[g].nodes()), 
                                    int(self.l.pflip*self.l.N)
                                ) for g in self.reprDict}
            #
            self.centedge = {g: self.l.get_central_edge(g) 
                             for g in self.reprDict}
            self['single'] = {g: [self.centedge[g]] for g in self.reprDict}
            self['singleZERR'] = {g: self.get_links_ZERR(
                self.centedge[g][0], g, self.l.geo) for g in self.reprDict}
            self['singleXERR'] = {g: self.get_links_XERR(
                self.centedge[g][0], g) for g in self.reprDict}
            self['rand'] = {g: [e for e in self.l.rEdgeFlip[g]] 
                            for g in self.reprDict}
            self['randZERR'] = {g: self.get_rand_pattern('ZERR', on_graph=g) 
                                   for g in self.reprDict}
            self['randXERR'] = {g: self.get_rand_pattern('XERR', on_graph=g) 
                                   for g in self.reprDict}
        #
        def get_links_XERR(self, node: Any, on_graph: str = L2D_ONREP):
            return [(node, nn) for nn in self.l.graph_neighbors(node, on_graph)]
        #
        def get_links_ZERR(self, node: Any, on_graph: str = L2D_ONREP, 
                           geometry: str = L2D_GEO):
            dd = {'triangular': self.get_links_triangle,
             'squared': self.get_links_square,
             'hexagonal': self.get_links_hexagon}
            return dd[geometry](node, on_graph)
        #
        def get_links_triangle(self, node: Any, on_graph = L2D_ONREP):
            node2 = list(self.l.graph_neighbors(node, on_graph))[0]
            common_neighbors = list(nx.common_neighbors(
                self.l.GraphReprDict[on_graph], node, node2))
            try:
                node3 = common_neighbors[0]
                links = [(node, node2), (node2, node3), (node, node3)]
            except IndexError:
                links = [(node, node2)]
            return links
        #
        def get_links_square(self, node: Any, on_graph = L2D_ONREP):
            g = self.l.GraphReprDict[on_graph]
            neighbors = list(g.neighbors(node))
            for i in range(1, len(neighbors)):
                first_neighbor = neighbors[0]  # Always the first neighbor
                second_neighbor = neighbors[i]  # Iterating through the rest
                # Find common neighbors excluding the start_node
                common_neighbors = set(g.neighbors(first_neighbor)) \
                    & set(g.neighbors(second_neighbor))
                common_neighbors.discard(node)
                # If there is a common neighbor, we found a square
                if common_neighbors:
                    common_neighbor = common_neighbors.pop()
                    # Now, extract the links forming the square
                    links = [(node, first_neighbor),
                            (node, second_neighbor),
                            (first_neighbor, common_neighbor),
                            (second_neighbor, common_neighbor)]
                    return links
            links = [(node, first_neighbor),
                     (node, second_neighbor)]
            return links
        #
        def get_links_hexagon(self, node: int, 
                              on_graph: str = L2D_ONREP):
            graph = self.l.GraphReprDict[on_graph]
            nodes_in_cycle = [node]
            node_nn = list(self.l.GraphReprDict[on_graph].neighbors(node))

            samp_node_nn_1 = node_nn[0]
            node_nn.remove(samp_node_nn_1)
            # #
            node_nn_1 = list(graph.neighbors(samp_node_nn_1))
            node_nn_1.remove(node)
            samp_node_nn_2 = node_nn_1[0]
            node_nn_1.remove(samp_node_nn_2)
            #
            flag = True
            for nn in node_nn:
                if flag:
                    common_nn = None
                    node_nn_1b = list(graph.neighbors(nn))
                    node_nn_1b.remove(node)
                    for i in node_nn_1b:
                        common_neighs = list(
                            nx.common_neighbors(graph, samp_node_nn_2, i))
                        if common_neighs != []:
                            nodes_in_cycle.extend([nn, i, common_neighs[0]])
                            flag = False
                            break
            try:
                nodes_in_cycle.extend([samp_node_nn_1, samp_node_nn_2])
            except IndexError:
                pass
            subH = graph.subgraph(nodes_in_cycle)
            links = [tuple(sorted(edge)) for edge in subH.edges()]
            return links
        #
        def get_rand_pattern(self, mode: str, 
                             on_graph: str = L2D_ONREP):
            if mode == "ZERR":
                if self.l.geo == 'squared':
                    mode = "square"
                elif self.l.geo == 'triangular':
                    mode = "triangle"
                elif self.l.geo == 'hexagonal':
                    mode = "hexagon"
            if mode == "XERR":
                patternList = [k for i in self.rNodeFlip[on_graph] 
                               for k in self.get_links_XERR(i, on_graph)]
            elif mode == "hexagon":
                patternList = [k for i in self.rNodeFlip[on_graph]
                                for k in self.get_links_hexagon(i, on_graph)]
            elif mode == "square":
                patternList = [k for i in self.rNodeFlip[on_graph] 
                               for k in self.get_links_square(i, on_graph)]
            elif mode == "triangle":
                patternList = [k for i in self.rNodeFlip[on_graph] 
                               for k in self.get_links_triangle(i, on_graph)]
            return list(set(patternList))
        #
        def get_links_rball(self, R: int = 1, center: Any = None, 
                            on_graph: str = L2D_ONREP):
            graph = self.l.GraphReprDict[on_graph]
            if not center:
                center = self.centedge[on_graph][0]
            neighs_to_flip = get_neighbors_within_distance(graph, center, R)
            links = {(node, neighbor) for node in neighs_to_flip 
                     for neighbor in graph.neighbors(node)}
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