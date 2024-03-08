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
            nxfunc = nx.grid_2d_graph
            self.syshape = (self.side1, self.side2)
        elif self.geo == DEFLattice2D_geodictabb['hex']:
            nxfunc = nx.hexagonal_lattice_graph
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
        self.G = nxfunc(self.side1, self.side2, periodic=pbc, 
                        **DEFLattice2D_kwfuncdict[self.geo])
        if self.geo == DEFLattice2D_geo:
            nx.set_node_attributes(self.G, values=dict(zip(self.G, self.G)), 
                                   name="pos")

    #
    class neg_weights_dicts_container:
        NEG_WEIGHTS_DICT_H_PFLIP = {}
        NEG_WEIGHTS_DICT_G_PFLIP = {}
        def __init__(self, lattice: SignedGraph):
            self.lattice = lattice
            self.midway_e = lattice.Ne//2
            self.midway_H = lattice.N//2 + lattice.side2//2
            self.H_cent_edge = self.get_central_edge_H()
            self.NEG_WEIGHTS_DICT_H_PFLIP = {e: -1 for e in random.sample(self.lattice.esetH, self.lattice.nflip)}
            self.WEIGHTS_DICT_H_PFLIP = {e: 1 for e in self.lattice.esetH}
            #
            # self.NEG_WEIGHTS_DICT_H_2ADJ = {lattice.esetH[self.midway_e]: -1, 
            #                                 lattice.esetH[self.midway_e+4]: -1}
            self.NEG_WEIGHTS_DICT_H_2CONT = {lattice.esetH[self.midway_e]: -1, 
                                            lattice.esetH[self.midway_e+2]: -1}
            #
            self.NEG_WEIGHTS_DICT_H_CROSS = self.get_neg_weights_dict_h_cross(self.H_cent_edge[0])
            self.NEG_WEIGHTS_DICT_H_SQUARE = self.get_neg_weights_dict_h_square(self.H_cent_edge[0])
            #
            self.node_flip_selection = np.random.choice(lattice.N, int(lattice.pflip * lattice.N))
            self.NEG_WEIGHTS_DICT_H_PCROSS = self.get_neg_weights_dict_h_pattern('cross')
            if lattice.geo == 'squared':
                try: 
                    self.DEFAULT_NEG_WEIGHTS_DICT_H = {self.H_cent_edge: -1}
                    self.DEFAULT_NEG_WEIGHTS_DICT_G = {lattice.invedge_map[self.H_cent_edge]: -1}
                except KeyError:
                    self.H_cent_edge = (self.H_cent_edge[0]-1, self.H_cent_edge[1]-1)
                    self.DEFAULT_NEG_WEIGHTS_DICT_H = {self.H_cent_edge: -1}
                    self.DEFAULT_NEG_WEIGHTS_DICT_G = {lattice.invedge_map[self.H_cent_edge]: -1}
                    pass
                self.NEG_WEIGHTS_DICT_H_PSQUARE = self.get_neg_weights_dict_h_pattern('square')
            elif lattice.geo == 'triangular':
                self.NEG_WEIGHTS_DICT_H_PTRIA = self.get_neg_weights_dict_h_pattern('triangle')
                self.G_cent_edge = ((lattice.side1//2-1,lattice.side2//2), (lattice.side1//2, lattice.side2//2))
                self.DEFAULT_NEG_WEIGHTS_DICT_G = {self.G_cent_edge: -1}
                self.DEFAULT_NEG_WEIGHTS_DICT_H = {lattice.edge_map[self.G_cent_edge]: -1}
            elif lattice.geo == 'hexagonal':
                self.NEG_WEIGHTS_DICT_H_PHEXA = self.get_neg_weights_dict_h_pattern('hexagon')
                

        #
        def get_central_edge_H(self):
            if self.lattice.geo == 'squared':
                return self.lattice.esetH[self.midway_e +1]
            elif self.lattice.geo == 'triangular':
                return (self.midway_H-1, self.midway_H)
            elif self.lattice.geo == 'hexagonal':
                target_node = (self.lattice.side2//2, self.lattice.side1)
                edge_with_change = [edge for edge in self.lattice.esetG if target_node in edge and
                                    any(other_node[0] == target_node[0] + 1 or other_node[0] == target_node[0] - 1
                                        for other_node in edge if other_node != target_node)][0]
                return self.lattice.edge_map[edge_with_change]
        #
        def get_neg_weights_dict_h_right(self, node: int):
            edge = (node, node+1)
            if edge not in self.lattice.esetH:
                edge = (node-1, node)
            dictH = {
                edge: -1,
            }
            return dictH
        #
        def get_neg_weights_dict_h_square(self, node: int):
            dictH = {
                (node, node+1): -1,
                (node, node+self.lattice.side1): -1,
                (node + 1, node+self.lattice.side1 + 1): -1,
                (node + self.lattice.side1, node+self.lattice.side1 + 1): -1,
            }
            return dictH
        #
        def get_neg_weights_dict_h_rann(self, node: int):
            random_neighbor = random.choice(list(self.lattice.H.neighbors(node)))
            dictH = {
                (node, random_neighbor): -1
            }
            return dictH
        #
        def get_neg_weights_dict_h_cross(self, node: int):
            dictH = {(node, nn): -1 
                     for nn in list(self.lattice.H.neighbors(node))}
            return dictH

        def get_neg_weights_dict_h_triangle(self, node: int):
            node1 = node
            node2 = list(self.lattice.H.neighbors(node1))[0]
            common_neighbors = list(nx.common_neighbors(self.lattice.H, node1, node2))
            node3 = common_neighbors[0]

            dictH = {
                (node1, node2): -1,
                (node2, node3): -1,
                (node1, node3): -1
            }
            return dictH
        def get_neg_weights_dict_h_hexagon(self, node: int):
            nodes_in_cycle = []
            nodes_in_cycle.append(node)
            node_nn = list(self.lattice.H.neighbors(node))
            samp_node_nn_1 = random.sample(node_nn, 1)[0]
            node_nn.remove(samp_node_nn_1)
            nodes_in_cycle.append(samp_node_nn_1)
            # #
            node_nn_1 = list(self.lattice.H.neighbors(samp_node_nn_1))
            node_nn_1.remove(node)
            samp_node_nn_2 = random.sample(node_nn_1, 1)[0]
            nodes_in_cycle.append(samp_node_nn_2)
            #
            flag = True
            for nn in node_nn:
                if flag:
                    common_nn = None
                    node_nn_1b = list(self.lattice.H.neighbors(nn))
                    node_nn_1b.remove(node)
                    for i in node_nn_1b:
                        common_neighs = list(nx.common_neighbors(self.lattice.H, samp_node_nn_2, i))
                        if common_neighs != []:
                            common_nn = i
                            nodes_in_cycle.append(nn)
                            flag = False
                            break
            nodes_in_cycle.append(common_nn)
            nodes_in_cycle.append(common_neighs[0])
            # # nodes_in_cycle.append()
            # listofGnodes = [self.lattice.invnode_map[k] for k in nodes_in_cycle]
            subH = self.lattice.H.subgraph(nodes_in_cycle)
            links = [tuple(sorted(edge)) for edge in subH.edges()]

            dictH = {
                link: -1 for link in links
            }
            return dictH
        #
        def get_neg_weights_dict_h_pattern(self, mode: str):
            if mode == "cross":
                merged_list = [item for i in self.node_flip_selection for item in self.get_neg_weights_dict_h_cross(i).items()]
            elif mode == "square":
                merged_list = [item for i in self.node_flip_selection for item in self.get_neg_weights_dict_h_square(i).items()]
            elif mode == "triangle":
                merged_list = [item for i in self.node_flip_selection for item in self.get_neg_weights_dict_h_triangle(i).items()]
            elif mode == "hexagon":
                merged_list = [item for i in self.node_flip_selection for item in self.get_neg_weights_dict_h_hexagon(i).items()]
            return dict(merged_list)
        #
        def get_neg_weights_dict_h_rball(self, R: int = 5):
            neighs_to_flip = get_neighbors_within_distance(self.lattice.H, self.midway_H, R)
            dictH = {(node, neighbor): -1 for node in neighs_to_flip for neighbor in self.lattice.H.neighbors(node)}
            return dictH
        #
        def get_neg_weights_dict_h_rball_EXT(self, R: int = 5):
            neighs_to_flip = get_neighbors_within_distance(self.lattice.H, self.midway_H, R)

            # Get the boundary nodes
            boundary_nodes = set()
            for node in neighs_to_flip:
                for neighbor in self.lattice.H.neighbors(node):
                    if neighbor not in neighs_to_flip:
                        boundary_nodes.add(node)
                        break

            # Create the dictionary with modified outer links
            dictH = {}
            for node in boundary_nodes:
                for neighbor in self.lattice.H.neighbors(node):
                    if neighbor not in neighs_to_flip:
                        dictH[(node, neighbor)] = -1
                    else:
                        dictH[(node, neighbor)] = 0

            return dictH
    #
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