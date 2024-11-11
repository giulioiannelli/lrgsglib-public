from .common import *
from ..config.funcs import project_3d_to_2d
from .funcs import *
from .SignedGraph import SignedGraph

class Lattice3D(SignedGraph):
    def __init__(
        self,
        dim: Union[int, Tuple[int, int, int]] = L3D_DIM,
        geo: str = L3D_GEO,
        pbc: bool = L3D_PBC,
        fbc_val: float = L3D_FBCV,
        stdFnameSFFX: str = L3D_STDFN,
        sgpathn: str = L3D_SGPATH,
        with_positions: bool = L3D_WITH_POS,
        theta: float = L3D_THETA,
        phi: float = L3D_PHI,
        pdil: float = L3D_PDIL,
        **kwargs,
    ) -> None:
        if isinstance(dim, (int, np.integer)):
            self.dim = (dim, dim, dim)
        elif isinstance(dim, tuple) and len(dim) == 3 and all(isinstance(x, int) for x in dim):
            self.dim = sorted(dim, reverse=True)
        else:
            raise ValueError("dim must be an integer or a tuple of three integers")
        self.dimL = list(self.dim)
        self.pbc = pbc
        self.pdil = pdil
        self.fbc_val = fbc_val
        self.geo = L3D_GEO_DICT[geo]
        self.theta = theta
        self.phi = phi
        self.__init_stdFname__(stdFnameSFFX)
        _ = L3D_PATH_DICT[self.geo]
        self.sgpathn = sgpathn + _ if sgpathn else _
        self.with_positions = with_positions
        self.__init_lattice__()
        super(Lattice3D, self).__init__(self.G, **kwargs)

    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = L3D_STDFN + SFFX

    def __init_lattice__(self) -> None:
        if self.geo == L3D_GEO_SC:
            if self.pdil == 0.:
                nxfunc = LatticeND_graph_FastPatch
            else:
                nxfunc = compose(LatticeND_graph_FastPatch, 
                                 remove_edges, 
                                 g_kwargs={'pdil': self.pdil})
        elif self.geo == L3D_GEO_BCC:
            nxfunc = self._generate_bcc_lattice
        elif self.geo == L3D_GEO_FCC:
            nxfunc = self._generate_fcc_lattice
        else:
            raise ValueError(f"Unsupported geometry '{self.geo}'.")
        self.syshape = self.dim
        if all(x == self.dim[0] for x in self.dim):
            self.syshapePth = f"N={self.dim[0]**len(self.dim)}"
        else:
            self.syshapePth = '_'.join([f"L{i}={side}" 
                                        for i, side in enumerate(self.dim)])        
        
        self.G = nxfunc(self.dim, periodic=self.pbc)
        if self.with_positions:
            self._set_positions()


    def _set_positions(self):
        pos = {node: project_3d_to_2d(*node, self.theta, self.phi) 
               for node in self.G.nodes()}
        nx.set_node_attributes(self.G, pos, 'pos')


    def _generate_bcc_lattice(self, dim):
        G = Graph()
        offsets = [(0, 0, 0), (0.5, 0.5, 0.5)]
        offsets2 = [(0.5, 0.5, 0.5), (0.5, 0.5, -0.5),
                    (0.5, -0.5, 0.5), (-0.5, 0.5, 0.5),
                    (0.5, -0.5, -0.5), (-0.5, -0.5, -0.5),
                    (-0.5, 0.5, -0.5), (-0.5, -0.5, 0.5)]
        range_adjust = 0 if self.pbc else -1
        nodes = [(x + ox, y + oy, z+ oz) 
                            for ox, oy, oz in offsets
                            for x, y, z in cProd_Iter_adj(dim, range_adjust)]
        G.add_nodes_from(nodes)

        edges = [((x + dx, y + dy, z + dz), (x + ddx + dx, y + ddy + dy, z + ddz + dz), {'type': 'link'})
                for x, y, z in cProd_Iter(dim) 
                for dx, dy, dz in offsets
                for ddx, ddy, ddz in offsets2
                if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

        G.add_edges_from(edges)
        edges = [((x + dx, y + dy, z + dz), (x + ddx + dx, y + ddy + dy, z + ddz + dz), {'type': 'box'})
                for x, y, z in cProd_Iter(dim) 
                for dx, dy, dz in offsets
                for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
                if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

        G.add_edges_from(edges)
        if self.pbc:
            # For BCC lattice PBC
            G.add_edges_from([((x, y, 0), (x+0.5 , y +0.5, dim[2] - 0.5), {'type': 'link'}) 
                            for x, y in cProdSel_Iter(dim, (0, 1))])
            G.add_edges_from([((x, 0, 0), (x+0.5 , dim[1]-0.5, dim[2] - 0.5), {'type': 'link'}) 
                            for x in range(dim[0])])
            G.add_edges_from([((0, y, 0), (dim[0]-0.5 , y+0.5, dim[2] - 0.5), {'type': 'link'}) 
                            for y in range(dim[1])])
            G.add_edges_from([((0, 0, z), (dim[0]-0.5, dim[1]-0.5, z+0.5), {'type': 'link'}) 
                            for z in range(dim[2])])
            G.add_edges_from([((x, 0, z), (x+0.5 , dim[1]-0.5, z+0.5), {'type': 'link'}) 
                            for x, z in cProdSel_Iter(dim, (0, 2))])
            G.add_edges_from([((0, y, z), (dim[0]-0.5 , y+0.5, z+0.5), {'type': 'link'}) 
                            for y, z in cProdSel_Iter(dim, (1, 2))])
            # G.add_edges_from([((x + ox, oy, z + oz), (x + ox, self.dim[1] - 1 + oy, z + oz))
            #                   for ox, oy, oz in offsets
            #                   for x, z in cProdSel_Iter(dim, (0, 2))])
            # G.add_edges_from([((ox, y + oy, z + oz), (self.dim[0] - 1 + ox, y + oy, z + oz))
            #                 for ox, oy, oz in offsets
            #                 for y, z in cProdSel_Iter(dim, (1, 2))])
        return G
    
    def _generate_fcc_lattice(self, dim):
        G = Graph()
        offsets = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        range_adjust = 0 if self.pbc else -1
        nodes = [(x + ox, y + oy, z+ oz) 
                            for ox, oy, oz in offsets
                            for x, y, z in cProd_Iter_adj(dim, range_adjust)]
        G.add_nodes_from(nodes)

        edges = [((x + dx, y + dy, z + dz), (x + ddx + dx, y + ddy + dy, z + ddz + dz))
                for x, y, z in cProd_Iter(dim) 
                for dx, dy, dz in offsets
                for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
                if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

        G.add_edges_from(edges)
        if self.pbc:
            # For BCC lattice PBC
            G.add_edges_from([((x + ox, y + oy, oz), (x + ox, y + oy, dim[2] - 1 + oz)) 
                            for ox, oy, oz in offsets
                            for x, y in cProdSel_Iter(dim, (0, 1))])
            G.add_edges_from([((x + ox, oy, z + oz), (x + ox, self.dim[1] - 1 + oy, z + oz))
                              for ox, oy, oz in offsets
                              for x, z in cProdSel_Iter(dim, (0, 2))])
            G.add_edges_from([((ox, y + oy, z + oz), (self.dim[0] - 1 + ox, y + oy, z + oz))
                            for ox, oy, oz in offsets
                            for y, z in cProdSel_Iter(dim, (1, 2))])
        return G
    # def _generate_fcc_lattice(self, dim):
    #     G = Graph()
        
    #     # Offsets for FCC lattice to include face-centered points
    #     fcc_offsets = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        
    #     # Adjusted range to ensure we don't go beyond the lattice dimensions
    #     range_adjusted = [range(d) for d in self.dim]

    #     # Generate all nodes, including shifted nodes for FCC structure
    #     nodes = [(x + ox, y + oy, z + oz) 
    #              for x, y, z in product(*range_adjusted)
    #              for ox, oy, oz in fcc_offsets]
        
    #     # Normalize nodes to ensure they are within bounds (for non-PBC case)
    #     # This step also removes duplicates from adding offsets that exceed dimensions
    #     normalized_nodes = set()
    #     for node in nodes:
    #         if not self.pbc:  # Without PBC, ensure nodes are within the lattice bounds
    #             nx, ny, nz = node
    #             if 0 <= nx < self.dim[0] and 0 <= ny < self.dim[1] and 0 <= nz < self.dim[2]:
    #                 normalized_nodes.add(node)
    #         else:  # With PBC, we might normalize or wrap the nodes differently
    #             normalized_nodes.add(node)

    #     G.add_nodes_from(normalized_nodes)

    #     node_list = list(G.nodes())
    #     # Convert node coordinates to a numpy array for vectorized operations
    #     coordinates = np.array(node_list)

    #     # Calculate Euclidean distances using broadcasting
    #     # The result is a matrix of shape (len(coordinates), len(coordinates))
    #     distances = np.sqrt(((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]) ** 2).sum(axis=2))
    #     for i,nd in enumerate(node_list):
    #         idx = np.where(distances[i] < 1.)[0]
    #         neighs = np.array(node_list, dtype=object)[idx]
    #         neighs = [tuple(ne) for ne in neighs]
    #         edges = ((nd, ne) for ne in neighs if not G.has_edge(ne, nd) and ne != nd)
    #         G.add_edges_from(edges)
    #     # Now, G contains all FCC lattice nodes, and you can proceed to add edges based on your criteria.

    #     return G
    # def _generate_fcc_lattice(self, dim):
    #     G = Graph()
        
    #     # FCC offsets for directly adjacent neighbors
    #     fcc_offsets = [(0.5, 0.5, 0), (0.5, -0.5, 0), (0.5, 0, 0.5), 
    #                    (0.5, 0, -0.5), (-0.5, 0.5, 0), (-0.5, -0.5, 0), 
    #                    (-0.5, 0, 0.5), (-0.5, 0, -0.5), (0, 0.5, 0.5), 
    #                    (0, 0.5, -0.5), (0, -0.5, 0.5), (0, -0.5, -0.5)]

    #     # Generate nodes
    #     nodes = list(product(range(self.dim[0]), range(self.dim[1]), range(self.dim[2])))
    #     G.add_nodes_from(nodes)

    #     node_list = list(G.nodes())
    #     # Convert node coordinates to a numpy array for vectorized operations
    #     coordinates = np.array(node_list)

    #     # Calculate Euclidean distances using broadcasting
    #     # The result is a matrix of shape (len(coordinates), len(coordinates))
    #     distances = np.sqrt(((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]) ** 2).sum(axis=2))
    #     for i,nd in enumerate(node_list):
    #         idx = np.where(distances[i] < 1.25)[0]
    #         neighs = np.array(node_list, dtype=object)[idx]
    #         neighs = [tuple(ne) for ne in neighs]
    #         edges = ((nd, ne) for ne in neighs if not G.has_edge(ne, nd))
    #         G.add_edges_from(edges)

    #     # # Add edges based on FCC connectivity
    #     # for x, y, z in nodes:
    #     #     for dx, dy, dz in fcc_offsets:
    #     #         nx, ny, nz = x + dx, y + dy, z + dz
                
    #     #         # Handle wrapping for PBC; otherwise, ensure the neighbor is within bounds
    #     #         if self.pbc:
    #     #             nx %= self.dim[0]
    #     #             ny %= self.dim[1]
    #     #             nz %= self.dim[2]
    #     #         elif not (0 <= nx < self.dim[0] and 0 <= ny < self.dim[1] and 0 <= nz < self.dim[2]):
    #     #             continue  # Skip adding this edge if it's out of bounds without PBC
                
    #     #         G.add_edge((x, y, z), (nx, ny, nz))

    #     return G

    # def _generate_fcc_lattice(self, dim):
    #     G = Graph()
    #     offsets = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    #     range_adjust = 0 if self.pbc else -1
    #     nodes = [(x + ox, y + oy, z+ oz) 
    #              for ox, oy, oz in offsets
    #              for x, y, z in cProd_Iter_adj(dim, range_adjust)]
    #     G.add_nodes_from(nodes)

    #     for x, y, z in cProd_Iter(dim):
    #         for dx, dy, dz in offsets:
    #             nx, ny, nz = x + dx, y + dy, z + dz
    #             for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5), (0.5, -0.5, 0), (0.5, 0, -0.5)]:
    #                 neighbor = (nx + ddx, ny + ddy, nz + ddz)
    #                 if neighbor in G.nodes() and (nx, ny, nz) in G.nodes() and (nx, ny, nz) != neighbor:
    #                     G.add_edge((nx, ny, nz), neighbor)
    #     # if self.pbc:
    #         # For FCC lattice PBC
    #         # for x in range(-1, self.dim[0], 1):
    #         #     for y in range(self.dim[1]):
    #         #         G.add_edge((x+0.5, y+0.5, 0), (x+0.5, y+0.5, self.dim[2]-0.5))
    #         #         G.add_edge((x, y, 0), (x, y, self.dim[2]-0.5))
    #         # for x in range(self.dim[0]):
    #         #     for z in range(-1, self.dim[2], 1):
    #         #         G.add_edge((x+0.5, self.dim[1]-0.5, z+0.5), (x+0.5, -0.5, z+0.5))
    #         #         G.add_edge((x, self.dim[1]-1, z), (x, 0, z))
    #         # for y in range(-1, self.dim[1], 1):
    #         #     for z in range(self.dim[2]):
    #         #         G.add_edge((self.dim[0]-0.5, y+0.5, z+0.5), (-0.5, y+0.5, z+0.5))
    #         #         G.add_edge((self.dim[0]-1, y, z), (0, y, z))
    #     return G
    def get_central_edge(self, on_g: str = L3D_ONREP):
        cnode = (self.dimL[0]//2-1, self.dimL[1]//2, self.dimL[2]//2)
        cnode_t = (self.dimL[0]//2, self.dimL[1]//2, self.dimL[2]//2)
        edge_t = (cnode, cnode_t)
        if on_g == 'G':
            return edge_t
        elif on_g == 'H':
            return self.edgeMap['H']['G'][edge_t]
    class nwContainer(dict):
            def __init__(self, l: SignedGraph, iterable=[], constant=None, 
                        **kwargs):
                super().__init__(**kwargs)
                self.update((key, constant) for key in iterable)
                self.l = l
                self.rd = self.l.GraphReprs
                #
                self.centedge = {g: self.l.get_central_edge(g) 
                             for g in self.rd}
                self.reprDict = list(self.l.GraphReprDict.keys())
                self['single'] = {g: [self.centedge[g]] for g in self.rd}
                self['rand'] = {g: [e for e in self.l.fleset[g]] 
                            for g in self.rd}