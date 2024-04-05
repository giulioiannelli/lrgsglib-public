from .objects import *

# class Lattice3D(SignedGraph):
#     #
#     def __init__(
#         self,
#         dim: tuple = DEFLattice3D_dim,
#         pbc: bool = DEFLattice3D_pbc,
#         fbc_val: float = DEFLattice3D_fbcv,
#         stdFnameSFFX: str = DEFLattice3D_stdFnsfx,
#         sgpath: str = DEFLattice3D_sgpath,
#         with_positions: bool = False,
#         **kwargs,
#     ) -> None:
#         self.dim = sorted(dim, reverse=True)
#         self.diml = list(dim)
#         self.pbc = pbc
#         self.fbc_val = fbc_val
#         self.__init_stdFname__(stdFnameSFFX)
#         #
#         _ = '['+('{},'*(len(dim)-1)+'{}').format(*dim)+']'
#         self.sgpath = sgpath + _ if sgpath else _
#         #
#         self.__init_lattice__()
#         super(Lattice3D, self).__init__(self.G, **kwargs)
#     #
#     def __init_stdFname__(self, SFFX: str = "") -> None:
#         self.stdFname = DEFLattice3D_stdFn + SFFX
#     #
#     def __init_lattice__(self) -> None:
#         #
#         nxfunc = nx.grid_graph
#         #
#         if all(dim[i] == dim[i+1] for i,dim in enumerate(self.diml[:-1])):
#             self.syshapeStr = f"N={np.prod(self.dim)}"
#         else:
#             self.syshapeStr = '_'.join('L{}={}'.format(i+1, self.dim[i]) 
#                                        for i in range(len(self.dim)))
#         self.syshapePth = f"{self.syshapeStr}/"
#         #
#         self.p_c = 0.222
#         self.r_c = np.sqrt(0.3418/(np.pi*self.p_c))
#         #
#         self.G = nxfunc(dim=self.dim, periodic=self.pbc)
    #

from networkx import Graph, set_node_attributes
import numpy as np

class Lattice3D(SignedGraph):
    def __init__(
        self,
        dim: tuple = L3D_DIM,
        pbc: bool = L3D_PBC,
        geo: str = L3D_GEO,
        fbc_val: float = L3D_FBCV,
        stdFnameSFFX: str = L3D_STDFN,
        sgpath: str = L3D_SGPATH,
        with_positions: bool = L3D_WITH_POS,
        theta: float = L3D_THETA,
        phi: float = L3D_PHI,
        **kwargs,
    ) -> None:
        self.dim = sorted(dim, reverse=True)
        self.dimL = list(self.dim)
        self.pbc = pbc
        self.fbc_val = fbc_val
        self.geo = L3D_GEO_DICT[geo]
        self.theta = theta
        self.phi = phi
        self.__init_stdFname__(stdFnameSFFX)
        _ = L3D_PATH_DICT[self.geo]
        self.sgpath = sgpath + _ if sgpath else _
        self.with_positions = with_positions
        self.__init_lattice__()
        super(Lattice3D, self).__init__(self.G, **kwargs)

    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = L3D_STDFN + SFFX

    def __init_lattice__(self) -> None:
        if self.geo == L3D_GEO_SC:
            nxfunc = self._generate_cubic_lattice
        elif self.geo == L3D_GEO_BCC:
            nxfunc = self._generate_bcc_lattice
        elif self.geo == L3D_GEO_FCC:
            nxfunc = self._generate_fcc_lattice
        else:
            raise ValueError(f"Unsupported geometry '{self.geo}'.")
        if all(x == self.dim[0] for x in self.dim):
            self.syshapeStr = f"N={self.dim[0]**len(self.dim)}"
        else:
            self.syshapeStr = '_'.join(["L{i}=side" for i, side in enumerate(self.dim)])        
        self.syshapePth = f"{self.syshapeStr}/"
        
        self.G = nxfunc(self.dim)
        if self.with_positions:
            self._set_positions()

    def _project_3d_to_2d(self, x, y, z, theta = None, phi = None):
        if theta == None:
            theta = self.theta
        if phi == None:
            phi = self.phi
        
        # Rotation matrix around the y-axis (theta)
        R_theta = np.array([
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)]
        ])
        
        # Rotation matrix around the x-axis (phi)
        R_phi = np.array([
            [1, 0, 0],
            [0, np.cos(phi), -np.sin(phi)],
            [0, np.sin(phi), np.cos(phi)]
        ])
        
        # Initial position vector
        position = np.array([x, y, z])
        
        # Apply rotations
        position_rotated = R_phi @ R_theta @ position  # Order matters
        
        # Project onto 2D plane (ignore z after rotation)
        x2, y2 = position_rotated[0], position_rotated[1]
        
        return x2, y2

    def _set_positions(self, theta = None, phi = None):
        pos = {node: self._project_3d_to_2d(*node, theta, phi) for node in self.G.nodes()}
        set_node_attributes(self.G, pos, 'pos')

    def _generate_cubic_lattice(self, dim):
        G = Graph()
        e_i = [tuple(1 if i == j else 0 for j in range(len(dim)))
               for i in range(len(dim))]
        nodes = list(cProd_Iter(dim))
        G.add_nodes_from(nodes)

        edges = [(pt, pte) for pt in nodes for drt in e_i
            if (pte := tuple(d + p for d, p in zip(drt, pt))) in G.nodes()]
        G.add_edges_from(edges)
        
        if self.pbc:
            G.add_edges_from([((x, y, 0), (x, y, dim[2]-1)) 
                              for x, y in cProdSel_Iter(dim, (0, 1))])
            G.add_edges_from([((x, 0, z), (x, dim[1]-1, z)) 
                              for x, z in cProdSel_Iter(dim, (0, 2))])
            G.add_edges_from([((0, y, z), (dim[0]-1, y, z))
                              for y, z in cProdSel_Iter(dim, (1, 2))])
        return G

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
    class nwContainer(dict):
            def __init__(self, l: SignedGraph, iterable=[], constant=None, 
                        **kwargs):
                super().__init__(**kwargs)
                self.update((key, constant) for key in iterable)
                self.l = l
                self.reprDict = list(self.l.GraphReprDict.keys())
                self['rand'] = {g: [e for e in self.l.rEdgeFlip[g]] 
                            for g in self.reprDict}