from .nx_objects import *

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
        dim: tuple = DEFLattice3D_dim,
        pbc: bool = DEFLattice3D_pbc,
        geo: str = 'cubic',
        fbc_val: float = DEFLattice3D_fbcv,
        stdFnameSFFX: str = DEFLattice3D_stdFnsfx,
        sgpath: str = DEFLattice3D_sgpath,
        with_positions: bool = False,
        theta: float = np.pi/2,
        phi: float = np.pi/2,
        **kwargs,
    ) -> None:
        self.dim = sorted(dim, reverse=True)
        self.diml = list(dim)
        self.pbc = pbc
        self.geo = geo
        self.theta = theta
        self.phi = phi
        self.fbc_val = fbc_val
        self.__init_stdFname__(stdFnameSFFX)
        _ = L3D_PATH_DICT[self.geo]
        self.sgpath = sgpath + _ if sgpath else _
        self.with_positions = with_positions
        self.__init_lattice__()
        super(Lattice3D, self).__init__(self.G, **kwargs)

    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = DEFLattice3D_stdFn + SFFX

    def __init_lattice__(self) -> None:
        if self.geo == 'cubic':
            nxfunc = self._generate_cubic_lattice
        elif self.geo == 'bcc':
            nxfunc = self._generate_bcc_lattice
        elif self.geo == 'fcc':
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

    def _project_3d_to_2d(self, x, y, z):
        # Rotation matrix around the y-axis (theta)
        R_theta = np.array([
            [np.cos(self.theta), 0, np.sin(self.theta)],
            [0, 1, 0],
            [-np.sin(self.theta), 0, np.cos(self.theta)]
        ])
        
        # Rotation matrix around the x-axis (phi)
        R_phi = np.array([
            [1, 0, 0],
            [0, np.cos(self.phi), -np.sin(self.phi)],
            [0, np.sin(self.phi), np.cos(self.phi)]
        ])
        
        # Initial position vector
        position = np.array([x, y, z])
        
        # Apply rotations
        position_rotated = R_phi @ R_theta @ position  # Order matters
        
        # Project onto 2D plane (ignore z after rotation)
        x2, y2 = position_rotated[0], position_rotated[1]
        
        return x2, y2

    def _set_positions(self):
        pos = {node: self._project_3d_to_2d(*node) for node in self.G.nodes()}
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
        range_adjust = 0 if self.pbc else -1
        nodes = list(cProd_Iter(dim)) + \
                [(x + 0.5, y + 0.5, z + 0.5) 
                        for x, y, z in cProd_Iter_adj(dim, range_adjust)]
        G.add_nodes_from(nodes)

        edges = [((x + dx, y + dy, z + dz), (x + ddx + dx, y + ddy + dy, z + ddz + dz))
                for x, y, z in cProd_Iter(dim)  # Assuming self.cprod includes all relevant nodes
                for dx, dy, dz in offsets
                for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
                if (x + ddx + dx, y + ddy + dy, z + ddz + dz) in G.nodes()]

        G.add_edges_from(edges)

        # for x in range(self.dim[0]):
        #     for y in range(self.dim[1]):
        #         for z in range(self.dim[2]):
        #             for dx, dy, dz in offsets:
        #                 nx, ny, nz = x + dx, y + dy, z + dz
        #                 for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        #                     neighbor = (nx + ddx, ny + ddy, nz + ddz)
        #                     if neighbor in G.nodes():
        #                         G.add_edge((nx, ny, nz), neighbor)
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
            # for x in range(self.dim[0]-1):
            #     for y in range(self.dim[1]-1):
            #         G.add_edge((x+0.5, y+0.5, 0.5), (x+0.5, y+0.5, self.dim[2]-1+0.5))
            #         G.add_edge((x, y, 0), (x, y, self.dim[2]-1))
            # for x in range(self.dim[0]):
            #     for z in range(self.dim[2]):
            #         G.add_edge((x+0.5, 0.5, z+0.5), (x+0.5, self.dim[1]-1+0.5, z+0.5))
            #         G.add_edge((x, 0, z), (x, self.dim[1]-1, z))
            # for y in range(self.dim[1]):
            #     for z in range(self.dim[2]):
            #         G.add_edge((0.5, y+0.5, z+0.5), (self.dim[0]-1+0.5, y+0.5, z+0.5))
            #         G.add_edge((0, y, z), (self.dim[0]-1, y, z))
        return G

    def _generate_fcc_lattice(self):
        G = Graph()
        offsets = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        for x in range(self.dim[0]):
            for y in range(self.dim[1]):
                for z in range(self.dim[2]):
                    for dx, dy, dz in offsets:
                        nx, ny, nz = x + dx, y + dy, z + dz
                        G.add_node((nx, ny, nz))
                        for ddx, ddy, ddz in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]:
                            neighbor = (nx + ddx, ny + ddy, nz + ddz)
                            if neighbor in G.nodes():
                                G.add_edge((nx, ny, nz), neighbor)
        if self.pbc:
            # For FCC lattice PBC
            for x in range(-1, self.dim[0], 1):
                for y in range(self.dim[1]):
                    G.add_edge((x+0.5, y+0.5, 0), (x+0.5, y+0.5, self.dim[2]-0.5))
                    G.add_edge((x, y, 0), (x, y, self.dim[2]-0.5))
            for x in range(self.dim[0]):
                for z in range(-1, self.dim[2], 1):
                    G.add_edge((x+0.5, self.dim[1]-0.5, z+0.5), (x+0.5, -0.5, z+0.5))
                    G.add_edge((x, self.dim[1]-1, z), (x, 0, z))
            for y in range(-1, self.dim[1], 1):
                for z in range(self.dim[2]):
                    G.add_edge((self.dim[0]-0.5, y+0.5, z+0.5), (-0.5, y+0.5, z+0.5))
                    G.add_edge((self.dim[0]-1, y, z), (0, y, z))
        return G
    class nwContainer(dict):
            def __init__(self, l: SignedGraph, iterable=[], constant=None, 
                        **kwargs):
                super().__init__(**kwargs)
                self.update((key, constant) for key in iterable)
                self.l = l
                self.reprDict = list(self.l.GraphReprDict.keys())