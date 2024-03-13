from .nx_objects import *

class Lattice3D(SignedGraph):
    #
    def __init__(
        self,
        dim: tuple = DEFLattice3D_dim,
        pbc: bool = DEFLattice3D_pbc,
        fbc_val: float = DEFLattice3D_fbcv,
        stdFnameSFFX: str = DEFLattice3D_stdFnsfx,
        sgpath: str = DEFLattice3D_sgpath,
        with_positions: bool = False,
        **kwargs,
    ) -> None:
        self.dim = sorted(dim, reverse=True)
        self.diml = list(dim)
        self.pbc = pbc
        self.fbc_val = fbc_val
        self.__init_stdFname__(stdFnameSFFX)
        #
        _ = '['+('{},'*(len(dim)-1)+'{}').format(*dim)+']'
        self.sgpath = sgpath + _ if sgpath else _
        #
        self.__init_lattice__()
        super(Lattice3D, self).__init__(self.G, **kwargs)
    #
    def __init_stdFname__(self, SFFX: str = "") -> None:
        self.stdFname = DEFLattice3D_stdFn + SFFX
    #
    def __init_lattice__(self) -> None:
        #
        nxfunc = nx.grid_graph
        #
        if all(dim[i] == dim[i+1] for i,dim in enumerate(self.diml[:-1])):
            self.syshapeStr = f"N={np.prod(self.dim)}"
        else:
            self.syshapeStr = '_'.join('L{}={}'.format(i+1, self.dim[i]) 
                                       for i in range(len(self.dim)))
        self.syshapePth = f"{self.syshapeStr}/"
        #
        self.p_c = 0.222
        self.r_c = np.sqrt(0.3418/(np.pi*self.p_c))
        #
        self.G = nxfunc(dim=self.dim, periodic=self.pbc)
    #