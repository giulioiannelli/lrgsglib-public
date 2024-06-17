from ..dyns import *
class BinDynSys:
    s_t = []
    dynpath = ""

    def __init__(
        self,
        sg: SignedGraph,
        ic: str = "uniform",
        runlang: str = "py",
        simpref: int = 1,
        savedyn: bool = False,
        seed: int = None,
    ):
        if seed:
            random.seed(seed)
            np.random.seed(seed)
        self.sg = sg
        self.ic = ic
        self.runlang = runlang
        self.simtime = simpref * sg.N
        self.savedyn = savedyn
        self.init_s()

    #
    def export_sinit(self):
        outf = open(f"{self.dynpath}s_{self.sg.stdFname}.bin", "wb")
        self.s.astype("int8").tofile(outf)

    def init_s(self):
        if self.ic == "uniform":
            self.s = np.random.choice([-1, 1], size=self.sg.N)
        elif self.ic.startswith("ground_state"):
            neig = int(self.ic.split("_")[-1])
            self.s = self.sg.bin_eigV(which=neig)
        if self.runlang != "py" and not self.sg.import_on:
            self.export_sinit()

    def ds1step(self, nd: int):
        raise NotImplementedError("Subclasses must implement this method")

    def dsNstep(self):
        return np.vectorize(self.ds1step, excluded="self")

    def run(self):
        raise NotImplementedError("Subclasses must implement this method")

    def run_py(self):
        raise NotImplementedError("Subclasses must implement this method")

    def run_c(self):
        raise NotImplementedError("Subclasses must implement this method")
