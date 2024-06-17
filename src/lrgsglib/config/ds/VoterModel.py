from ..dyns import *

class VoterModel(BinDynSys):
    dyn_UVclass = "voter_model"
    eqSTEP_def = 10
    id_string_voterdyn = ""

    def __init__(self, sg: SignedGraph = Lattice2D, **kwargs) -> None:
        self.sg = sg
        self.dynpath = f"{self.sg.DEFAULT_VOTERDIR}{self.sg.syshapePth}"
        super(VoterModel, self).__init__(self.sg, **kwargs)

    def ds1step(self, nd: int):
        nodedict = dict(self.sg.H[nd])
        neighs = list(nodedict.keys())
        nn = np.random.choice(neighs)
        self.s[nd] = nodedict[nn]["weight"] * self.s[nn]

    def run_py(self):
        dsNstep = self.dsNstep()
        nodes = list(self.sg.H.nodes())
        for _ in range(self.simtime):
            smp = sample(nodes, self.sg.N)
            self.s_t.append(self.s.copy())
            # self.ds1step(np.random.randint(self.sg.N))
            dsNstep(smp)
    
    def run_c(self,
              adjfname: str = "",
              out_suffix: str = "",
              eqSTEP: int = 0):
        if eqSTEP:
            self.eqSTEP_def = eqSTEP
        if adjfname == "":
            adjfname = self.sg.stdFname
        out_suffix = out_suffix + self.id_string_voterdyn
        if out_suffix == "":
            out_suffix = '\'\''
        self.cprogram = [
            # pth_join(DIR_SRC, DIR_PCK, self.dyn_UVclass),
            # f"{DIR_SRC_DEFAULT}{DIR_PCK_DEFAULT}{self.dyn_UVclass}",
            f"{self.sg.N}",
            f"{self.sg.pflip}",
            f"{self.eqSTEP_def}",
            self.sg.datPath,
            adjfname,
            out_suffix
        ]
        subprocess.call(self.cprogram)

    def run(self, **kwargs):
        if self.runlang.startswith("py"):
            self.run_py()
        elif self.runlang.startswith("C"):
            self.run_c(**kwargs)



