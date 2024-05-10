from ..shared import *
from ..config.const import *
from ..config.utils import *
from ..nx_patches.funcs import *
from ..nx_patches.objects import *


class BinStocProc:
    def __init__(
        self,
        sg: SignedGraph,
        simTime: int = 1,
        initMode: str = "",
        runMode: str = BSP_RUN_MODE,
        storeMode: str = BSP_STORE_MODE,
        storeFreq: int = BSP_STORE_FREQ,
    ):
        self.s_t = []
        self.s_tList = []
        self.s_tCntr = Counter()
        self.sg = sg
        self.initMode = initMode
        if runMode in BSP_RUN_MODES:
            self.runMode = BSP_RUN_MODES_DICT[runMode]
        else:
            raise ValueError(
                f'runMode "{runMode}" is not in the allowed list of strings'
            )
        if storeMode in BSP_STORE_MODES:
            self.storeMode = storeMode
            if storeFreq > 0:
                self.storeFreq = storeFreq
            else:
                raise ValueError("storeFreq must be a positive integer")
        else:
            raise ValueError(
                f'storeMode "{storeMode}" is not in the allowed list of strings'
            )
        self.simTime = simTime

    def __run_py_1__(self):
        raise NotImplementedError("Subclasses must implement this method")

    def __run_py_N__(self):
        raise NotImplementedError("Subclasses must implement this method")

    def __run_py__(self):
        raise NotImplementedError("Subclasses must implement this method")

    def run(self):
        raise NotImplementedError("Subclasses must implement this method")


class SignedRW(BinStocProc):
    def __init__(self, sg: SignedGraph, **kwargs):
        super().__init__(sg, **kwargs)
        if self.initMode:
            self.__init_proc__()
        if self.storeMode:
            if self.storeMode == "persistent":
                self.store_state = self.store_state_per
            elif self.storeMode == "sequential":
                self.store_state = self.store_state_seq
        else:
            self.store_state = do_nothing
        self.node_index = {node: i for i,node in enumerate(self.sg.nodeList)}
        self.index_node = {i: node for i,node in enumerate(self.sg.nodeList)}


    def __init_random__(self):
        pos = random.choice(self.sg.nodeList)
        val = random.choice(BSP_VAL)
        self.s_0 = [pos, val]

    def __init_proc__(self):
        if self.initMode == "random":
            self.__init_random__()
        else:
            raise ValueError("Invalid initMode")
        self.s_t = self.s_0
        self.s_tList.append(self.s_t)

    def __run_py_1__(self):
        adj_matrix = self.sg.Adj
        ndidx = self.node_index[self.s_t[0]]
        start, end = adj_matrix.indptr[ndidx], adj_matrix.indptr[ndidx+1]
        neighbors_indices = adj_matrix.indices[start:end]
        cnn = np.random.choice(neighbors_indices)
        self.s_t = [self.index_node[cnn], self.s_t[1] * adj_matrix[ndidx, cnn]]

    def __run_py_N__(self):
        for t in range(self.sg.N):
            self.__run_py_1__()
            self.store_state(t)

    def __run_py__(self):
        for _ in range(self.simTime):
            self.__run_py_N__()

    def run(self):
        if self.runMode in BSP_RUN_MODES_PY_LIST:
            self.__run_py__()
        else:
            raise ValueError("Invalid runMode")

    def store_state_seq(self, t: int):
        if not (t % self.storeFreq):
            self.s_tList.append(self.s_t)
    def store_state_per(self, t: int):
        if not (t % self.storeFreq):
            self.s_tCntr[self.s_t[0]] += self.s_t[1]
