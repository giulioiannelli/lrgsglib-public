from LRGSG_package.LRGSG import *
import warnings
warnings.filterwarnings("ignore")
#
STEPS = 1500
SPIKE_THRESHOLD = 0.1
#
pow2_m, pow2_M, pow2_A = 6, 14, 15
ssize_range = range(pow2_m, pow2_M, 2)
asize_range = range(pow2_m, pow2_M, 1)[:len(ssize_range)]
#
lsN = np.array([2**i for i in ssize_range])
lsL = np.sqrt(lsN).astype('int')
lsA = np.array([2**(pow2_A-i) for i in asize_range])
lsp = np.round([0.001, 0.01, 0.025, 0.05, 0.075, 0.085, 0.095, 0.098, 0.100, 
                0.102, 0.104, 0.106, 0.11, 0.2, 0.3, 0.4, 0.5], 3)
#
for iL, iN, iA in zip(lsL, lsN, lsA):
    #
    path = f"{datPath_l2d_sq}N={iN}_navg={iA}/"
    if not os.path.isdir(path):
        os.makedirs(path)
    #
    sqLattice = Lattice2D(#
        side1 = iL,
        geometry = 'squared'
    )
    for pflip in lsp:
        savename = lambda idstr : f"{path}p={pflip:{pflip_fmt}}_{idstr}{eBIN}"
        if os.path.exists(savename('Sm1')):
            continue
        Sm1File = open(savename('Sm1'), "wb")
        CspeFile = open(savename('Cspe'), "wb")
        slspecFile = open(savename('slspec'), "wb")
        #
        for nr in tqdm(range(iA), desc=f"replicas for L={iL}, p={pflip}"):
            while True:
                SLRG_obj = SignedLaplacianAnalysis(#
                    system = sqLattice,
                    pflip = pflip,
                    t2 = 4,
                    steps = STEPS
                )
                SLRG_obj.flip_random_fract_edges()
                SLRG_obj.compute_entropy()
                index = np.where(np.diff(SLRG_obj.Cspe) > SPIKE_THRESHOLD)
                if not index[0].size:
                    break
            #
            Sm1BytesArray = bytes(SLRG_obj.Sm1)
            CspeBytesArray = bytes(SLRG_obj.Cspe)
            slspecBytesArray = bytes(SLRG_obj.slspectrum)
            #
            Sm1File.write(Sm1BytesArray)
            CspeFile.write(CspeBytesArray)
            slspecFile.write(slspecBytesArray)
        Sm1File.close()
        CspeFile.close()
        slspecFile.close()
