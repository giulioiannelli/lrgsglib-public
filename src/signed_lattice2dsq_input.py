from LRGSG_package.LRGSG import *
warnings.filterwarnings("ignore")
#
description = """
    Compute entropy of the density matrix pdf and the Signed Laplacian spectrum
    associated to a Random Bond Ising Model.
    Must be executed from parent folder "/path/to/LRG-Signed", i.e. called as
    python3 src/signed_lattice2dsq_input.py ARGUMENTS
"""
hPars_L = """
    (int) system linear size (L = sqrt(N))
"""
hPars_p = """
    (float) fraction of spins to flip 
"""
default_nA = "| default=1000"
hPars_nA = f"""
    (int) number of averages {default_nA:->10}
"""
default_st = "| default=1500"
hPars_st = f"""
    (int) number of steps in the specific heat sampling {default_st:->10}
"""
default_geom = "| default=squared"
hPars_geom = f"""
    (str) geometry of the 2d lattice {default_geom:->10}
"""
default_t1 = "| default=-2"
hPars_t1 = f"""
    (float) lower exponent of the time scale: {default_t1:->10}
"""
default_t2 = "| default=5"
hPars_t2 = f"""
    (float) lower exponent of the time scale: {default_t2:->10}
"""
prsr = argparse.ArgumentParser(description=description)
prsr.add_argument('L', type=int, help=hPars_L)
prsr.add_argument('p', type=float, help=hPars_p)
prsr.add_argument('-t1', '--lowert_exponent', type=float, help=hPars_t1, default=-2)
prsr.add_argument('-t2', '--highert_exponent', type=float, help=hPars_t2, default=5)
prsr.add_argument('-nA', '--number_of_averages',type=int, help=hPars_nA, default=1000)
prsr.add_argument('-st', '--sampling_steps', type=int, help=hPars_st, default=1500)
prsr.add_argument('-geo', '--geometry', type=str, help=hPars_geom, default='squared')
prsr.set_defaults(feature=True)
args = prsr.parse_args()
#
STEPS = args.sampling_steps
GEOMETRY = args.geometry
SPIKE_THRESHOLD = 0.05
#
sideL = args.L
N = sideL * sideL
nA = args.number_of_averages
pflip = args.p
#
if not os.path.isdir((path := f"{datPath_l2d(GEOMETRY)}N={N}_navg={nA}/")):
    os.makedirs(path)
#
lattice = Lattice2D(side1 = sideL, geometry = GEOMETRY)
#
savename = lambda idstr: f"{path}p={pflip:{pflip_fmt}}_{idstr}{eBIN}"
if os.path.exists(savename('Sm1')) and os.path.exists(savename('slspec')):
    exit()
Sm1File = open(savename('Sm1'), "wb")
slspecFile = open(savename('slspec'), "wb")
#
for nr in range(nA):
    random.seed(nr)
    while True:
        SLRG_obj = SignedLaplacianAnalysis(#
            system = lattice,
            pflip = pflip,
            t1 = args.lowert_exponent,
            t2 = args.highert_exponent,
            steps = STEPS
        )
        SLRG_obj.flip_random_fract_edges()
        SLRG_obj.compute_entropy()
        index = np.where(np.diff(SLRG_obj.Cspe) > SPIKE_THRESHOLD)
        if not index[0].size:
            break
    #
    Sm1File.write(bytes(SLRG_obj.Sm1))
    slspecFile.write(bytes(SLRG_obj.slspectrum))
Sm1File.close()
slspecFile.close()
