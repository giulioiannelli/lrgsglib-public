import numpy as np
import matplotlib.cm as cm
from matplotlib.cm import twilight
from matplotlib.colors import LinearSegmentedColormap   


twilight_lim_low = 0.2
twilight_lim_high = 0.8
twilight_lim_blu = 0.65
cred, cblu = twilight(twilight_lim_low), twilight(twilight_lim_blu)
restr_twilight_vals = twilight(
    np.linspace(twilight_lim_low, twilight_lim_high)
)
restr_twilight = LinearSegmentedColormap.from_list(
    "restr_twilight", restr_twilight_vals
)
cm.register_cmap(name="restr_twilight", cmap=restr_twilight)
cm.register_cmap(name='restr_twilight_r', cmap=restr_twilight.reversed())

PLT_SL2DSQ_SIDE1 = 7
PLT_SL2DSQ_SIDE2 = 7
PLT_SL2DSQ_UNIDS = 1.
PLT_SL2DSQ_LNEXTL = .75
PLT_SL2DSQ_KWNODE = dict(marker='o', ms=20, mec='k', mfc='w')
PLT_SL2DSQ_KWEXTL = dict(marker=None, ls=':', zorder=0)
PLT_SL2DSQ_KWLINE = dict(lw=3)
PLT_SL2DSQ_PEC = cblu
PLT_SL2DSQ_CPEC = cred
PLT_SL2DSQ_MODE = 'rand'
PLT_SL2DSQ_KWTXT = dict(fontsize=20, c='k', ha='center', va='center')
PLT_SL2DSQ_VDCPEC = [(0, 2), (2, 2), (3, 2), (1, 1), (1, 0)]
PLT_SL2DSQ_HDCPEC = [(3, 2), (3, 1), (1, 0), (2, 0)]
PLT_SL2DSQ_PFLIP = 0.2