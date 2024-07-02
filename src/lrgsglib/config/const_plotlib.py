from ..shared import *
#
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.colors as mplc
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
#
from matplotlib import gridspec, rc_context
from matplotlib.axes import Axes
from matplotlib.cm import hsv, twilight, ScalarMappable
from matplotlib.colors import Colormap, ListedColormap, BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle, Rectangle, Ellipse, PathPatch
from matplotlib.path import Path
from matplotlib.text import Text
from matplotlib.ticker import ScalarFormatter, MultipleLocator
#
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from mpl_toolkits.mplot3d import Axes3D

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
credcblu = ListedColormap([cred, cblu])
cm.register_cmap(name="restr_twilight", cmap=restr_twilight)
cm.register_cmap(name='restr_twilight_r', cmap=restr_twilight.reversed())

red_blue = LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])
cm.register_cmap(name="red_blue", cmap=red_blue)
cm.register_cmap(name="red_blue_r", cmap=red_blue.reversed())

PLT_SL2DSQ_SIDE1 = 7
PLT_SL2DSQ_SIDE2 = 7
PLT_SL2DSQ_UNIDS = 1.
PLT_SL2DSQ_LNEXTL = .75
PLT_SL2DSQ_KWNODE = dict(marker='o', ms=20, mec='k', mfc='w')
PLT_SL2DSQ_KWEXTL = dict(marker='', ls=':', zorder=0)
PLT_SL2DSQ_KWLINE = dict(lw=3)
PLT_SL2DSQ_PEC = cblu
PLT_SL2DSQ_CPEC = cred
PLT_SL2DSQ_MODE = 'rand'
PLT_SL2DSQ_KWTXT = dict(fontsize=24, c='k', ha='center', va='center')
PLT_SL2DSQ_VDCPEC = [(0, 2), (2, 2), (3, 2), (1, 1), (1, 0)]
PLT_SL2DSQ_HDCPEC = [(3, 2), (3, 1), (1, 0), (2, 0)]
PLT_SL2DSQ_PFLIP = 0.2