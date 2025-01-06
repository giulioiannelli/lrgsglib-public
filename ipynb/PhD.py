#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
get_ipython().run_line_magic('matplotlib', 'inline')
#
from lrgsglib.core import *
#
move_to_rootf(print_tf=True)
plt.style.use('ipynb/nb_plotsheet.mplstyle')
#
path_plot = PATHPLOT / Path('PhD')


# ### Anderson Localization

# In[2]:


side = 128
geo = 'squared'
pflip_l = [0.0001, 0.001, 0.01, 0.1]


# In[5]:


from matplotlib.colors import SymLogNorm
norm = SymLogNorm(linthresh=1e-3, linscale=1, base=10)

fig, ax = plt.subplots(ncols=len(pflip_l), figsize=(len(pflip_l)*10, 10))
flatax = ax.flatten()[::-1]

for i,pflip in enumerate(pflip_l[::-1]):
    l = Lattice2D(side, geo=geo, pflip=pflip, seed=3)
    l.flip_random_fract_edges()
    l.compute_k_eigvV()
    l.eigV[0] = flip_to_positive_majority_adapted(l.eigV[0])
    # l.eigv[0] /= np.min(np.abs(l.eigv[0]))
    flatax[i].axis('off')
    div = make_axes_locatable(flatax[i])
    cax = div.append_axes("right", "5%", "5%")
    im = flatax[i].imshow(l.eigV[0].reshape(*l.syshape), cmap='plasma')
    cb = fig.colorbar(im, cax=cax)
fig.tight_layout()
fig.savefig(path_plot / Path('anderson_pflip.pdf'), bbox_inches='tight', transparent=True)


# In[82]:


l.seed


# In[68]:


fig, ax = plt.subplots(ncols=len(pflip_l), figsize=(6*10, 10))
flatax = ax.flatten()

from matplotlib.ticker import LogLocator


for i,pflip in enumerate(pflip_l):
    l = Lattice2D(side, geo=geo, pflip=pflip)
    l.flip_random_fract_edges()
    l.compute_k_eigvV()
    l.eigV[0] = flip_to_positive_majority_adapted(l.eigV[0])
    flatax[i].axis('off')
    div = make_axes_locatable(flatax[i])
    cax = div.append_axes("right", "5%", "5%")
    im = flatax[i].imshow(l.eigV[0].reshape(*l.syshape), cmap=restr_twilight, norm='log')
    cb = fig.colorbar(im, cax=cax)
    # cb.ax.yaxis.set_major_formatter()
    minor_locator = LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10)  # Minor ticks
    cb.ax.yaxis.set_minor_locator(minor_locator)
    cb.ax.minorticks_on()
fig.tight_layout()


# In[61]:


fig, ax = plt.subplots(ncols=len(pflip_l), figsize=(3*10, 10))
flatax = ax.flatten()

from matplotlib.ticker import FuncFormatter
from matplotlib.colors import SymLogNorm

def log_scientific_formatter(x, pos):
    if x > 0:
        log_value = int(np.log10(x))  # Compute the base-10 logarithm
        return f"$10^{{{log_value:d}}}$"  # Format in scientific notation
    elif x == 0:
        return "0"  # Special case for zero
    else:
        return f"$-10^{{{int(np.log10(-x)):d}}}$"  # For negative values

for i,pflip in enumerate(pflip_l):
    l = Lattice2D(side, geo=geo, pflip=pflip)
    l.flip_random_fract_edges()
    l.compute_k_eigvV()
    l.eigV[0] = flip_to_positive_majority_adapted(l.eigV[0])
    flatax[i].axis('off')
    div = make_axes_locatable(flatax[i])
    cax = div.append_axes("right", "5%", "5%")
    norm = SymLogNorm(linthresh=np.min(np.abs(l.eigV[0])), vmin=np.min(l.eigV[0]), vmax=np.max(l.eigV[0]))
    im = flatax[i].imshow(l.eigV[0].reshape(*l.syshape), cmap=restr_twilight, norm=norm)
    cb = fig.colorbar(im, cax=cax)
    formatter = FuncFormatter(log_scientific_formatter)
    cb.ax.yaxis.set_major_formatter(formatter)
fig.tight_layout()


# In[ ]:




