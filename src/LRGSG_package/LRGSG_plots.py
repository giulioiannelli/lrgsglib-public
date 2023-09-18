import matplotlib.pyplot as plt
#
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, Normalize

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider

from typing import Any, Optional, Union

def create_custom_colormap(c1="#0000ff", c2="#fc0303"):
    """
    Create a custom colormap transitioning between two specified colors.

    Parameters:
    -----------
    c1 : str, optional
        The starting color in hexadecimal format (e.g., "#0000ff" for blue).
        Default is "#0000ff".

    c2 : str, optional
        The ending color in hexadecimal format (e.g., "#fc0303" for red).
        Default is "#fc0303".

    Returns:
    --------
    matplotlib.colors.LinearSegmentedColormap
        A custom colormap transitioning from 'c1' to 'c2'.

    Notes:
    ------
    - The function creates a custom colormap by defining a transition from one color to another.
    - It returns the colormap object that can be used for coloring plots.

    Example:
    --------
    custom_cmap = create_custom_colormap(c1="#00ff00", c2="#ff0000")
    # The result is a custom colormap transitioning from green to red.

    """
    colors = [c1, c2]  # Red to black
    cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)
    return cmap

def imshow_colorbar_caxdivider(mappable, ax, position="right", size="5%", 
                               pad=0.05) -> [AxesDivider, Axes, Any]:
    """
    Display colorbar in a specified position relative to a given axis.

    Parameters:
    -----------
    im : matplotlib.image.AxesImage
        The image to display.

    ax : matplotlib.axes.Axes
        The axis where the image will be displayed.

    position : str, optional
        The position of the colorbar relative to the 'ax' parameter.
        Default is 'right', which means the colorbar will be placed to the right of 'ax'.
        Other options include 'top', 'bottom', 'left', and 'right'.

    size : str or float, optional
        The size of the colorbar relative to 'ax'. Can be specified as a percentage (e.g., "5%")
        or as an absolute size in points (e.g., 100). Default is "5%".

    pad : float, optional
        The padding between the colorbar and the 'ax' in relative units.
        Default is 0.05, which corresponds to 5% of the axis size.

    Returns:
    --------
    div : mpl_toolkits.axes_grid1.mpl_axes.Axes
        The axes divider created for positioning the colorbar.

    cax : matplotlib.axes.Axes
        The colorbar axis where the colorbar is displayed.

    clb : matplotlib.colorbar.Colorbar
        The colorbar object created.

    Examples:
    ---------
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import numpy as np

    # Create a sample image and axis
    img = np.random.rand(10, 10)
    fig, ax = plt.subplots()

    # Display the image with a colorbar on the right side
    imshow_colorbar_caxdivider(ax.imshow(img), ax, position="right")

    plt.show()

    """
    list_caxdivider = ['position', 'size', 'pad']
    kwargs_caxdivider = dict(zip(list_caxdivider, [position, size, pad]))
    #
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(**kwargs_caxdivider)
    #
    list_colorbar = ['mappable', 'cax']
    kwargs_colorbar = dict(zip(list_colorbar, [mappable, cax]))
    #
    clb = plt.colorbar(**kwargs_colorbar)
    return divider, cax, clb