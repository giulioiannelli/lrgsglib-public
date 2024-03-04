import random
#
import numpy as np
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.colors as mplc
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
#
from cycler import cycler
from matplotlib import gridspec
from matplotlib.axes import Axes
from matplotlib.cm import hsv
from matplotlib.colors import Colormap, ListedColormap, LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle, Rectangle, Ellipse, PathPatch
from matplotlib.path import Path
from matplotlib.text import Text
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#
from .const import *







def set_alpha_torgb(rgbacol, alpha=0.5):
    """
    Sets the alpha (transparency) channel of an RGBA color tuple and returns a new RGBA color tuple.

    Parameters:
        rgbacol (tuple): A tuple representing an RGBA color in the format (R, G, B, A), where R, G, B, and A
                         are integers between 0 and 255.
        alpha (float, optional): The alpha (transparency) value to set for the RGBA color. Default is 0.5.

    Returns:
        tuple: A new RGBA color tuple in the format (R, G, B, alpha), where R, G, B are the same as in the
               input tuple, and alpha is the specified transparency value.

    Example:
        >>> set_alpha_torgb((255, 0, 0, 255), 0.2)
        (255, 0, 0, 0.2)

    Note:
        - This function creates a new RGBA color tuple with the specified alpha value while preserving the
          original RGB color components.
        - The alpha value should be a float between 0.0 (fully transparent) and 1.0 (fully opaque).
    """
    return (rgbacol[0], rgbacol[1], rgbacol[2], alpha)


def create_custom_colormap(c1="#0000ff", c2="#fc0303", nc: int = 0):
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
    nocol = dict(N=nc) if nc else dict()
        
    cmap = LinearSegmentedColormap.from_list("custom_colormap", colors, **nocol)
    return cmap


def generate_maxpercdiff_colormap(
    number_of_distinct_colors: int = 80, number_of_shades: int = 7
):
    """
    Generates a perceptually distinct colormap.

    Parameters:
    -----------
    number_of_distinct_colors : int, optional
        The number of distinct colors in the colormap. Default is 80.

    Returns:
    --------
    ListedColormap
        A ListedColormap object representing the generated colormap.

    Notes:
    ------
    This function generates a perceptually distinct colormap using a saw-tooth pattern in the HSV color space.
    The number of distinct colors can be customized by providing the `number_of_distinct_colors` parameter.

    Reference:
    -----------
    - Based on the "saw-like" pattern technique for generating distinct colors.
    - HSV colormap is used to create cyclic color variations.

    Example:
    --------
    >>> colormap = generate_colormap(100)
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> data = np.random.rand(10, 10)
    >>> plt.imshow(data, cmap=colormap)
    >>> plt.colorbar()
    >>> plt.show()
    """

    no_distinct_colors_wmultshades = int(
        np.ceil(number_of_distinct_colors / number_of_shades) * number_of_shades
    )

    # Create an array with uniformly drawn floats taken from <0, 1) partition
    linearly_distributed_nums = (
        np.arange(no_distinct_colors_wmultshades)
        / no_distinct_colors_wmultshades
    )

    # We are going to reorganize monotonically growing numbers in such a way that there will be a single array with a saw-like pattern
    #     but each sawtooth is slightly higher than the one before
    # First divide linearly_distributed_nums into number_of_shades sub-arrays containing linearly distributed numbers
    arr_by_shade_rows = linearly_distributed_nums.reshape(
        number_of_shades, no_distinct_colors_wmultshades // number_of_shades
    )

    # Transpose the above matrix (columns become rows) - as a result, each row contains a sawtooth with values slightly higher than the row above
    arr_by_shade_columns = arr_by_shade_rows.T

    # Keep the number of sawtooths for later
    number_of_partitions = arr_by_shade_columns.shape[0]

    # Flatten the above matrix - join each row into a single array
    nums_distributed_like_rising_saw = arr_by_shade_columns.reshape(-1)

    # HSV colormap is cyclic (https://matplotlib.org/tutorials/colors/colormaps.html#cyclic), we'll use this property
    initial_cm = hsv(nums_distributed_like_rising_saw)

    lower_partitions_half = number_of_partitions // 2
    upper_partitions_half = number_of_partitions - lower_partitions_half

    # Modify the lower half in such a way that colors towards the beginning of the partition are darker
    # First colors are affected more, colors closer to the middle are affected less
    lower_half = lower_partitions_half * number_of_shades
    for i in range(3):
        initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8 / lower_half)

    # Modify the second half in such a way that colors towards the end of the partition are less intense and brighter
    # Colors closer to the middle are affected less, colors closer to the end are affected more
    for i in range(3):
        for j in range(upper_partitions_half):
            modifier = (
                np.ones(number_of_shades)
                - initial_cm[
                    lower_half
                    + j * number_of_shades : lower_half
                    + (j + 1) * number_of_shades,
                    i,
                ]
            )
            modifier = j * modifier / upper_partitions_half
            initial_cm[
                lower_half
                + j * number_of_shades : lower_half
                + (j + 1) * number_of_shades,
                i,
            ] += modifier

    return ListedColormap(initial_cm)


def imshow_colorbar_caxdivider(
    mappable, ax, position="right", size="5%", pad=0.05
) -> Union[AxesDivider, Axes, Any]:
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
    list_caxdivider = ["position", "size", "pad"]
    kwargs_caxdivider = dict(zip(list_caxdivider, [position, size, pad]))
    #
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(**kwargs_caxdivider)
    #
    list_colorbar = ["mappable", "cax"]
    kwargs_colorbar = dict(zip(list_colorbar, [mappable, cax]))
    #
    clb = plt.colorbar(**kwargs_colorbar)
    return divider, cax, clb

def get_complementary_color(color_name):
    """
    Calculates the complementary color for a given color.

    The function takes a color name or hex value and returns the complementary color in hex format.
    It uses the matplotlib library for color conversions and calculations.

    Parameters:
    -----------
    color_name : str
        The name of the color or its hex value. This can be a common color name (like 'red')
        or a hex code (like '#FF0000'). If a common color name is provided, the function first
        converts it to its hex equivalent.

    Returns:
    --------
    str
        The hex code of the complementary color. The complementary color is calculated by inverting
        the RGB values of the original color. For instance, the complementary color of 'red' ('#FF0000')
        would be cyan ('#00FFFF').

    Example:
    --------
    comp_color = get_complementary_color('red')  # Returns '#00FFFF', which is cyan

    Notes:
    ------
    - The function assumes that the input color is either a valid named color in the matplotlib
    color maps or a valid hex color code.
    - If the color name is not recognized in the matplotlib color maps, it's treated as a hex code.
    - The RGB values are assumed to be in the range [0, 1], as per matplotlib's color representation.
    """

    from matplotlib import colors as mplcol
    if color_name in plt.colormaps():
        hex_color = mplcol.to_hex(color_name)
    else:
        hex_color = color_name

    # Convert hex to RGB
    rgb = mplcol.to_rgb(hex_color)

    # Calculate complementary color
    comp_rgb = (1 - rgb[0], 1 - rgb[1], 1 - rgb[2])

    # Convert back to hex
    comp_hex = mplcol.to_hex(comp_rgb)

    return comp_hex

def plot_square_lattice(
    ax,
    size: int = 7,
    kwargs_nodes: dict = dict(marker="o", ms=20, mec="k", mfc="w"),
    kwargs_extl: dict = dict(ls=":"),
    etxl_len: float = 0.75,
    kwargs_lines: dict = dict(lw=3),
    pec="blue",
    cpec="red",
):
    """
    Plots a 2D square lattice with customizable nodes, lines, and external lines,
    each with their own styling parameters.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object on which the lattice will be plotted.
    size : int, optional
        The size of the lattice (number of nodes in one dimension). Default is 7.
    kwargs_nodes : dict, optional
        Keyword arguments for styling the nodes in the plot.
        Includes 'marker' (shape of the marker), 'ms' (marker size),
        'mec' (marker edge color), and 'mfc' (marker face color). Default is a white filled circle.
    kwargs_extl : dict, optional
        Keyword arguments for styling the external (boundary) lines of the lattice.
        Includes 'ls' (line style). Default is dotted lines.
    extl_len : float, optional
        Length of the external lines extending from the lattice boundaries. Default is 0.75.
    kwargs_lines : dict, optional
        Keyword arguments for styling the internal lines of the lattice.
        Includes 'lw' (line width). Default is a line width of 3.
    pec : str, optional
        Primary edge color for the lattice lines. Default is 'blue'.
    cpec : str, optional
        Complementary edge color for the lattice lines. Default is 'red'.

    Notes:
    ------
    The function randomly assigns the primary or complementary edge color to each line
    in the lattice. This includes both the internal lines and the external boundary lines.
    The nodes are plotted over the lines for a clear visualization.

    Example:
    --------
    fig, ax = plt.subplots()
    plot_square_lattice(ax, size=5, pec='green', cpec='magenta')
    plt.show()

    This will plot a 5x5 lattice with green and magenta lines, default node style, and default external line style.
    """
    x, y = np.meshgrid(range(size), range(size))
    # Plot each point in the lattice
    for i in range(size):
        for j in range(size):
            kwargs_lines["color"] = random.choice([pec, cpec])
            import matplotlib as mpl
            if kwargs_lines["color"] == pec:
                if i < size - 1:
                    ax.plot(
                        [x[i, j], x[i + 1, j]],
                        [y[i, j], y[i + 1, j]],
                        zorder=1,
                        **kwargs_lines,
                    )  # Vertical
                if j < size - 1:
                    ax.plot(
                        [x[i, j], x[i, j + 1]],
                        [y[i, j], y[i, j + 1]],
                        zorder=1,
                        **kwargs_lines,
                    )  # Horizontal
            else:
                with mpl.rc_context({'path.sketch': (5, 15, 1)}):
                    if i < size - 1:
                        ax.plot(
                            [x[i, j], x[i + 1, j]],
                            [y[i, j], y[i + 1, j]],
                            zorder=1,
                            **kwargs_lines,
                        )  # Vertical
                    if j < size - 1:
                        ax.plot(
                            [x[i, j], x[i, j + 1]],
                            [y[i, j], y[i, j + 1]],
                            zorder=1,
                            **kwargs_lines,
                        )  # Horizontal
            ax.plot(x[i, j], y[i, j], zorder=2, **kwargs_nodes)  # Nodes

    # Adding dashed lines on the boundaries
    for i in range(size):
        kwargs_extl["color"] = random.choice([pec, cpec])
        # Left and right boundaries
        ax.plot(
            [x[i, 0], x[i, 0] - etxl_len],
            [y[i, 0], y[i, 0]],
            zorder=0,
            **kwargs_extl
        )
        ax.plot(
            [x[i, -1], x[i, -1] + etxl_len],
            [y[i, -1], y[i, -1]],
            zorder=0,
            **kwargs_extl
        )

        # Top and bottom boundaries
        ax.plot(
            [x[0, i], x[0, i]],
            [y[0, i], y[0, i] - etxl_len],
            zorder=0,
            **kwargs_extl
        )
        ax.plot(
            [x[-1, i], x[-1, i]],
            [y[-1, i], y[-1, i] + etxl_len],
            zorder=0,
            **kwargs_extl
        )

    # Remove axes
    ax.axis("off")


def pos_evolving_grid_rw(n_steps: int, p: float, init: str = 'random'):
    pos_state = np.zeros((n_steps, 3))  # Now also includes state in the last dimension]
    if init == 'random':
        pos_state[0, 2] = np.random.choice([-1, +1])  # Initial state
    elif init == 'fixed':
        pos_state[0, 2] = +1
    global edge_signs 
    edge_signs = {}

    # Function to update position and state
    def update_position_state(current_pos_state, direction, p):
        x, y, state = current_pos_state
        key = (x, y, direction)
        if key not in edge_signs:
            edge_signs[key] = -1 if np.random.rand() < p else 1
        sign = edge_signs[key]

        # Update position based on direction
        if direction == 0:  # Move up
            y += 1
        elif direction == 1:  # Move down
            y -= 1
        elif direction == 2:  # Move left
            x -= 1
        elif direction == 3:  # Move right
            x += 1
        
        # Flip state if necessary
        if sign == -1:
            state *= -1

        return np.array([x, y, state])

    # Generate the walk
    for i in range(1, n_steps):
        direction = np.random.randint(0, 4)  # Choose direction
        pos_state[i] = update_position_state(pos_state[i-1], direction, p)
    return pos_state

def average_evolving_rw(replica: int = 10**2, n_steps: int = 10**4, p=0.3, init: str = 'random'):
    cumulative_walk = {}
    for _ in range(replica):
        walk = pos_evolving_grid_rw(n_steps, p,  init = 'fixed')
        for x, y, value in walk:
            key = (int(x), int(y))
            if key in cumulative_walk:
                cumulative_walk[key] += value
            else:
                cumulative_walk[key] = value

    # Find extents for the 2D array
    min_x = min(key[0] for key in cumulative_walk)
    max_x = max(key[0] for key in cumulative_walk)
    min_y = min(key[1] for key in cumulative_walk)
    max_y = max(key[1] for key in cumulative_walk)

    # Create and fill the array
    array_shape = (int(max_y - min_y + 1), int(max_x - min_x + 1))
    walk_array = np.zeros(array_shape)
    for (x, y), value in cumulative_walk.items():
        walk_array[y - min_y, x - min_x] = value / replica  # Average the states

    return walk_array

def plot_evolving_grid_rw(n_steps: int = 10**4, p: float = 0.103, ax: Axes = plt.Axes, col1: Any = 'red', col2: Any = 'blue', init: str = 'random'):
    # Initialize a 3D array: (n_steps, 2 positions, 1 state)
    colors = {1: col1, -1: col2}
    pos_state = pos_evolving_grid_rw(n_steps, p, init = init)
    # Plot the steps with appropriate color based on the state
    for i in range(1, n_steps):
        ax.plot(pos_state[i-1:i+1, 0], pos_state[i-1:i+1, 1], color=colors[pos_state[i, 2]], linewidth=2)




def perform_random_walks(G, steps, N):
    """
    Perform N random walks on graph G, each with a given number of steps,
    and update node states based on edge weights.
    """
    node_states = {node: 0 for node in G.nodes()}  # Initial node states
    
    for _ in range(N):
        current_state = np.random.choice([-1, 1])
        current_node = list(G.nodes())[np.random.randint(len(G))]  # Start at the origin for each walk
        node_states[current_node] += current_state
        for _ in range(steps):
            neighbors = list(G.neighbors(current_node))
            next_node = neighbors[np.random.randint(len(neighbors))]
            edge_weight = G.edges[current_node, next_node]['weight']
            # Update the state of the next node based on edge weight
            current_state = current_state*edge_weight
            node_states[next_node] += current_state
            current_node = next_node
            
    return node_states

def visualize_node_states(node_states, n, m, ax: Axes = plt.Axes):
    """
    Visualize the final states of nodes as a heatmap.
    """
    # Convert node states to a 2D array
    state_array = np.zeros((n, m))
    for (x, y), state in node_states.items():
        state_array[x, y] = state
    
    ax.imshow(state_array)
    ax.axis('off')  # Hide the axes


def set_ax_ratio_1_withlim(ax):
    # Calculate the ranges and centers
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    x_range = x_max - x_min
    y_range = y_max - y_min

    # Determine the largest range
    max_range = max(x_range, y_range) / 2

    # Set the new limits to ensure the plot is square and centered
    ax.set_xlim(x_center - max_range, x_center + max_range)
    ax.set_ylim(y_center - max_range, y_center + max_range)