#
from ..shared import *
from .const import *
from .const_plotlib import *
from .funcs import *
from .errwar import *
#
def set_new_lower_ybound(ax, new_lower_bound):
    """
    Sets a new lower bound for the y-axis of a matplotlib axes object,
    while keeping the upper bound unchanged.
    
    Parameters:
    ----------
    ax : matplotlib.axes.Axes
        The axes object to modify the y-axis limits of.
    new_lower_bound : float
        The new lower limit for the y-axis.
    """
    if new_lower_bound is not None:
        current_ylim = ax.get_ylim()  # Get current y-axis limits
        ax.set_ylim(new_lower_bound, current_ylim[1])  # Update the lower limit

def get_opposite_color(color: ColorType, output_format: str = 'hex'):
    """
    Calculate the opposite color for a given color in various formats.

    Parameters:
    -----------
    color : ColorType
        The input color which can be in any of the following formats:
        - Hexadecimal string (e.g., '#6496c8')
        - Named color (e.g., 'red')
        - RGB tuple with values in the range [0, 1] (e.g., (0.1, 0.5, 0.8))
        - RGB tuple with values in the range [0, 255] (e.g., (100, 150, 200))
        - RGBA tuple with values in the range [0, 1] or [0, 255] (e.g., (0.1, 0.5, 0.8, 1.0) or (100, 150, 200, 255))
    output_format : str, optional
        The desired output format of the opposite color. Options are:
        - 'hex': Returns the color as a hexadecimal string (default)
        - 'rgb': Returns the color as an RGB tuple with values in the range [0, 1]

    Returns:
    --------
    str or tuple
        The opposite color represented in the specified format.

    Example:
    --------
    >>> opposite_color((0.1, 0.5, 0.8))
    '#e66a33'
    >>> opposite_color('#6496c8', output_format='rgb')
    (0.6078431372549019, 0.4117647058823529, 0.21568627450980393)
    >>> opposite_color('red')
    '#00ffff'
    """
    # Convert the input color to an RGB tuple (values between 0 and 1)
    try:
        if isinstance(color, tuple):
            if len(color) == 3:
                if max(color) > 1:  # Assuming the values are in [0, 255]
                    rgb = tuple(c / 255 for c in color)
                else:
                    rgb = color
            elif len(color) == 4:  # Ignore alpha channel for the opposite color
                if max(color) > 1:  # Assuming the values are in [0, 255]
                    rgb = tuple(c / 255 for c in color[:3])
                else:
                    rgb = color[:3]
            else:
                raise ValueError("Unsupported tuple length for color. Must be 3 (RGB) or 4 (RGBA).")
        else:
            rgb = to_rgb(color)
    except ValueError:
        raise ValueError("Unsupported color format. Please provide a valid color.")

    # Convert RGB to the opposite color
    rgb_array = np.array(rgb)
    opposite_rgb = 1 - rgb_array  # Find the opposite for each channel

    # Return the opposite color in the desired output format
    if output_format == 'hex':
        return to_hex(opposite_rgb)
    elif output_format == 'rgb':
        return tuple(opposite_rgb)
    else:
        raise ValueError("Unsupported output format. Please use 'hex' or 'rgb'.")

def set_alpha_torgb(rgbacol, alpha=0.5):
    """
    Sets the alpha (transparency) channel of an RGBA color tuple and returns a 
    new RGBA color tuple.

    Parameters:
        rgbacol (tuple): A tuple representing an RGBA color in the format 
            (R, G, B, A), where R, G, B, and A are integers between 0 and 255.
        alpha (float, optional): The alpha (transparency) value to set for the 
            RGBA color. Default is 0.5.

    Returns:
        tuple: A new RGBA color tuple in the format (R, G, B, alpha), where 
            R, G, B are the same as in the input tuple, and alpha is the 
            specified transparency value.

    Example:
        >>> set_alpha_torgb((255, 0, 0, 255), 0.2)
        (255, 0, 0, 0.2)

    Note:
        - This function creates a new RGBA color tuple with the specified alpha 
            value while preserving the original RGB color components.
        - The alpha value should be a float between 0.0 (fully transparent) 
            and 1.0 (fully opaque).
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
    - The function creates a custom colormap by defining a transition from 
        one color to another.
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
    This function generates a perceptually distinct colormap using a saw-tooth 
        pattern in the HSV color space. The number of distinct colors can be ]
        customized by providing the `number_of_distinct_colors` parameter.

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
        Default is 'right', which means the colorbar will be placed to the 
        right of 'ax'. Other options include 'top', 'bottom', 'left', 
        and 'right'.

    size : str or float, optional
        The size of the colorbar relative to 'ax'. Can be specified as a 
        percentage (e.g., "5%") or as an absolute size in points (e.g., 100). 
        Default is "5%".

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

    The function takes a color name or hex value and returns the complementary 
    color in hex format. It uses the matplotlib library for color conversions 
    and calculations.

    Parameters:
    -----------
    color_name : str
        The name of the color or its hex value. This can be a common color 
        name (like 'red') or a hex code (like '#FF0000'). If a common color 
        name is provided, the function first converts it to its hex equivalent.

    Returns:
    --------
    str
        The hex code of the complementary color. The complementary color is 
        calculated by inverting the RGB values of the original color. For 
        instance, the complementary color of 'red' ('#FF0000')
        would be cyan ('#00FFFF').

    Example:
    --------
    comp_color = get_complementary_color('red')  # Returns '#00FFFF', which is cyan

    Notes:
    ------
    - The function assumes that the input color is either a valid named color
        in the matplotlib color maps or a valid hex color code.
    - If the color name is not recognized in the matplotlib color maps, it's 
        treated as a hex code.
    - The RGB values are assumed to be in the range [0, 1], as per 
        matplotlib's color representation.
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

def scheme_Lattice2DSquared(
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
    scheme_Lattice2DSquared(ax, size=5, pec='green', cpec='magenta')
    plt.show()

    This will plot a 5x5 lattice with green and magenta lines, default node style, and default external line style.
    """
    x, y = np.meshgrid(range(size), range(size))
    # Plot each point in the lattice
    for i in range(size):
        for j in range(size):
            kwargs_lines["color"] = random.choice([pec, cpec])
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
                with rc_context({'path.sketch': (5, 15, 1)}):
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

def plot_line_with_style(ax, x_coords, y_coords, color, kwargs_lines, cpec):
    """
    Helper function to plot a line with a specific color and optional style.
    """
    kwargs_lines["color"] = color
    if color == cpec:
        with rc_context({'path.sketch': (5, 15, 1)}):
            ax.plot(x_coords, y_coords, zorder=1, **kwargs_lines)
    else:
        ax.plot(x_coords, y_coords, zorder=1, **kwargs_lines)
#
def scheme_Lattice2DSquared(
        ax: Axes, 
        side1: int = PLT_SL2DSQ_SIDE1,
        side2: int = PLT_SL2DSQ_SIDE2,
        mode: str = PLT_SL2DSQ_MODE,
        kwargNodes: dict = PLT_SL2DSQ_KWNODE,
        kwargsExtl: dict = PLT_SL2DSQ_KWEXTL,
        lenExtl: float = PLT_SL2DSQ_LNEXTL,
        kwargsLines: dict = PLT_SL2DSQ_KWLINE,
        pec: ColorType = PLT_SL2DSQ_PEC,
        cpec: ColorType = PLT_SL2DSQ_CPEC,
        pflip: float = PLT_SL2DSQ_PFLIP,
        kwargsTxt: dict = PLT_SL2DSQ_KWTXT
        ) -> None:
    """
        Function to plot a square lattice where the color and style of each link
        can be controlled individually. This includes the ability to create defects
        or randomly alter the appearance of links and nodes within the lattice.

        Parameters
        ----------
        ax : Axes
            The matplotlib axes to plot on.
        side1 : int, optional
            Number of nodes on one side of the square lattice.
        side2 : int, optional
            Number of nodes on the other side of the square lattice.
        mode : str, optional
            The mode of coloring the lattice links. Can be 'rand' for random 
            colors or 'defects' to specify certain links with different colors.
        kwargNodes : dict, optional
            Keyword arguments for plotting the nodes.
        kwargsExtl : dict, optional
            Keyword arguments for plotting the extended links.
        lenExtl : float, optional
            Length of the extended links.
        kwargsLines : dict, optional
            Keyword arguments for plotting the lattice links.
        pec : ColorType, optional
            Primary color for the lattice links.
        cpec : ColorType, optional
            Color for the specified or random links based on the mode.
        pflip : float, optional
            Probability of flipping the color of a link in 'rand' mode.
        kwargsTxt : dict, optional
            Keyword arguments for plotting text labels on the lattice.

        Returns
        -------
        None
            This function does not return a value. It modifies the given Axes object
            in-place, adding a visual representation of a 2D squared lattice.

        Notes
        -----
        The 'defects' mode allows for specifying particular links to have a different
        appearance (e.g., color) to represent defects or special cases in the lattice.
        The 'rand' mode uses the `pflip` parameter to randomly change the appearance of
        links, simulating randomness or disorder within the lattice.

        This function is designed to be flexible, with many aspects of the visualization
        customizable through keyword arguments. This allows for a wide range of visual
        styles and representations to suit different requirements or preferences.
    """
    if side1 != PLT_SL2DSQ_SIDE1 and side2 == PLT_SL2DSQ_SIDE1:
        side2 = side1
    if kwargsExtl == PLT_SL2DSQ_KWEXTL:
        kwargsExtl.update(dict(c=pec))
    x, y = np.meshgrid(np.arange(side1), np.arange(side2), indexing='ij')
    #
    if mode == 'defects':
        def determine_line_color(i, j, direction):
            if direction == "vertical":
                if (i, j) in PLT_SL2DSQ_VDCPEC:
                    return cpec
                return pec
            # Conditions for horizontal lines
            elif direction == "horizontal":
                if (i, j) in PLT_SL2DSQ_HDCPEC:
                    return cpec
                return pec
            else: 
                return pec
    elif mode == 'rand':
        def determine_line_color(*args):
            return random.choices([pec, cpec], [1-pflip, pflip], k=1)[0]
    #
    for i in range(side1):
        for j in range(side2):
            # Vertical lines
            if i < side1 - 1:
                xPt = [x[i, j], x[i+1, j]]
                yPt = [y[i, j], y[i+1, j]]
                lc = determine_line_color(i, j, "vertical")
                plot_line_with_style(ax, xPt, yPt, lc, kwargsLines, cpec)
            # Horizontal lines
            if j < side2 - 1:
                xPt = [x[i, j], x[i, j+1]]
                yPt = [y[i, j], y[i, j+1]]
                lc = determine_line_color(i, j, "horizontal")
                plot_line_with_style(ax, xPt, yPt, lc, kwargsLines, cpec)
            # Nodes
            ax.plot(x[i, j], y[i, j], zorder=2, **kwargNodes)
    
    for j in range(side2):
        # Extend left from the left nodes
        ax.plot([x[0, j], x[0, j]- lenExtl], [y[0, j], y[0, j]], **kwargsExtl)
        # Extend right from the right nodes
        ax.plot([x[-1, j], x[-1, j]+ lenExtl], 
                [y[-1, j], y[-1, j]], **kwargsExtl)

    for i in range(side1):
        # Extend downwards from the bottom nodes
        ax.plot([x[i, 0], x[i, 0]], [y[i, 0], y[i, 0] - lenExtl], **kwargsExtl)
        # Extend upwards from the top nodes
        ax.plot([x[i, -1], x[i, -1]], 
                [y[i, -1], y[i, -1] + lenExtl], **kwargsExtl)


    if mode == 'defects':
        ax.text(x[0, 2] + .5, y[0, 2] + 0.3, rf'$S$', color=cpec, **kwargsTxt)
        ax.text(x[3, 2] + .5, y[3, 2] + 0.3, r"$X$", color=cpec, **kwargsTxt)
        ax.text(x[2, 0] + .3, y[2, 0] + 0.5, r"$Z$", color=cpec, **kwargsTxt)


def defects_on_lattice_plot(sizes, lattices, ax, direction: str = 'parallel', 
                            geometry: str = 'squared', cell: str = 'single', 
                            fit_mode: str = 'lmfit'):
    #
    from .funcs import flip_to_positive_majority
    from scipy.optimize import curve_fit
    #
    newLowerBound = None
    xShiftConst = 0
    kwlogfit = dict(marker='', c='red')# label=r'$a\log(x) + b$'
    if direction == 'perpendicular':
        ylabel = r'${\phi(\bar{x}_\perp,\, y)}/{\phi_{\min}}$'
        ax.set_xlabel(r'$y$')
    else:
        ylabel = r'${\phi(x,\, \bar{y}_\parallel)}/{\phi_{\min}}$'
        ax.set_xlabel(r'$x$')
    if geometry == 'squared':
        if cell == 'single':
            if direction == 'parallel':
                xShiftConst = 0
                newLowerBound = .5
                slice_cut = lambda side: np.s_[:, lattices[side].side2//2]
            else:
                xShiftConst = -.5
                newLowerBound = 0
                slice_cut = lambda side: np.s_[lattices[side].side1//2, :]
        if cell == 'singleZERR':
            if direction == 'parallel':
                xShiftConst = 0
                newLowerBound = -.5
                slice_cut = lambda side: np.s_[:, lattices[side].side2//2-1]
            else:
                xShiftConst = 1
                newLowerBound = -.5
                slice_cut = lambda side: np.s_[lattices[side].side1//2, :]

        if cell == 'singleXERR':
            if direction == 'parallel':
                xShiftConst = -.5
                slice_cut = lambda side: np.s_[:, lattices[side].side2//2-1]
            else:
                xShiftConst = +.5
                slice_cut = lambda side: np.s_[lattices[side].side1//2, :]
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV = np.max(eigV) - eigV
                eigV /= np.max(eigV)
                return eigV
        else:
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV = np.max(eigV) - eigV
                eigV /= np.max(eigV)
                return eigV
    elif geometry == 'triangular':
        xShiftConst = 0
        if cell != 'singleXERR':
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV /= np.min(eigV)
                return eigV
        else:
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV /= np.max(eigV)
                return eigV
        if direction == 'parallel':
            xShiftConst = -1
            if cell == 'singleZERR':
                xShiftConst = 0
            elif cell == 'singleXERR':
                xShiftConst = -.5
            slice_cut = lambda side: np.s_[lattices[side].side1//2, :]
            if cell != 'singleXERR':
                newLowerBound = .5 
        else:
            xShiftConst = -.5
            slice_cut = lambda side: np.s_[:, lattices[side].side2//2]

    ax.set_ylabel(ylabel, labelpad=10)
    ax.set_xscale('symlog')
    ax.set_yscale('log')
    #
    cmapv = restr_twilight(np.linspace(0, 1, len(sizes)))
    lista = []
    for side, c in zip(sizes[::-1], cmapv):
        kwdict = {"ls": '-',
                'c': c, 
                'marker': 'H', 
                'ms': 10, 
                'mfc': set_alpha_torgb(c, 0.75), 
                'label': fr"$N={side**2}$"} 
        eigV = unravel_1d_to_2d_nodemap(lattices[side].eigV[0], lattices[side].invnode_map, lattices[side].syshape)
        eigen_state = func(eigV)
        phi_plot = eigen_state[slice_cut(side)]
        if side == sizes[-1]:
            phi_plot0 = phi_plot
        # print(np.min(eigen_state), min(phi_plot), slice_cut(side))
        if geometry == 'squared':
            sideop = side
        elif geometry == 'triangular':
            if direction == 'parallel':
                sideop = lattices[side].side2 
            elif direction == 'perpendicular':
                sideop = lattices[side].side1
        x = np.rint(np.linspace(-sideop//2, sideop//2, num=sideop)+xShiftConst)
        ax.plot(x, phi_plot, **kwdict)
        ax.plot(np.argmax(np.abs(x))+2, phi_plot[np.argmax(np.abs(x))]+2, 'or')
        lista.append(x[np.argmax(np.abs(x))+2]+sideop//2)
        print(np.argmax(np.abs(x)), x[np.argmax(np.abs(x)):np.argmax(np.abs(x))+2], )
    #
    if fit_mode and cell != 'singleXERR':
        # idx = (x < -side//3) | (x > side//3)
        x_values = np.linspace(-sizes[-1]//2, sizes[-1]//2, num=sizes[-1])+xShiftConst
        y_values = phi_plot0
        
        x_plot = np.concatenate(
            [
                np.linspace(-sizes[-1]//2, -1, num=100),
                np.linspace(-1, 1, num=2000),
                np.linspace(1, sizes[-1]//2, num=100)
            ]
        )
        if fit_mode == 'lmfit':
            true_params = dict(a=1, b=2, c=0.5, d=0.1)
            # Create a Model with the function and fit to the data
            regressor = lmfit.Model(sym_log_func, nan_policy='omit')
            params = regressor.make_params(a=1, b=1.5, c=1, d=0)
            params['b'].set(min=1e-10)  # Prevent b from being zero or negative
            params['d'].set(min=-np.inf, max=np.inf) 
            result = regressor.fit(y_values, params, x=x_values)
            # y_fit = regressor.eval(params=result.params, x=x_plot)
            popt = np.array([result.params[key].value for key in result.params])
        elif fit_mode == 'scipy':
            weights = np.abs(x)
            popt, _ = curve_fit(sym_log_func, x, phi_plot0, sigma=1/weights, 
                                p0=[1, 1.5, 1, 0])
        y_fit = sym_log_func_unsafe(x_plot, *popt)
        ax.plot(x_plot, y_fit, **kwlogfit)
    set_new_lower_ybound(ax, newLowerBound)
    #
    kwvlines = {'ls':'--', 'lw':1, 'c':'k'}
    for i in [-1, +1, -lattices[sizes[0]].r_c, +lattices[sizes[0]].r_c]:
        ax.axvline(i, **kwvlines)
    return lista
#
def plot_log_distribution(data, fig_ax=None, binnum=20, 
                          xlabel="X-Axis (log scale)", ylabel="Frequency", 
                          title="Log-Scale Distribution", grid=True, 
                          log_scale=True, **kwargs):
    """
    Plots the distribution of the given data using a logarithmic scale on the x-axis.

    Parameters:
    -----------
    data : np.ndarray or list
        The data to be plotted in the distribution.

    fig_ax : tuple, optional
        Tuple of (fig, ax) to use an existing figure and axis. If None, a new 
        figure and axis are created. Default is None.

    binnum : int, optional
        Number of bins for the histogram. Default is 20.

    xlabel : str, optional
        Label for the x-axis. Default is "X-Axis (log scale)".

    ylabel : str, optional
        Label for the y-axis. Default is "Frequency".

    title : str, optional
        Title for the plot. Default is "Log-Scale Distribution".

    grid : bool, optional
        Whether to display a grid. Default is True.

    log_scale : bool, optional
        Whether to apply a logarithmic scale to the x-axis. Default is True.

    **kwargs : optional
        Additional keyword arguments for the `plt.bar` function to customize the
        plot appearance (e.g., `color`, `alpha`, etc.).

    Returns:
    --------
    None
        Displays the plot using the provided or created figure and axis.
    """
    # Logarithmic binning of data
    bin_centers, hist, bin_w = log_binning(data, binnum=binnum)
    
    # Check if custom fig and ax are provided, otherwise create them
    if fig_ax is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = fig_ax
    
    # Plot the histogram using a bar chart with customizable kwargs
    ax.bar(bin_centers, hist, width=bin_w, alpha=0.75, ec='black', **kwargs)
    
    # Apply logarithmic scale if specified
    if log_scale:
        ax.set_xscale('log')
    
    # Set labels, title, and grid as per the provided arguments
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    if grid:
        ax.grid(True, which="both", ls="--")  # Grid for both major and minor ticks
    
    # # Show the plot if no custom fig_ax is provided
    # if fig_ax is None:
    #     plt.show()
def set_color_cycle(
    arr: Sequence,
    ax: Optional[Axes] = None,
    my_cmap: Optional[Colormap] = None
) -> None:
    """
    Sets the color cycle of the given Axes based on the provided colormap and array length.

    Parameters
    ----------
    arr : Sequence
        An array-like object that determines the number of distinct colors needed. The length of `arr` 
        dictates how many colors will be generated from the colormap.
        
    ax : matplotlib.axes.Axes, optional
        The Matplotlib Axes object to which the color cycle will be applied. This is where 
        the color properties will be set, affecting subsequent plot elements added to this Axes.
        If not provided, the current axes (`plt.gca()`) will be used.
    
    my_cmap : matplotlib.colors.Colormap, optional
        A Matplotlib colormap instance used to map normalized values to colors. This colormap 
        defines the range and variation of colors in the cycle. If not provided, defaults to 
        the custom colormap `'restr_twilight'`.
        

    Returns
    -------
    None
        This function does not return any value. It modifies the `ax` object in place by setting 
        its color property cycle.

    Raises
    ------
    ValueError
        If `arr` is empty, as at least one color is required to set the color cycle.
    LookupError
        If the default colormap `'restr_twilight'` is not found and `my_cmap` is not provided.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from matplotlib import cm
    >>> # Assume 'restr_twilight' colormap is defined elsewhere
    >>> # If not, define it here or replace with an existing colormap
    >>> # Sample data
    >>> arr = [1, 2, 3, 4, 5]
    >>> x = np.linspace(0, 10, 100)
    >>> # Create a figure and axes
    >>> fig, ax = plt.subplots()
    >>> # Apply the custom color cycle with default colormap
    >>> set_color_cycle(arr, ax=ax)
    >>> # Plot multiple lines; each line will use the next color in the cycle
    >>> for i, val in enumerate(arr):
    >>>     y = np.sin(x + i)
    >>>     ax.plot(x, y, label=f'Line {i+1}')
    >>> # Add legend and show plot
    >>> ax.legend()
    >>> plt.show()

    Notes
    -----
    - Ensure that the length of `arr` corresponds to the number of distinct elements you plan to plot 
      to avoid color repetition.
    - The colormap (`my_cmap`) can be any Matplotlib colormap. You can create custom colormaps or use 
      predefined ones like `'viridis'`, `'plasma'`, `'inferno'`, `'magma'`, etc.
    - If using the default `'restr_twilight'` colormap, ensure it is registered with Matplotlib using 
      `plt.register_cmap()` before calling this function.
    """
    if not arr:
        raise ValueError("The input array 'arr' must contain at least one element to set the color cycle.")
    # Set default axes to current axes if not provided
    if ax is None:
        ax = plt.gca()
    # Set default colormap to 'restr_twilight' if not provided
    if my_cmap is None:
        try:
            my_cmap = plt.get_cmap('restr_twilight')
        except ValueError:
            raise LookupError("The default colormap 'restr_twilight' is not found. Please provide a valid colormap.")
    # Generate evenly spaced values between 0 and 1 based on the length of arr
    normalized_values = np.linspace(0.0, 1.0, len(arr))
    # Map the normalized values to colors using the provided colormap
    colors = my_cmap(normalized_values)
    # Create a cycler with the generated colors
    color_cycle = cycler(color=colors)
    # Apply the color cycle to the Axes
    ax.set_prop_cycle(color_cycle)

def plot_honeycomb_grid(data: NDArray, fig: Any, ax: Any, triangle_size: float = 1.0) -> None:
    """
    Plot a triangular grid representing the honeycomb structure based on the provided 2D array data.

    Parameters
    ----------
    data : np.ndarray
        A 2D numpy array representing values at each point in the honeycomb structure.
    fig : matplotlib.figure.Figure
        A matplotlib figure object to plot on.
    ax : matplotlib.axes._axes.Axes
        A matplotlib axes object to plot on.
    triangle_size : float, optional
        The side length of the equilateral triangles in the grid. Default is 1.0.
    
    Returns
    -------
    None
    """
    import matplotlib.patches as patches
    # Determine number of rows and columns from data
    n_rows, n_cols = data.shape

    # Height of an equilateral triangle
    triangle_height = (np.sqrt(3) / 2) * triangle_size

    # Loop over each row and column to create triangular patches
    for row in range(n_rows):
        for col in range(n_cols):
            # Determine the x and y positions for the center of each triangle
            x_base = col * triangle_size / 2
            y_base = row * triangle_height

            # Define the vertices of the triangles (upward and downward)
            if (row + col) % 2 == 0:
                # Upward-pointing triangle
                vertices = [
                    (x_base, y_base + (triangle_height / 2)),
                    (x_base - (triangle_size / 2), y_base - (triangle_height / 2)),
                    (x_base + (triangle_size / 2), y_base - (triangle_height / 2)),
                ]
            else:
                # Downward-pointing triangle
                vertices = [
                    (x_base, y_base - (triangle_height / 2)),
                    (x_base - (triangle_size / 2), y_base + (triangle_height / 2)),
                    (x_base + (triangle_size / 2), y_base + (triangle_height / 2)),
                ]

            # Get color based on data
            color = credcblu(data[row, col])

            # Create and add the triangle patch to the axis
            triangle = patches.Polygon(vertices, facecolor=color)
            ax.add_patch(triangle)

    # Set the limits and aspect ratio to ensure proper display
    ax.set_xlim(-triangle_size / 2, (n_cols / 2 + 0.5) * triangle_size)
    ax.set_ylim(-triangle_height / 2, (n_rows + 0.5) * triangle_height)
    ax.set_aspect('equal')

    # Hide the axes for better visualization
    ax.axis('off')
