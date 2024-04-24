import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from typing import Any

def adjust_to_even(x):
    """
    Rounds the input number to the nearest even integer.

    If the input is exactly halfway between two even numbers, it rounds up to the
    higher even number.

    Parameters
    ----------
    x : float
        The input number to be rounded.

    Returns
    -------
    int
        The nearest even integer to the input number.

    Examples
    --------
    >>> round_to_nearest_even(128 * np.sqrt(3))
    222
    >>> round_to_nearest_even(5.5)
    6
    >>> round_to_nearest_even(2.1)
    2
    """
    lower_even = int(x) - int(x) % 2
    upper_even = lower_even + 2
    return lower_even if x - lower_even < upper_even - x else upper_even

def triangular_lattice_graph_FastPatch(m: int, n: int, periodic: bool = False, with_positions: bool = True, create_using: Any = None) -> nx.Graph:
    """
    Generates a triangular lattice graph with optional periodic boundary conditions (PBC). 
    This function creates a graph representing a triangular lattice with `m` rows and `n` columns.
    Nodes in the lattice are connected in a manner that forms a pattern of equilateral triangles.
    When periodic boundary conditions are enabled, the lattice simulates a toroidal surface
    where edges wrap around the opposite side, creating a continuous pattern.
    
    Parameters
    ----------
    m : int
        The number of rows in the lattice.
    n : int
        The number of columns in the lattice.
    periodic : bool, optional (default=False)
        If True, applies periodic boundary conditions, simulating a toroidal topology.
        Requires `m >= 3` and `n >= 5`.
    with_positions : bool, optional (default=True)
        If True, calculates and stores the positions of each node in the node attribute 'pos',
        arranging nodes in equilateral triangles.

    Returns
    -------
    NetworkX graph
        A graph object representing the m by n triangular lattice, optionally with PBC.

    Raises
    ------
    NetworkXError
        If periodic is True and m < 3 or n < 5.
    """
    from networkx import empty_graph, NetworkXError, set_node_attributes
    G = empty_graph(0, create_using)

    # if periodic:
    #     if n < 5 or m < 3:
    #         msg = f"m > 2 and n > 4 required for periodic. m={m}, n={n}"
    #         raise NetworkXError(msg)
    rows = range(m)
    cols = range(n)
    # identify boundary nodes if periodic
    if periodic:
        for i in range(n+1):
            for j in range(m+1):
                G.add_node((i % n, j % m))  # Add node with PBC
                # Add horizontal edges within the grid, with PBC for the last column
                if i < n or j % 2 == 0:  # For even rows, wrap horizontally
                    G.add_edge((i % n, j % m), ((i + 1) % n, j % m))
                # Add vertical and diagonal edges, with PBC for the last row
                if j < m:
                    G.add_edge((i % n, j % m), (i % n, (j + 1) % m))
                    if j % 2:  # Diagonal for even rows
                        G.add_edge((i % n, j % m), ((i + 1) % n, (j + 1) % m))
                    else:  # Diagonal for odd rows, wrapping if at the edge
                        G.add_edge((i % n, j % m), ((i - 1) % n, (j + 1) % m))
    else:
        # Make grid
        G.add_edges_from(((i, j), (i + 1, j)) for j in rows for i in cols[:n-1])
        G.add_edges_from(((i, j), (i, j + 1)) for j in rows[:m-1] for i in cols)
        # add diagonals
        G.add_edges_from(((i, j), (i + 1, j + 1)) for j in rows[1:m-1:2] for i in cols[:n-1])
        G.add_edges_from(((i + 1, j), (i, j + 1)) for j in rows[:m-1:2] for i in cols[:n-1])

    # Add position node attributes
    if with_positions:
        ii = (i for i in rows for j in cols)
        jj = (j for i in rows for j in cols)
        xx = (0.5 * (i % 2) + j for i in rows for j in cols)
        h = np.sqrt(3) / 2
        if periodic:
            yy = (h * i + 0.01 * j * j for i in rows for j in cols)
        else:
            yy = (h * i for i in rows for j in cols)
        pos = {(j, i): (x, y) for i, j, x, y in zip(ii, jj, xx, yy) if (j, i) in G}
        set_node_attributes(G, pos, "pos")
    return G



def hexagonal_lattice_graph_FastPatch(n: int, m: int, periodic: bool = False, with_positions: bool = True, create_using: Any = None) -> nx.Graph:
    """
    Generate a hexagonal lattice graph with optional periodic boundary conditions.
    
    Parameters:
    m : int
        Number of rows of hexagons.
    n : int
        Number of columns of hexagons.
    periodic : bool, optional
        If True, applies periodic boundary conditions to create a toroidal topology.
    with_positions : bool, optional
        If True, calculates and assigns positions to nodes to represent a hexagonal lattice.

    Returns:
    nx.Graph
        A hexagonal lattice graph.
    """
    from networkx import set_node_attributes, empty_graph, NetworkXError
    G = empty_graph(0, create_using)
    if m == 0 or n == 0:
        return G
    # if periodic and (n % 2 == 1 or m < 2 or n < 2):
    #     msg = "periodic hexagonal lattice needs m > 1, n > 1 and even n"
    #     raise NetworkXError(msg)
    rows = range(m)
    cols = range(n)
    # Adjusted for clarity: In a hexagonal lattice, 'm' rows and 'n' columns describe a different layout.
    col_edges = (((i, j), (i, j + 1)) for i in cols for j in rows[:m-1])
    row_edges = (((i, j), (i + 1, j)) for i in cols[:n-1] for j in rows if i % 2 == j % 2)
    G.add_edges_from(col_edges)
    G.add_edges_from(row_edges)
    if periodic:
        more_row_edges = (((cols[0], j), (cols[n-1], j)) for j in rows if j % 2)
        more_col_edges = (((i, rows[0]), (i, rows[m-1])) for i in cols if i+1 % 2)
        G.add_edges_from(more_row_edges)
        G.add_edges_from(more_col_edges)
        
    # Optionally assign positions for visualization
    if with_positions:

        # calc position in embedded space
        ii = (i for i in cols for j in rows)
        jj = (j for i in cols for j in rows)
        xx = (0.5 + i + i // 2 + (j % 2) * ((i % 2) - 0.5) for i in cols for j in rows)
        h = np.sqrt(3) / 2
        if periodic:
            yy = (h * j + 0.01 * i * i for i in cols for j in rows)
        else:
            yy = (h * j for i in cols for j in rows)
        # exclude nodes not in G
        pos = {(i, j): (x, y) for i, j, x, y in zip(ii, jj, xx, yy) if (i, j) in G}
        set_node_attributes(G, pos, "pos")

    return G

side1 = 10
tri = triangular_lattice_graph_FastPatch(side1, side1, periodic=False, with_positions=True)

fig, ax = plt.subplots(1, 1, figsize=(6, 6))
ax.set_aspect('equal')
nx.draw(tri, ax=ax, pos=nx.get_node_attributes(tri, 'pos'), node_size=10, with_labels=False)
plt.show()

side2 = 10
side1 = adjust_to_even(side2/np.sqrt(3))
hex = hexagonal_lattice_graph_FastPatch(side1, side2, periodic=False, with_positions=True)

fig, ax = plt.subplots(1, 1, figsize=(6, 6))
ax.set_aspect('equal')
nx.draw(hex, ax=ax, pos=nx.get_node_attributes(hex, 'pos'), node_size=10, with_labels=False)
plt.show()

