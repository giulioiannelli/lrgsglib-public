import scipy.sparse as scsp
#
from ..const import *
from ..utils import *
from networkx.drawing.layout import _process_params, rescale_layout
#
def get_kth_order_neighbours(G: nx.Graph, node: int, order: int = 1) -> list:
    """
    Returns the kth-order neighbors of a given node in a networkx graph.

    Parameters
    ----------
    G : Graph
        A NetworkX graph.
    node : int
        The node for which kth-order neighbors are to be found.
    order : int, optional
        The order of neighbors to be considered (default is 1).

    Returns
    -------
    kth_order_neighbors : list
        A list of kth-order neighbors of the given node.

    Notes
    -----
    - The graph G should be connected.
    - The function uses the single_source_shortest_path_length function of networkx to compute the shortest path lengths from the given node to all other nodes in the graph.
    - The kth-order neighbors are selected based on the condition d == order.

    Example
    -------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
    >>> get_kth_order_neighbours(G, 1, 2)
    [4]
    """
    neighbor_dict = nx.single_source_shortest_path_length(G, node, cutoff=order)
    kth_order_neighbors = [n for n, d in neighbor_dict.items() if d == order]
    return kth_order_neighbors
#
def get_neighbors_at_distance(G: nx.Graph, node: int, distance: int) -> list:
    """
    Returns the neighbors of a given node at a specific distance in a networkx graph.

    Parameters
    ----------
    G : Graph
        A NetworkX graph.
    node : int
        The node for which neighbors are to be found.
    distance : int
        The specific distance from the node.

    Returns
    -------
    neighbors : list
        A list of neighbors at the given distance from the given node.

    Notes
    -----
    - The graph G should be connected.
    - The function uses the single_source_shortest_path_length function of networkx to compute the shortest path lengths from the given node to all other nodes in the graph.
    - The neighbors are selected based on the condition d == distance.

    Example
    -------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
    >>> get_neighbors_at_distance(G, 1, 2)
    [4]
    """
    neighbor_dict = nx.single_source_shortest_path_length(G, node)
    neighbors_at_distance = [n for n, d in neighbor_dict.items() if d == distance]
    return neighbors_at_distance
#
def get_neighbors_within_distance(G: nx.Graph, node: Any, distance: int) -> list:
    """
    Returns the neighbors of a given node within a specific distance in a networkx graph.

    Parameters
    ----------
    G : Graph
        A NetworkX graph.
    node : Any
        The node for which neighbors are to be found.
    distance : int
        The maximum distance from the node.

    Returns
    -------
    neighbors : list
        A list of neighbors within the given distance from the given node.

    Notes
    -----
    - The graph G should be connected.
    - The function uses the single_source_shortest_path_length function of networkx to compute the shortest path lengths from the given node to all other nodes in the graph.
    - The neighbors are selected based on the condition d <= distance.

    Example
    -------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
    >>> get_neighbors_within_distance(G, 1, 2)
    [2, 3, 4]
    """
    neighbor_dict = nx.single_source_shortest_path_length(G, node)
    neighbors_within_distance = [n for n, d in neighbor_dict.items() if d <= distance]
    return neighbors_within_distance
#
def slaplacian_matrix_fromA(G: Graph, nodelist: list = None, weight: str = "weight"
                      ) -> csr_array:
    """Returns the signed Laplacian matrix of G.

    The graph Laplacian is the matrix L = |D| - A, where
    A is the adjacency matrix and |D| is the diagonal matrix of absolute values
    of node degrees.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : list, optional
       The rows and columns are ordered according to the nodes in nodelist.
       If nodelist is None, then the ordering is produced by G.nodes().

    weight : string or None, optional (default='weight')
       The edge data key used to compute each value in the matrix.
       If None, then each edge has weight 1.

    Returns
    -------
    L : SciPy sparse array
      The Laplacian matrix of G.
    """
    if nodelist is None:
        nodelist = list(G)
    A = nx.to_scipy_sparse_array(G, nodelist=nodelist, weight=weight,
                                 format="csr")
    D = csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 0, *A.shape,
                                    format="csr"))
    return D - A
#
def slaplacian_matrix(G: Graph, nodelist: list = None, weight: str = "weight"
                      ) -> csr_array:
    """Returns the signed Laplacian matrix of G.

    The graph Laplacian is the matrix L = |D| - A, where
    A is the adjacency matrix and |D| is the diagonal matrix of absolute values
    of node degrees.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : list, optional
       The rows and columns are ordered according to the nodes in nodelist.
       If nodelist is None, then the ordering is produced by G.nodes().

    weight : string or None, optional (default='weight')
       The edge data key used to compute each value in the matrix.
       If None, then each edge has weight 1.

    Returns
    -------
    L : SciPy sparse array
      The Laplacian matrix of G.
    """
    if nodelist is None:
        nodelist = list(G)
    A = nx.to_scipy_sparse_array(G, nodelist=nodelist, weight=weight,
                                 format="csr")
    D = csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 0, *A.shape,
                                    format="csr"))
    return D - A

def signed_spectral_layout(G, weight="weight", scale=1, center=None, dim=2):
    """Position nodes using the eigenvectors of the graph signed Laplacian.

    Using the unnormalized Laplacian, the layout shows possible clusters of
    nodes which are an approximation of the ratio cut. If dim is the number of
    dimensions then the positions are the entries of the dim eigenvectors
    corresponding to the ascending eigenvalues starting from the second one.

    Parameters
    ----------
    G : NetworkX graph or list of nodes
        A position will be assigned to every node in G.

    weight : string or None   optional (default='weight')
        The edge attribute that holds the numerical value used for
        the edge weight.  If None, then all edge weights are 1.

    scale : number (default: 1)
        Scale factor for positions.

    center : array-like or None
        Coordinate pair around which to center the layout.

    dim : int
        Dimension of layout.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> pos = nx.spectral_layout(G)

    Notes
    -----
    Directed graphs will be considered as undirected graphs when
    positioning the nodes.

    For larger graphs (>500 nodes) this will use the SciPy sparse
    eigenvalue solver (ARPACK).
    """
    # handle some special cases that break the eigensolvers
    import numpy as np

    G, center = _process_params(G, center, dim)

    if len(G) <= 2:
        if len(G) == 0:
            pos = np.array([])
        elif len(G) == 1:
            pos = np.array([center])
        else:
            pos = np.array([np.zeros(dim), np.array(center) * 2.0])
        return dict(zip(G, pos))
    try:
        # Sparse matrix
        if len(G) < 500:  # dense solver is faster for small graphs
            raise ValueError
        A = nx.to_scipy_sparse_array(G, weight=weight, dtype="d")
        # Symmetrize directed graphs
        if G.is_directed():
            A = A + np.transpose(A)
        pos = _sparse_spectral_signed(A, dim)
    except (ImportError, ValueError):
        # Dense matrix
        A = nx.to_numpy_array(G, weight=weight)
        # Symmetrize directed graphs
        if G.is_directed():
            A += A.T
        pos = _spectral_signed(A, dim)

    pos = rescale_layout(pos, scale=scale) + center
    pos = dict(zip(G, pos))
    return pos



def _spectral_signed(A, dim=2):
    # Input adjacency matrix A
    # Uses dense eigenvalue solver from numpy
    import numpy as np

    try:
        nnodes, _ = A.shape
    except AttributeError as err:
        msg = "spectral() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from err

    # form Laplacian matrix where D is diagonal of degrees
    D = np.identity(nnodes, dtype=A.dtype) * np.sum(np.abs(A), axis=1)
    L = D - A

    eigenvalues, eigenvectors = np.linalg.eig(L)
    # sort and keep smallest nonzero
    index = np.argsort(eigenvalues)[1 : dim + 1]  # 0 index is zero eigenvalue
    return np.real(eigenvectors[:, index])


def _sparse_spectral_signed(A, dim=2):
    # Input adjacency matrix A
    # Uses sparse eigenvalue solver from scipy
    # Could use multilevel methods here, see Koren "On spectral graph drawing"
    try:
        nnodes, _ = A.shape
    except AttributeError as err:
        msg = "sparse_spectral() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from err

    # form Laplacian matrix
    # TODO: Rm csr_array wrapper in favor of spdiags array constructor when available
    D = csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 0, nnodes, nnodes))
    L = D - A

    k = dim + 1
    # number of Lanczos vectors for ARPACK solver.What is the right scaling?
    ncv = max(2 * k + 1, int(np.sqrt(nnodes)))
    # return smallest k eigenvalues and eigenvectors
    eigenvalues, eigenvectors = scsp.linalg.eigsh(L, k, which="SM", ncv=ncv)
    index = np.argsort(eigenvalues)[1:k]  # 0 index is zero eigenvalue
    return np.real(eigenvectors[:, index])

def signedlaplacian_spectrum(G, weight="weight"):
    """Returns eigenvalues of the signed Laplacian of G

    Parameters
    ----------
    G : graph
       A NetworkX graph

    weight : string or None, optional (default='weight')
       The edge data key used to compute each value in the matrix.
       If None, then each edge has weight 1.

    Returns
    -------
    evals : NumPy array
      Eigenvalues

    Notes
    -----
    For MultiGraph/MultiDiGraph, the edges weights are summed.
    See to_numpy_array for other options.

    See Also
    --------
    laplacian_matrix
    """
    from scipy.linalg import eigvalsh
    return eigvalsh(slaplacian_matrix(G, weight=weight).todense())


def triangular_lattice_graph_modified(
    m, n, periodic=False, with_positions=True, create_using=None
):
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
    if n == 0 or m == 0:
        return G
    if periodic:
        if n < 5 or m < 3:
            msg = f"m > 2 and n > 4 required for periodic. m={m}, n={n}"
            raise NetworkXError(msg)

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
        # # add diagonals
        G.add_edges_from(((i, j), (i + 1, j + 1)) for j in rows[1:m-1:2] for i in cols[:n-1])
        G.add_edges_from(((i + 1, j), (i, j + 1)) for j in rows[:m-1:2] for i in cols[:n-1])

    # Add position node attributes
    if with_positions:
        ii = (i for i in cols for j in rows)
        jj = (j for i in cols for j in rows)
        xx = (0.5 * (j % 2) + i for i in cols for j in rows)
        h = np.sqrt(3) / 2
        if periodic:
            yy = (h * j + 0.01 * i * i for i in cols for j in rows)
        else:
            yy = (h * j for i in cols for j in rows)
        pos = {(i, j): (x, y) for i, j, x, y in zip(ii, jj, xx, yy) if (i, j) in G}
        set_node_attributes(G, pos, "pos")
    return G


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

def squared_lattice_graph_FastPatch(m, n, periodic=False, create_using=None, with_positions: bool = True):
    """Returns the two-dimensional grid graph.

    The grid graph has each node connected to its four nearest neighbors.

    Parameters
    ----------
    m, n : int 
        If an integer, nodes are from `range(n)`.

    periodic : bool or iterable
        If `periodic` is True, both dimensions are periodic. If False, none
        are periodic.  If `periodic` is iterable, it should yield 2 bool
        values indicating whether the 1st and 2nd axes, respectively, are
        periodic.

    create_using : NetworkX graph constructor, optional (default=nx.Graph)
        Graph type to create. If graph instance, then cleared before populated.

    Returns
    -------
    NetworkX graph
        The (possibly periodic) grid graph of the specified dimensions.

    """
    from networkx import set_node_attributes, empty_graph, NetworkXError
    from networkx.utils import pairwise
    G = empty_graph(0, create_using)
    rows = range(m)
    cols = range(n)
    G.add_nodes_from((i, j) for i in rows for j in cols)
    G.add_edges_from(((i, j), (pi, j)) for pi, i in pairwise(rows) for j in cols)
    G.add_edges_from(((i, j), (i, pj)) for i in rows for pj, j in pairwise(cols))

    try:
        periodic_r, periodic_c = periodic
    except TypeError:
        periodic_r = periodic_c = periodic

    if periodic_r and len(rows) > 2:
        first = rows[0]
        last = rows[-1]
        G.add_edges_from(((first, j), (last, j)) for j in cols)
    if periodic_c and len(cols) > 2:
        first = cols[0]
        last = cols[-1]
        G.add_edges_from(((i, first), (i, last)) for i in rows)
    # both directions for directed
    if G.is_directed():
        G.add_edges_from((v, u) for u, v in G.edges())
    if with_positions:
        # calc position in embedded space
        ii = (i for i in rows for j in cols)
        jj = (j for i in rows for j in cols)
        xx = (i for i in rows for j in cols)  # x position matches column index
        yy = (j for i in rows for j in cols)
        if periodic:
            xx = (i + 0.02 * j * j for i in rows for j in cols)  # x position matches column index
            yy = (j + 0.02 * i * i for i in rows for j in cols)
        # exclude nodes not in G
        pos = {(i, j): (x, y) for i, j, x, y in zip(ii, jj, xx, yy) if (i, j) in G}
        set_node_attributes(G, pos, "pos")
    return G
