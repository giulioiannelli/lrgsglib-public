import networkx as nx
import numpy as np
import scipy as sp
#
import scipy.sparse as scsp
#
from networkx.classes.graph import Graph
from networkx.drawing.layout import _process_params, rescale_layout
from scipy.sparse import csr_array
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
def get_neighbors_within_distance(G: nx.Graph, node: int, distance: int) -> list:
    """
    Returns the neighbors of a given node within a specific distance in a networkx graph.

    Parameters
    ----------
    G : Graph
        A NetworkX graph.
    node : int
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
    D = scsp.csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 0, *A.shape,
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
    D = scsp.csr_array(scsp.spdiags(np.abs(A).sum(axis=1), 0, *A.shape,
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
    import numpy as np
    import scipy as sp
    import scipy.sparse  # call as sp.sparse
    import scipy.sparse.linalg  # call as sp.sparse.linalg

    try:
        nnodes, _ = A.shape
    except AttributeError as err:
        msg = "sparse_spectral() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from err

    # form Laplacian matrix
    # TODO: Rm csr_array wrapper in favor of spdiags array constructor when available
    D = sp.sparse.csr_array(sp.sparse.spdiags(np.abs(A).sum(axis=1), 0, nnodes, nnodes))
    L = D - A

    k = dim + 1
    # number of Lanczos vectors for ARPACK solver.What is the right scaling?
    ncv = max(2 * k + 1, int(np.sqrt(nnodes)))
    # return smallest k eigenvalues and eigenvectors
    eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(L, k, which="SM", ncv=ncv)
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
    return sp.linalg.eigvalsh(slaplacian_matrix(G, weight=weight).todense())