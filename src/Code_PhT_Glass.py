import os
import sys
import random
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
from typing import Any


def flip_sel_edges(G, links):
    """Flips specific edges of a graph G."""
    #
    neg_weights_dict = {}
    if links:
        neg_weights_dict = {
            (u, v): -1
            * G.get_edge_data(u, v)["weight"]
            for u, v in links
        }
    nx.set_edge_attributes(
        G, values=neg_weights_dict, name="weight"
    )


def get_links_XERR(node, G):
    return [(node, nn) for nn in G.neighbors(node)]

    
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

def main():
    geo = 'hex'
    L=int(sys.argv[1])
    p1=float(sys.argv[2])
    navrg=100
    smax=np.zeros((navrg))
    scnd=np.zeros((navrg))
    trhd=np.zeros((navrg))
    four=np.zeros((navrg))
    var=0.0
    if geo == 'hex':
        side1 = adjust_to_even(L/np.sqrt(3))
        H = hexagonal_lattice_graph_FastPatch(side1, L, periodic=True, with_positions=True)
    else:
        side1 = L
        H = nx.grid_2d_graph(side1, side1, periodic=True)
    mapping = {old_label:new_label for new_label, old_label in enumerate(H.nodes())}
    G = nx.relabel_nodes(H, mapping)
    # Calculate adjacency matrix
    adj_matrix = nx.adjacency_matrix(G)
    # Compute the vector of degrees
    diabs = np.array(adj_matrix.sum(axis=1)).flatten()
    N=len(G.nodes())
    if geo == 'hex':
        file1='data/test_hexerr/PhaseT_Hex_L_'+str(L)+'_p1_'+str(p1)+'_'
    else:
        file1='data/test_sqerr/PhaseT_Sq_L_'+str(L)+'_p1_'+str(p1)+'_'
    cont=0
    for avrg in range(navrg):
        # Iterate over all edges and set weights to 1
        nx.set_edge_attributes(G, values=1, name="weight")
        node_list_flip = random.sample(list(G.nodes), int(G.number_of_nodes()*p1)+1)
        patternList = [k for i in node_list_flip 
                    for k in get_links_XERR(i, G)]
        flip_sel_edges(G, patternList)
        # Calculate adjacency matrix
        adj = nx.adjacency_matrix(G, weight='weight',dtype=np.float64)
        # Laplacian matrix
        slapl = diags(diabs, 0, format='csr', dtype=np.float64) - adj
        # Find the first eigenvector
        _, eigV = eigsh(slapl, k=1, which='SM')
        # Count the number of negative and positive elements
        negative_count = np.sum(eigV < 0)
        positive_count = np.sum(eigV > 0)
        # Determine majority sign
        minority_sign = -1 if negative_count > positive_count else 1
        # Find positions of elements with majority sign
        if minority_sign == -1:
            min_pos = np.where(eigV < 0)[0]  # Get positions of negative values
        else:
            min_pos = np.where(eigV > 0)[0]  # Get positions of positive values

        # Extract subgraph with nodes at min_pos
        subG = G.subgraph(min_pos)
        # Find the connected components of the subgraph
        clust = np.array([len(component) for component in nx.connected_components(subG)])
        srtclust=sorted(clust)
        if len(clust) > 1:
            smax[cont] = np.max(clust) / (1.0 * N)
            scnd[cont] = srtclust[-2] / (1.0 * N)
            if len(clust) > 2:
                trhd[cont] = srtclust[-3] / (1.0 * N)
            if len(clust) > 3:
                four[cont] = srtclust[-4] / (1.0 * N)
            var += ((np.sum(np.multiply(clust,clust))- np.max(clust)*np.max(clust))/(np.sum(clust) - np.max(clust)))
        else:
            smax[cont] = 1.0*len(clust)
            var += 0.0
        if (avrg % 100) == 0:
            x = [np.mean(smax[0:cont]), np.std(smax[0:cont]), np.mean(scnd[0:cont]), np.std(scnd[0:cont]), np.mean(trhd[0:cont]), np.std(trhd[0:cont]),
            np.mean(four[0:cont]), np.std(four[0:cont]), var / (avrg + 1), p1, int(avrg + 1)]
            np.savetxt(file1, np.atleast_2d(x), delimiter=',')
        cont=cont+1
    x = [np.mean(smax[0:cont]), np.std(smax[0:cont]), np.mean(scnd[0:cont]), np.std(scnd[0:cont]), np.mean(trhd[0:cont]), np.std(trhd[0:cont]),
            np.mean(four[0:cont]), np.std(four[0:cont]), var / (avrg + 1), p1, int(avrg + 1)]
    np.savetxt(file1, np.atleast_2d(x), delimiter=',')
    # Close the file
if __name__ == "__main__":
    main()