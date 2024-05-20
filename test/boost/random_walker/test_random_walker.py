import random_walker
import networkx as nx
import numpy as np

def generate_edge_list(lattice_size):
    G = nx.grid_2d_graph(lattice_size, lattice_size)
    H = nx.convert_node_labels_to_integers(G)
    edges = list(H.edges())
    return edges

def main():
    lattice_size = 10
    p = 0.5
    num_disorder_realizations = 1
    num_walks_per_realization = 1000
    
    edge_list = generate_edge_list(lattice_size)
    # print(edge_list)
    walker = random_walker.RandomWalker(edge_list, p)
    
    final_positions = []
    
    for _ in range(num_disorder_realizations):
        walker.initialize()
        for _ in range(num_walks_per_realization):
            start_node = 50
            final_positions.append(walker.walk(start_node))
    
    print("Average final position:", np.mean(final_positions))

if __name__ == "__main__":
    main()
