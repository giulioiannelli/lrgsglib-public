import simple_graph_ext

def test_graph():
    sg = simple_graph_ext.SimpleGraph()
    graph = sg.get_graph()
    print("Graph created with", graph.num_vertices(), "vertices and", graph.num_edges(), "edges.")

if __name__ == "__main__":
    test_graph()