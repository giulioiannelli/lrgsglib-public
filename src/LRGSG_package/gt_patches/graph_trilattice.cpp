#include "graph-tool/src/graph/graph.hh"

void create_triangular_lattice(graph_tool::GraphInterface& gi, size_t width, size_t height)
{
    auto g = gi.get_graph();
    // Assuming a straightforward triangular lattice
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            int current = i * width + j;
            if (j < width - 1) // Connect right
                add_edge(g, current, current + 1);
            if (i < height - 1) // Connect down
                add_edge(g, current, current + width);
            if (i < height - 1 && j < width - 1) // Connect diagonal
                add_edge(g, current, current + width + 1);
        }
    }
}