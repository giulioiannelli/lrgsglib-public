#include "simple_graph.hpp"

SimpleGraph::SimpleGraph() : graph(false, true) {}

GraphInterface& SimpleGraph::get_graph() {
    return graph;
}
