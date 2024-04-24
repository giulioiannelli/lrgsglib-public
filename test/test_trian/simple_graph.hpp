#ifndef SIMPLE_GRAPH_HPP
#define SIMPLE_GRAPH_HPP

#include <graph.hh>

using namespace graph_tool;

class SimpleGraph {
public:
    SimpleGraph(); // Constructor to initialize an empty graph
    GraphInterface& get_graph(); // Getter for the internal graph

private:
    GraphInterface graph;
};

#endif // SIMPLE_GRAPH_HPP
