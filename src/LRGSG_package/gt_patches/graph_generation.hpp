// Include guard to prevent multiple inclusion
#ifndef GRAPH_TOOL_GENERATION_TRIANGULAR_LATTICE_HPP
#define GRAPH_TOOL_GENERATION_TRIANGULAR_LATTICE_HPP

#include "graph-tool/src/graph/graph.hh"


namespace graph_tool {

    void create_triangular_lattice(GraphInterface& gi, size_t width, size_t height);

} // namespace graph_tool

#endif // GRAPH_TOOL_GENERATION_TRIANGULAR_LATTICE_HPP