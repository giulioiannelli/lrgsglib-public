#include <boost/python.hpp>
#include "simple_graph.hpp"

SimpleGraph::SimpleGraph() : graph(false, true) {}

GraphInterface& SimpleGraph::get_graph() {
    return graph;
}

BOOST_PYTHON_MODULE(simple_graph_ext) {
    using namespace boost::python;
    class_<SimpleGraph>("SimpleGraph", init<>())
        .def("get_graph", &SimpleGraph::get_graph, return_internal_reference<>());
}
