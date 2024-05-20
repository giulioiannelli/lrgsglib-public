#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/list.hpp>
#include <vector>
#include <random>
#include <unordered_set>

using namespace boost::python;

struct Edge {
    int u, v;
    bool sign; // true for positive, false for negative
};

class RandomWalker {
public:
    RandomWalker(const std::vector<std::pair<int, int>>& edges, double p)
        : edges_(edges), p_(p), rng_(std::random_device{}()) {
        initialize();
    }

    // Conversion constructor: takes a Python list
    RandomWalker(boost::python::list& edge_list, double p)
        : p_(p), rng_(std::random_device{}()) {
        for (int i = 0; i < len(edge_list); ++i) {
            boost::python::tuple edge = extract<boost::python::tuple>(edge_list[i]);
            int u = extract<int>(edge[0]);
            int v = extract<int>(edge[1]);
            edges_.push_back(std::make_pair(u, v));
        }
        initialize();
    }

    void initialize() {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        signed_edges_.clear();
        for (const auto& edge : edges_) {
            signed_edges_.push_back({edge.first, edge.second, dist(rng_) >= p_});
        }
    }

    std::pair<int, bool> step(int current_node) {
        std::unordered_set<int> neighbors;
        for (const auto& edge : signed_edges_) {
            if (edge.u == current_node && edge.sign) neighbors.insert(edge.v);
            if (edge.v == current_node && edge.sign) neighbors.insert(edge.u);
        }
        if (neighbors.empty()) return {current_node, false};

        std::uniform_int_distribution<int> dist(0, neighbors.size() - 1);
        auto it = neighbors.begin();
        std::advance(it, dist(rng_));
        return {*it, true};
    }

    int walk(int start_node) {
        int current_node = start_node;
        while (true) {
            auto result = step(current_node);
            if (!result.second || is_border_node(current_node)) break;
            current_node = result.first;
        }
        return current_node;
    }

private:
    bool is_border_node(int node) {
        int degree = 0;
        for (const auto& edge : signed_edges_) {
            if (edge.u == node || edge.v == node) degree++;
        }
        return degree != coordination_number_;
    }

    std::vector<std::pair<int, int>> edges_;
    std::vector<Edge> signed_edges_;
    double p_;
    int coordination_number_ = 4; // Example for a 2D lattice
    std::mt19937 rng_;
};

BOOST_PYTHON_MODULE(random_walker) {
    class_<RandomWalker>("RandomWalker", init<boost::python::list&, double>())
        .def("initialize", &RandomWalker::initialize)
        .def("walk", &RandomWalker::walk);
}
