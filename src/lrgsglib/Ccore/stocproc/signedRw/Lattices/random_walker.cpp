#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <cmath>

namespace py = boost::python;

// Converter from Python tuple to std::pair<int, int>
struct pair_to_python_tuple {
    static PyObject* convert(const std::pair<int, int>& p) {
        return py::incref(py::make_tuple(p.first, p.second).ptr());
    }
};

struct pair_from_python_tuple {
    pair_from_python_tuple() {
        py::converter::registry::push_back(
            &convertible,
            &construct,
            py::type_id<std::pair<int, int>>());
    }

    static void* convertible(PyObject* obj_ptr) {
        if (!PyTuple_Check(obj_ptr)) return 0;
        if (PyTuple_Size(obj_ptr) != 2) return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((py::converter::rvalue_from_python_storage<std::pair<int, int>>*)data)->storage.bytes;
        int first = py::extract<int>(PyTuple_GetItem(obj_ptr, 0));
        int second = py::extract<int>(PyTuple_GetItem(obj_ptr, 1));
        new (storage) std::pair<int, int>(first, second);
        data->convertible = storage;
    }
};

std::pair<int, int> normalize_edge(int a, int b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
}

class RandomWalker {
public:
    RandomWalker(py::list edge_list, py::dict sign_dict, int width, int height, int coordination_number)
        : width(width), height(height), coordination_number(coordination_number) {
        // srand(time(0));
        initializeNeighbors(edge_list);
        initializeSigns(sign_dict);
    }

    void walk(int start) {
        int current = start;
        int start_x = current / width;
        int start_y = current % width;
        int current_x, current_y;
        distance = 0;
        euclidean_distance = 0.0;
        stopped_node = current; // Initialize the stopping node

        while (true) {
            const auto& neighbors = this->neighbors[current];
            if (isBorder(current)) break;

            int next = neighbors[rand() % neighbors.size()];
            auto normalized_edge = normalize_edge(current, next);
            if (signs[normalized_edge] == -1) break;

            current = next;
            stopped_node = current; // Update the stopping node
            distance++;

        }
        current_x = current / width;
        current_y = current % width;
        euclidean_distance = std::sqrt(std::pow(current_x - start_x, 2) + std::pow(current_y - start_y, 2));
    }

    int get_step_distance() const {
        return distance;
    }

    double get_euclidean_distance() const {
        return euclidean_distance;
    }

    int get_stopped_node() const {
        return stopped_node;
    }

private:
    std::map<std::pair<int, int>, int> signs;
    std::map<int, std::vector<int>> neighbors;
    int distance = 0;
    double euclidean_distance = 0.0;
    int width;
    int height;
    int coordination_number;
    int stopped_node;

    void initializeNeighbors(py::list edge_list) {
        for (int i = 0; i < len(edge_list); ++i) {
            py::object item = edge_list[i];
            std::pair<int, int> edge = py::extract<std::pair<int, int>>(item);
            auto normalized_edge = normalize_edge(edge.first, edge.second);
            neighbors[normalized_edge.first].push_back(normalized_edge.second);
            neighbors[normalized_edge.second].push_back(normalized_edge.first);
        }
    }

    void initializeSigns(py::dict sign_dict) {
        py::list keys = sign_dict.keys();
        for (int i = 0; i < len(keys); ++i) {
            py::object key_obj = keys[i];
            std::pair<int, int> key = py::extract<std::pair<int, int>>(key_obj);
            auto normalized_key = normalize_edge(key.first, key.second);
            int sign = py::extract<int>(sign_dict[key_obj]);
            this->signs[normalized_key] = sign;
        }
    }

    bool isBorder(int node) {
        return neighbors[node].size() < coordination_number;
    }
};

BOOST_PYTHON_MODULE(random_walker) {
    py::to_python_converter<std::pair<int, int>, pair_to_python_tuple>();
    pair_from_python_tuple();

    py::class_<RandomWalker>("RandomWalker", py::init<py::list, py::dict, int, int, int>())
        .def("walk", &RandomWalker::walk)
        .def("get_step_distance", &RandomWalker::get_step_distance)
        .def("get_euclidean_distance", &RandomWalker::get_euclidean_distance)
        .def("get_stopped_node", &RandomWalker::get_stopped_node);
}
