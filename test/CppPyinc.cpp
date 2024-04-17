#include <boost/python.hpp>
#include "CppPyinc.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(CppPyinc)
{
    class_<World>("World")
        .def("greet", &World::greet)
        .def("set", &World::set)
    ;
}