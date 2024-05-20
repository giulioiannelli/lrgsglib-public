#include <boost/python.hpp>

std::string greet() {
    return "Hello, World!";
}

BOOST_PYTHON_MODULE(hello) {
    using namespace boost::python;
    def("greet", greet);
}
