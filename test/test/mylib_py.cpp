#include <boost/python.hpp>
#include "mylib.hpp"

BOOST_PYTHON_MODULE(mylib_ext) {
    using namespace boost::python;
    class_<MyLib>("MyLib")
        .def("add", &MyLib::add)
    ;
}
