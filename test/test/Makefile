CXX = g++
CXXFLAGS = -O2 -fPIC -shared -std=c++11
PYTHON_INC = $(shell python -m pybind11 --includes)
BOOST_INC = -I $(CONDA_PREFIX)/include

BOOST_LIB = -L $(CONDA_PREFIX)/lib -lboost_python$(shell python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')") -lboost_system

all: mylib_ext.so

mylib_ext.so: mylib.cpp mylib_py.cpp
	$(CXX) $(CXXFLAGS) $(PYTHON_INC) $(BOOST_INC) mylib.cpp mylib_py.cpp -o mylib_ext.so $(BOOST_LIB)

clean:
	rm -f mylib_ext.so
