CXX = g++
PYTHON_VERSION = 3.10
# PYTHON_INCLUDE = $(shell python$(PYTHON_VERSION)-config --includes)

PYTHON_INCLUDE = $(shell python3 -c "from sysconfig import get_paths as gp; print(gp()['include'])")
PYTHON_LIB = $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")

CXXFLAGS = -fPIC -shared -I$(PYTHON_INCLUDE)
LDFLAGS = -lboost_python310 -lpython3.10

CppPyinc.so: CppPyinc.cpp CppPyinc.hpp
	$(CXX) $(CXXFLAGS) CppPyinc.cpp -o CppPyinc.so $(LDFLAGS)

clean:
	rm -f *.so