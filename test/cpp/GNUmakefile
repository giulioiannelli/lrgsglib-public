export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

PYTHON_INCLUDE = $(shell python3 -c "from sysconfig import get_paths as gp; print(gp()['include'])")
PYTHON_LIB = $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")

CXX := g++
LDFLAGS := -Wl,--no-undefined
BOOST_LIB := -lboost_python39
CXXFLAGS := $(shell pkg-config --cflags --libs graph-tool-py3.9) $(BOOST_LIB) -fPIC -shared -I$(PYTHON_INCLUDE)
CXXFLAGS += -Wall -Wextra -v $(LDFLAGS)

all: triangular_lattice.so

triangular_lattice.so: triangular_lattice.cc
	$(CXX) $(CXXFLAGS) $< -o $@ 

clean:
	rm -f *.so
