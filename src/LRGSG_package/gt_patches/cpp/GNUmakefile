export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

CXX := g++
CXXFLAGS := $(shell pkg-config --cflags --libs graph-tool-py3.9) -fPIC -shared

all: triangular_lattice.so

triangular_lattice.so: triangular_lattice.cc
	$(CXX) $(CXXFLAGS) $< -o $@ 

clean:
	rm -f *.so
