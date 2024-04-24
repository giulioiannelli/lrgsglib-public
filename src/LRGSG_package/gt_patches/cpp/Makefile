export PKG_CONFIG_PATH:=$(CONDA_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)

BOOST_LIB := -lboost_python312

CXX := gcc
CXXFLAGS := $(shell pkg-config --cflags --libs graph-tool-py3.12) $(BOOST_LIB) -fPIC -shared

all: triangular_lattice.so

triangular_lattice.so: triangular_lattice.cc
	$(CXX) $(CXXFLAGS) $< -o $@ 

clean:
	rm -f *.so
