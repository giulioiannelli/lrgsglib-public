CXX := g++
CXXFLAGS := -std=c++17 -fPIC -shared -Wl --no-undefined
INCLUDES := -I$(CONDA_PREFIX)/include -I$(CONDA_PREFIX)/include/python3.12
LIBS := -L$(CONDA_PREFIX)/lib -lboost_python312

CppPyinc.so: CppPyinc.cpp CppPyinc.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) CppPyinc.cpp -o CppPyinc.so

clean:
	rm -f *.so

clean:
	rm -f *.so