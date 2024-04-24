CXX = g++
# Extract flags for graph-tool using pkg-config, ensuring we're targeting the correct Python version
GRAPH_TOOL_FLAGS = $(shell pkg-config --cflags --libs graph-tool-py3.12)
GRAPH_TOOL_INCLUDE = -I/home/opisthofulax/anaconda3/envs/lrgsgenv/lib/python3.12/site-packages/graph_tool/include
CXXFLAGS = -O2 -Wall -fPIC $(GRAPH_TOOL_FLAGS) $(GRAPH_TOOL_INCLUDE) -shared

TARGET = simple_graph_ext.so
SOURCES = simple_graph.cpp simple_graph_py.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@

clean:
	rm -f $(TARGET)
