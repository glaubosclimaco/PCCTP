# Makefile para PCCTP - projeto organizado em src/ e include/
# Dependencias: Gurobi, LEMON (liblemon-dev)

CXX = g++
TARGET = pcctp

# Fontes em src/ (main em hib.cpp)
SRCS = src/hib.cpp src/grasp.cpp src/genius.cpp src/instances.cpp src/instancia_loader.cpp \
       src/mineracao.cpp src/PR.cpp src/random_provider.cpp src/Utils.cpp src/Graph.cpp \
       src/MinCut.cpp src/Regras.cpp

OBJS = $(SRCS:.cpp=.o)

# Includes: cabeçalhos em include/
INCLUDES = -I include

# Gurobi
GUROBI_HOME ?= $(firstword $(wildcard /opt/gurobi*/linux64) /opt/gurobi752/linux64)
GUROBI_INC = -I$(GUROBI_HOME)/include
GUROBI_LIB_SO = $(firstword $(wildcard $(GUROBI_HOME)/lib/libgurobi*.so))
GUROBI_LIB = -L$(GUROBI_HOME)/lib -lgurobi_c++ $(if $(GUROBI_LIB_SO),-l$(patsubst lib%.so,%,$(notdir $(GUROBI_LIB_SO))),-lgurobi) -lm

# LEMON
LEMON_CFLAGS = $(shell pkg-config --cflags lemon 2>/dev/null || echo "-I/usr/include")
LEMON_LIBS   = $(shell pkg-config --libs lemon 2>/dev/null || echo "-lemon")

CXXFLAGS = -std=c++11 -Wall $(INCLUDES) $(GUROBI_INC) $(LEMON_CFLAGS)
LDFLAGS = $(GUROBI_LIB) $(LEMON_LIBS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

# Regra para src/*.cpp -> src/*.o
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<


clean:
	rm -f $(OBJS) $(TARGET)
	find src -name '*.o' -delete
	rm -f *.o
