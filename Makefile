CXX = g++
CXXFLAGS = -O3 -march=native -flto -fopenmp -std=c++17 -Wall -Wextra
CPPFLAGS += -Iinclude

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

SOURCES = main.cpp $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = picaso

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(ROOTLIBS)
	$(RM) $(OBJECTS)

main.o: main.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

src/%.o: src/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

debug:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="-O0 -g -fopenmp -std=c++17 -Wall -Wextra -DDEBUG_MODE" all

clean:
	$(RM) $(OBJECTS) $(TARGET)

.PHONY: all debug clean
