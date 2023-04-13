SHELL:=/bin/bash

BIN:=spmvrd

HEADERS:=$(wildcard *.h)
SOURCES:=$(wildcard *.cpp)
OBJECTS:=$(SOURCES:.cpp=.o)

CXXFLAGS?=-std=c++17 -g -Wall -Wextra -fopenmp -Ofast -march=native -mtune=native
CXXFLAGS+=-DNDEBUG
# CXXFLAGS?=-std=c++17 -g -Wall -Wextra -fopenmp -O0
CXX?=g++

$(BIN): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(BIN)
