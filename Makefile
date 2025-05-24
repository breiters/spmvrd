SHELL:=/bin/bash

BIN:=spmvrd

HEADERS:=$(wildcard *.h)
SOURCES:=$(wildcard *.cpp)
OBJECTS:=$(SOURCES:.cpp=.o)

CXXFLAGS?=-std=c++20 -fopenmp -Ofast -march=native -mtune=native -flto

CXXFLAGS+=-DNDEBUG
# CXXFLAGS+=-pg -fno-inline # for profiling / debugging
CXXFLAGS+=-g3 -Wall -Wextra -Wpedantic
# CXXFLAGS+=-Weffc++
CXXFLAGS+=-Wno-sign-compare -Wno-sign-conversion
CXXFLAGS+=-Wno-unused-parameter -Wno-unused-function
CXXFLAGS+=-Wno-unused-variable -Wno-unused-but-set-variable
CXXFLAGS+=-Wno-maybe-uninitialized
# CXXFLAGS+=-Wconversion -Wdouble-promotion

# CXXFLAGS+=-fsanitize=undefined #-fsanitize-trap
# CXXFLAGS+=-fsanitize=address
# CXXFLAGS+=-fsanitize=thread

CXXFLAGS+=-DMEMBLOCKLEN=$(shell getconf LEVEL1_DCACHE_LINESIZE)
CXXFLAGS+=-DCACHE_LINESIZE=$(shell getconf LEVEL1_DCACHE_LINESIZE)

CXX?=g++

$(BIN): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(BIN)
