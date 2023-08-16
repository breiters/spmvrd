SHELL:=/bin/bash

BIN:=spmvrd

HEADERS:=$(wildcard *.h)
SOURCES:=$(wildcard *.cpp)
OBJECTS:=$(SOURCES:.cpp=.o)

CXXFLAGS?=-std=c++17 -fconcepts-ts -fopenmp -Ofast -march=native -mtune=native -flto #-fwhole-program

CXXFLAGS+=-DNDEBUG
CXXFLAGS+=-g3 -Wall -Wextra
CXXFLAGS+=-Wno-sign-compare -Wno-sign-conversion
CXXFLAGS+=-Wno-unused-parameter -Wno-unused-function
CXXFLAGS+=-Wno-unused-variable -Wno-unused-but-set-variable
CXXFLAGS+=-Wno-maybe-uninitialized
# CXXFLAGS+=-Wconversion -Wdouble-promotion
# CXXFLAGS+=-fsanitize=undefined #-fsanitize-trap
# CXXFLAGS+=-fsanitize=address
# CXXFLAGS+=-fsanitize=thread

CXX?=g++

$(BIN): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(BIN)
