c++ = g++
dependency = matrix.hpp dmatrix.hpp complex.hpp cmatrix.hpp
example = ./examples/example.cpp
flags = -O2 -Wall -std=c++11
comp_only = -c
debug = -DDEBUG

example: $(example) $(dependency)
	$(c++) $< $(flags) -o $@

.PHONY: all
