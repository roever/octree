all: example

example: example.cpp octree.hpp
	clang++ -g --std=c++1z example.cpp -o example -O3 -W -Wall


