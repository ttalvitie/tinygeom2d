CXX ?= g++

.PHONY: clean

test: test.cpp $(wildcard ../tinygeom2d/*.hpp)
	$(CXX) test.cpp -o test -std=c++11 -Wall -Wextra -g -I../

clean:
	rm -f test
