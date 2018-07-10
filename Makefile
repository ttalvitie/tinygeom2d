.PHONY: clean

test: test.cpp $(wildcard *.hpp)
	$(CXX) test.cpp -o test -std=c++11 -Wall -Wextra -g

clean:
	rm -f test
