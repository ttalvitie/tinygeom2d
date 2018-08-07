CXX ?= g++
SRCS := $(wildcard *.cpp)
BINS := $(SRCS:%.cpp=%)

.PHONY: all clean

all: $(BINS)

define BIN_RULE
$(1): $(1).cpp svg.hpp $(wildcard tinygeom2d/*.hpp)
	$(CXX) $(1).cpp -o $(1) -std=c++11 -Wall -Wextra -O2 -I../
endef
$(foreach b,$(BINS),$(eval $(call BIN_RULE,$(b))))

clean:
	rm -f $(BINS)
