CXXFLAGS := -std=c++0x -Wall -Wextra -O3 -march=native -I. -Igoogletest/googletest/include -Igoogletest/googletest -Igoogletest/googletest/src -pthread
LDFLAGS := -pthread
ifeq ($(buildtype),coverage)
	CXXFLAGS += -coverage
	LDFLAGS += -coverage
endif

INCS := alglib.hpp
OBJS := sample.o googletest/googletest/src/gtest-all.o
PROG := alglib

$(PROG): $(INCS) $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(PROG)

all: clean $(PROG)

test: $(PROG)
	./$(PROG)

googletest/googletest/src/gtest-all.o: googletest/googletest/src/gtest-all.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	$(RM) $(PROG) $(OBJS)
	find . -name "*.gc??" -exec rm -rf {} \;

get-deps:
	pip install --user cpp-coveralls

coveralls:
	coveralls --exclude gtest --gcov-options '\-lp'
