CXXFLAGS := -std=c++0x -Wall -Wextra -O3 -march=native -I. -pthread
LDFLAGS := -pthread
ifeq ($(buildtype),coverage)
	CXXFLAGS += -coverage
	LDFLAGS += -coverage
endif

INCS := alglib.hpp
OBJS := sample.o gtest/gtest-all.o
PROG := alglib

$(PROG): $(INCS) $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(PROG)

all: clean $(PROG)

test: $(PROG)
	./$(PROG)

clean:
	$(RM) $(PROG) $(OBJS)
	find . -name "*.gc??" -exec rm -rf {} \;

get-deps:
	pip install --user cpp-coveralls

coveralls:
	coveralls --exclude gtest --gcov-options '\-lp'
