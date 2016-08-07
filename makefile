CXXFLAGS := -std=c++0x -Wall -Wextra -O3 -march=native -I. -pthread

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
