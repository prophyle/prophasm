CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lz

.PHONY: all clean

all: ../prophasm

../prophasm: prophasm.o
	$(CXX) $(CXXFLAGS) $(DFLAGS) $^ -o $@ -L. $(LIBS)

prophasm.o: prophasm.cpp $(wildcard *.h)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $<

clean:
	rm -f prophasm *.o
