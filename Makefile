
CXXFLAGS = -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native
CXX = clang++

benchmark: benchmark.cc reedsolomon.hh galoisfield.hh galoisfieldtables.hh
	$(CXX) $(CXXFLAGS) $< -o $@

generator: generator.cc
	$(CXX) $(CXXFLAGS) $< -o $@

galoisfieldtables.hh: generator
	./generator > $@

test: benchmark
	./benchmark

.PHONY: clean test

clean:
	rm -f benchmark generator galoisfieldtables.hh

