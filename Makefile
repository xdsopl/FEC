
CXXFLAGS = -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native
CXX = clang++

benchmark: benchmark.cc reedsolomon.hh galoisfield.hh
	$(CXX) $(CXXFLAGS) $< -o $@

test: benchmark
	./benchmark

.PHONY: clean test

clean:
	rm -f benchmark

