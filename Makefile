
CXXFLAGS = -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native
CXX = clang++

testbench: testbench.cc reedsolomon.hh galoisfield.hh galoisfieldtables.hh
	$(CXX) $(CXXFLAGS) -g $< -o $@

benchmark: testbench.cc reedsolomon.hh galoisfield.hh galoisfieldtables.hh
	$(CXX) $(CXXFLAGS) -DNDEBUG $< -o $@

tablesgenerator: tablesgenerator.cc
	$(CXX) $(CXXFLAGS) $< -o $@

galoisfieldtables.hh: tablesgenerator
	./tablesgenerator > $@

test: testbench
	./testbench

speed: benchmark
	./benchmark

.PHONY: clean test

clean:
	rm -f benchmark testbench tablesgenerator galoisfieldtables.hh

