
CXXFLAGS = -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native
CXX = clang++

testbench: testbench.cc reedsolomon.hh bose_chaudhuri_hocquenghem.hh berlekampmassey.hh chien.hh forney.hh correction.hh galoisfield.hh galoisfieldtables.hh
	$(CXX) $(CXXFLAGS) -g $< -o $@

benchmark: testbench.cc reedsolomon.hh bose_chaudhuri_hocquenghem.hh berlekampmassey.hh chien.hh forney.hh correction.hh galoisfield.hh galoisfieldtables.hh
	$(CXX) $(CXXFLAGS) -DNDEBUG $< -o $@

tablesgenerator: tablesgenerator.cc
	$(CXX) $(CXXFLAGS) $< -o $@

galoisfieldtables.hh: tablesgenerator
	./tablesgenerator > $@

test: testbench
	uname -p
	./testbench

speed: benchmark
	uname -p | tee RESULTS
	./benchmark | tee -a RESULTS

.PHONY: clean test

clean:
	rm -f benchmark testbench tablesgenerator galoisfieldtables.hh

