
CXXFLAGS = -stdlib=libc++ -std=c++11 -W -Wall -O3 -march=native
CXX = clang++

testbench: testbench.cc reed_solomon.hh bose_chaudhuri_hocquenghem.hh berlekamp_massey.hh chien.hh forney.hh correction.hh galois_field.hh galoisfieldtables.hh
	$(CXX) $(CXXFLAGS) -g $< -o $@

benchmark: testbench.cc reed_solomon.hh bose_chaudhuri_hocquenghem.hh berlekamp_massey.hh chien.hh forney.hh correction.hh galois_field.hh galoisfieldtables.hh
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

