/*
reedsolomon - Reed-Solomon error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <iostream>
#include <cassert>
#include "galoisfield.hh"
#include "reedsolomon.hh"

template <typename TYPE>
void print_table(TYPE *table, const char *name, int N)
{
	std::cout << name << "[" << N << "] = {";
	for (int i = 0; i < N; ++i) {
		std::cout << " " << (int)table[i];
		if (i < N - 1)
			std::cout << ",";
	}
	std::cout << " };" << std::endl;
}

template <int NR, int FR, int M, int P, typename TYPE>
void test(ReedSolomon<NR, FR, GF::Types<M, P, TYPE>> &rs, TYPE *code, TYPE *target)
{
	//print_table(rs.gf.exp_table, "exp_table", rs.gf.Q);
	//print_table(rs.gf.log_table, "log_table", rs.gf.Q);
#if 0
	std::cout << "g(x) = ";
	for (int i = NR; i > 0; --i) {
		if (rs.G[i] != 1)
			std::cout << (int)rs.G[i] << "*";
		std::cout << "x";
		if (i != 1)
			std::cout << "^" << i;
		std::cout << " + ";
	}
	std::cout << (int)rs.G[0] << std::endl;
#endif
	rs.encode(code);
	for (int i = 0; i < rs.N; ++i)
		assert(code[i] == target[i]);
	//print_table(code, "parity", NR);
}

int main()
{
	{ // BBC WHP031 RS(15, 11) T=2
		ReedSolomon<4, 0, GF::Types<4, 0b10011, uint8_t>> rs;
		uint8_t code[15] = { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
		uint8_t target[15] = { 3, 3, 12, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
		test(rs, code, target);
	}
	{ // DVB-T RS(255, 239) T=8
		ReedSolomon<16, 0, GF::Types<8, 0b100011101, uint8_t>> rs;
		uint8_t code[255];
		uint8_t target[255] = { 1, 126, 147, 48, 155, 224, 3, 157, 29, 226, 40, 114, 61, 30, 244, 75 };
		for (int i = 0; i < 239; ++i)
			target[16+i] = code[16+i] = i + 1;
		test(rs, code, target);
	}
}

