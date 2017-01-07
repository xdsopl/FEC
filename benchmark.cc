/*
reedsolomon - Reed-Solomon error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <iostream>
#include <cassert>
#include <random>
#include <chrono>
#include <functional>
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
void test(std::string name, ReedSolomon<NR, FR, GF::Types<M, P, TYPE>> &rs, TYPE *code, TYPE *target, std::vector<uint8_t> &data)
{
	std::cout << "testing: " << name << std::endl;
#if 0
	std::cout << "g(x) = x^" << NR << " + ";
	for (int i = NR-1; i > 0; --i) {
		if (value(rs.generator[i]) != 1)
			std::cout << (int)value(rs.generator[i]) << "*";
		std::cout << "x";
		if (i != 1)
			std::cout << "^" << i;
		std::cout << " + ";
	}
	std::cout << (int)value(rs.generator[0]) << std::endl;
#endif
	rs.encode(code);
	for (int i = 0; i < rs.N; ++i)
		assert(code[i] == target[i]);
	//print_table(code + rs.K, "parity", NR);
	for (int i = rs.N-1, n = 0; i >= 0 && n < 2; --i, ++n)
		code[i] ^= 1;
	assert(rs.decode(code));

	int blocks = (8 * data.size() + M * rs.K - 1) / (M * rs.K);
	uint8_t *tmp = new uint8_t[rs.N * blocks];
	unsigned acc = 0, bit = 0, pos = 0;
	for (unsigned byte : data) {
		acc |= byte << bit;
		bit += 8;
		while (bit >= M) {
			bit -= M;
			tmp[pos++] = rs.N & acc;
			acc >>= M;
			if (pos % rs.N >= rs.K)
				pos += NR;
		}
	}
	{
		auto start = std::chrono::system_clock::now();
		for (int i = 0; i < blocks; ++i)
			rs.encode(tmp + i * rs.N);
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int mbs = data.size() / (msec.count() * 1000);
		std::cout << "encoding of " << data.size() << " random bytes (" << blocks << " blocks) took " << msec.count() << " milliseconds (" << mbs << "MB/s)." << std::endl;
	}
	{
		auto start = std::chrono::system_clock::now();
		for (int i = 0; i < blocks; ++i)
			assert(rs.decode(tmp + i * rs.N));
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int mbs = data.size() / (msec.count() * 1000);
		std::cout << "decoding of " << data.size() << " random bytes (" << blocks << " blocks) took " << msec.count() << " milliseconds (" << mbs << "MB/s)." << std::endl;
	}
	delete[] tmp;
}

int main()
{
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<uint8_t> distribution(0, 255);
	std::vector<uint8_t> data(10000000);
	std::generate(data.begin(), data.end(), std::bind(distribution, generator));
	{
		ReedSolomon<4, 0, GF::Types<4, 0b10011, uint8_t>> rs;
		uint8_t code[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
		uint8_t target[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 3, 3, 12, 12 };
		test("BBC WHP031 RS(15, 11) T=2", rs, code, target, data);
	}
	{
		ReedSolomon<16, 0, GF::Types<8, 0b100011101, uint8_t>> rs;
		uint8_t code[255], target[255];
		for (int i = 0; i < 239; ++i)
			target[i] = code[i] = i + 1;
		uint8_t parity[16] = { 1, 126, 147, 48, 155, 224, 3, 157, 29, 226, 40, 114, 61, 30, 244, 75 };
		for (int i = 0; i < 16; ++i)
			target[239+i] = parity[i];
		test("DVB-T RS(255, 239) T=8", rs, code, target, data);
	}
}

