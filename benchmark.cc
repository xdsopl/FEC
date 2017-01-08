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
		if (!value(rs.generator[i]))
			continue;
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
	assert(!rs.decode(code));

	int blocks = (8 * data.size() + M * rs.K - 1) / (M * rs.K);
	TYPE *tmp = new TYPE[rs.N * blocks];
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
			assert(!rs.decode(tmp + i * rs.N));
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int bytes = (rs.N * blocks * M) / 8;
		int mbs = bytes / (msec.count() * 1000);
		std::cout << "checking of " << bytes << " encoded bytes (" << blocks << " blocks) took " << msec.count() << " milliseconds (" << mbs << "MB/s)." << std::endl;
	}
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<int> bit_dist(0, M-1), pos_dist(0, rs.N-1);
	auto rnd_bit = std::bind(bit_dist, generator);
	auto rnd_pos = std::bind(pos_dist, generator);
	int corrupt = 0;
	const int places = NR/2;
	for (int i = 0; i < blocks; ++i) {
		int pos[places];
		for (int j = 0; j < places; ++j) {
			pos[j] = rnd_pos();
			for (int k = 0; k < j;) {
				if (pos[k++] == pos[j]) {
					pos[j] = rnd_pos();
					k = 0;
				}
			}
			tmp[i * rs.N + pos[j]] ^= 1 << rnd_bit();
			++corrupt;
		}
	}
	{
		auto start = std::chrono::system_clock::now();
		int errors = 0;
		for (int i = 0; i < blocks; ++i)
			errors += rs.decode(tmp + i * rs.N);
		assert(corrupt == errors);
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int bytes = (rs.N * blocks * M) / 8;
		int mbs = bytes / (msec.count() * 1000);
		std::cout << "decoding of " << bytes << " encoded bytes with " << errors << " errors took " << msec.count() << " milliseconds (" << mbs << "MB/s)." << std::endl;
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
	{
		ReedSolomon<64, 0, GF::Types<16, 0b10001000000001011, uint16_t>> rs;
		uint16_t code[65535], target[65535];
		for (int i = 0; i < 65471; ++i)
			target[i] = code[i] = i + 1;
		uint16_t parity[64] = { 1466, 53891, 11202, 5869, 43221, 51318, 62552, 36624, 61542, 53556, 25812, 11238, 38013, 14282, 50002, 51910, 17214, 17728, 12240, 31701, 53101, 23446, 25703, 875, 48946, 154, 24094, 25945, 24127, 7224, 50886, 25737, 36097, 13487, 3169, 18228, 10429, 17419, 52734, 24366, 36761, 35318, 25193, 30701, 42547, 61112, 33698, 62580, 21986, 37239, 48458, 11335, 49801, 17076, 35041, 37113, 37451, 31305, 5758, 29560, 19653, 54009, 52622, 50101 };
		for (int i = 0; i < 64; ++i)
			target[65471+i] = parity[i];
		test("FUN RS(65535, 65471) T=32", rs, code, target, data);
	}
}

