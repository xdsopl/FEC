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
#include <algorithm>
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

template <int NR, int FCR, int M, int P, typename TYPE>
void test(std::string name, ReedSolomon<NR, FCR, GF::Types<M, P, TYPE>> &rs, TYPE *code, TYPE *target, std::vector<uint8_t> &data)
{
	std::cout << "testing: " << name << std::endl;

	{
		rs.encode(code);
		bool error = false;
		for (int i = 0; i < rs.N; ++i)
			error |= code[i] != target[i];
		if (error)
			std::cout << "encoder error!" << std::endl;
		assert(!error);
		//print_table(code + rs.K, "parity", NR);
		error = rs.decode(code);
		if (error)
			std::cout << "decoder error!" << std::endl;
		assert(!error);
		int corrupt = 0;
		for (int i = 0; i < rs.N && i < 1; ++i, ++corrupt)
			code[i] ^= 1;
		int corrected = rs.decode(code);
		if (corrupt != corrected)
			std::cout << "decoder error: expected " << corrupt << " but got " << corrected << std::endl;
		assert(corrupt == corrected);
		if (corrected >= 0 && rs.decode(code)) {
			std::cout << "decoder error: result of correction is not a codeword!" << std::endl;
			assert(false);
		}
		error = false;
		for (int i = 0; i < rs.N; ++i)
			error |= code[i] != target[i];
		if (error)
			std::cout << "decoder error: code doesnt match target" << std::endl;
		assert(!error);
	}

	int blocks = (8 * data.size() + M * rs.K - 1) / (M * rs.K);
	TYPE *coded = new TYPE[rs.N * blocks];
	{
		unsigned acc = 0, bit = 0, pos = 0;
		for (unsigned byte : data) {
			acc |= byte << bit;
			bit += 8;
			while (bit >= M) {
				bit -= M;
				coded[pos++] = rs.N & acc;
				acc >>= M;
				if (pos % rs.N >= rs.K)
					pos += NR;
			}
		}
	}
	{
		auto start = std::chrono::system_clock::now();
		for (int i = 0; i < blocks; ++i)
			rs.encode(coded + i * rs.N);
		auto end = std::chrono::system_clock::now();
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int mbs = (data.size() + msec.count() / 2) / msec.count();
		int bytes = (rs.N * blocks * M) / 8;
		int redundancy = (100 * (bytes-data.size()) + data.size() / 2) / data.size();
		std::cout << "encoding of " << data.size() << " random bytes into " << bytes << " codeword bytes (" << redundancy << "% redundancy) in " << blocks << " blocks took " << msec.count() << " milliseconds (" << mbs << "KB/s)." << std::endl;
	}
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<int> bit_dist(0, M-1), pos_dist(0, rs.N-1);
	auto rnd_bit = std::bind(bit_dist, generator);
	auto rnd_pos = std::bind(pos_dist, generator);
	std::vector<uint8_t> recovered(data.size());
	TYPE *tmp = new TYPE[rs.N * blocks];
	for (int i = 0; i < rs.N * blocks; ++i)
		tmp[i] = coded[i];
	for (int places = 0; places <= NR; ++places) {
		int corrupt = 0;
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
		int corrected = 0, wrong = 0;
		auto start = std::chrono::system_clock::now();
		for (int i = 0; i < blocks; ++i) {
			int result = rs.decode(tmp + i * rs.N);
			if (places > NR/2 && result >= 0)
				for (int j = i * rs.N; j < (i + 1) * rs.N; ++j)
					wrong += coded[j] != tmp[j];
			corrected += result;
		}
		auto end = std::chrono::system_clock::now();
		if (corrupt != corrected)
			std::cout << "decoder error: expected " << corrupt << " corrected errors but got " << corrected << " and " << wrong << " wrong corrections." << std::endl;
		assert(places > NR/2 || corrupt == corrected);
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		int bytes = (rs.N * blocks * M) / 8;
		int mbs = (bytes + msec.count() / 2) / msec.count();
		std::cout << "decoding of " << blocks << " blocks with " << places << " errors per block took " << msec.count() << " milliseconds (" << mbs << "KB/s)." << std::endl;
		unsigned acc = 0, bit = 0, pos = 0;
		for (uint8_t &byte: recovered) {
			while (bit < 8) {
				acc |= (unsigned)tmp[pos++] << bit;
				bit += M;
				if (pos % rs.N >= rs.K)
					pos += NR;
			}
			bit -= 8;
			byte = 255 & acc;
			acc >>= 8;
		}
		if (data != recovered) {
			std::cout << "decoder error: data could not be recovered from corruption" << std::endl;
			assert(places > NR/2);
		}
	}
	delete[] tmp;
	delete[] coded;
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
		ReedSolomon<64, 1, GF::Types<16, 0b10001000000001011, uint16_t>> rs;
		uint16_t code[65535], target[65535];
		for (int i = 0; i < 65471; ++i)
			target[i] = code[i] = i + 1;
		uint16_t parity[64] = { 25271, 26303, 22052, 31318, 31233, 6076, 40148, 29468, 47507, 32655, 12404, 13265, 23901, 38403, 50967, 50433, 40818, 226, 62296, 23636, 56393, 12952, 11476, 44416, 518, 50014, 10037, 57582, 33421, 42654, 54025, 7157, 4826, 52148, 17167, 23294, 6427, 40953, 11168, 35305, 18209, 1868, 39971, 54928, 27566, 1424, 4846, 25347, 34710, 42190, 56452, 21859, 49805, 28028, 41657, 25756, 22014, 24479, 28758, 17438, 12976, 61743, 46735, 1557 };
		for (int i = 0; i < 64; ++i)
			target[65471+i] = parity[i];
		test("FUN RS(65535, 65471) T=32", rs, code, target, data);
	}
}

