/*
reedsolomon - Reed-Solomon error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef REEDSOLOMON_HH
#define REEDSOLOMON_HH

#include "galoisfield.hh"

template <int NR, int FR, typename GF>
class ReedSolomon
{
public:
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, K = N - NR;
	IndexType generator[NR+1];
	ReedSolomon()
	{
		ValueType tmp[NR+1];
		IndexType root(FR), pe(1);
		for (int i = 0; i < NR; ++i) {
			tmp[i] = ValueType(1);
			for (int j = i; j > 0; --j)
				tmp[j] = fma(root, tmp[j], tmp[j-1]);
			tmp[0] *= root;
			root *= pe;
		}
		tmp[NR] = ValueType(1);
#if 0
		std::cout << "g(x) = ";
		for (int i = NR; i > 0; --i) {
			if (!tmp[i].v)
				continue;
			if (tmp[i].v != 1)
				std::cout << (int)tmp[i].v << "*";
			std::cout << "x";
			if (i != 1)
				std::cout << "^" << i;
			std::cout << " + ";
		}
		std::cout << (int)tmp[0].v << std::endl;
#endif
		for (int i = 0; i <= NR; ++i)
			generator[i] = index(tmp[i]);
	}
	void encode(ValueType *code)
	{
		for (int i = 0; i < NR; ++i)
			code[K+i] = ValueType(0);
		for (int i = 0; i < K; ++i) {
			ValueType feedback = code[i] + code[K];
			if (feedback) {
				IndexType fb = index(feedback);
				for (int j = 1; j < NR; ++j)
					code[K+j-1] = fma(fb, generator[NR-j], code[K+j]);
				code[N-1] = value(generator[0] * fb);
			} else {
				for (int j = 1; j < NR; ++j)
					code[K+j-1] = code[K+j];
				code[N-1] = ValueType(0);
			}
		}
	}
	int Berlekamp_Massey_algorithm(ValueType *s, ValueType *C)
	{
		ValueType B[NR+1];
		B[0] = C[0] = ValueType(1);
		for (int i = 1; i <= NR; ++i)
			B[i] = C[i] = ValueType(0);
		int L = 0;
		for (int n = 0, m = 1; n < NR; ++n) {
			ValueType d(s[n]);
			for (int i = 1; i <= L; ++i)
				d += C[i] * s[n-i];
			if (!d) {
				++m;
			} else {
				ValueType T[NR+1];
				for (int i = 0; i < m; ++i)
					T[i] = C[i];
				for (int i = m; i <= NR; ++i)
					T[i] = fma(d, B[i-m], C[i]);
				if (2 * L <= n) {
					L = n + 1 - L;
					for (int i = 0; i <= NR; ++i)
						B[i] = C[i] / d;
					m = 1;
				} else {
					++m;
				}
				for (int i = 0; i <= NR; ++i)
					C[i] = T[i];
			}
		}
		return L;
	}
	int Chien_search(ValueType *locator, ValueType *locations)
	{
		int degree = NR;
		while (!locator[degree])
			if (--degree < 0)
				return degree;
		ValueType tmp[degree+1];
		for (int i = 0; i <= degree; ++i)
			tmp[i] = locator[i];
		int count = 0;
		for (int i = 0; i < N; ++i) {
			ValueType sum(tmp[0]);
			for (int j = 1; j <= degree; ++j)
				sum += tmp[j] *= IndexType(j);
			if (!sum)
				locations[count++] = ValueType(i);
		}
		return count;
	}
	int correct(ValueType *code, ValueType *syndromes)
	{
		ValueType locator[NR+1];
		int errors = Berlekamp_Massey_algorithm(syndromes, locator);
		ValueType locations[NR];
		int count = Chien_search(locator, locations);
#if 0
		static int init;
		if (!init) {
			init = 1;
			std::cout << "s(x) = ";
			for (int i = NR; i > 0; --i) {
				if (!syndromes[i].v)
					continue;
				if (syndromes[i].v != 1)
					std::cout << (int)syndromes[i].v << "*";
				std::cout << "x";
				if (i != 1)
					std::cout << "^" << i;
				std::cout << " + ";
			}
			std::cout << (int)syndromes[0].v << std::endl;
			std::cout << "C(x) = ";
			for (int i = NR; i > 0; --i) {
				if (!locator[i].v)
					continue;
				if (locator[i].v != 1)
					std::cout << (int)locator[i].v << "*";
				std::cout << "x";
				if (i != 1)
					std::cout << "^" << i;
				std::cout << " + ";
			}
			std::cout << (int)locator[0].v << std::endl;
			std::cout << "locations =";
			for (int i = 0; i < count; ++i)
				std::cout << " " << (int)locations[i].v;
			std::cout << std::endl;
		}
#endif
		return errors;
	}
	int decode(ValueType *code)
	{
		ValueType syndromes[NR];
		for (int i = 0; i < NR; ++i)
			syndromes[i] = code[0];
		for (int j = 1; j < N; ++j) {
			IndexType root(FR), pe(1);
			for (int i = 0; i < NR; ++i) {
				syndromes[i] = fma(root, syndromes[i], code[j]);
				root *= pe;
			}
		}
		for (int i = 0; i < NR; ++i)
			if (syndromes[i])
				return correct(code, syndromes);
		return 0;
	}
	void encode(value_type *code)
	{
		encode(reinterpret_cast<ValueType *>(code));
	}
	int decode(value_type *code)
	{
		return decode(reinterpret_cast<ValueType *>(code));
	}
};

#endif
