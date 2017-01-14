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
#ifndef NDEBUG
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
	int Chien_search(ValueType *locator, int locator_degree, ValueType *locations)
	{
		ValueType tmp[locator_degree+1];
		for (int i = 0; i <= locator_degree; ++i)
			tmp[i] = locator[i];
		int count = 0;
		for (int i = 0; i < N; ++i) {
			ValueType sum(tmp[0]);
			for (int j = 1; j <= locator_degree; ++j)
				sum += tmp[j] *= IndexType(j);
			if (!sum)
				locations[count++] = ValueType(i);
		}
		return count;
	}
	int compute_evaluator(ValueType *syndromes, ValueType *locator, int locator_degree, ValueType *evaluator)
	{
		// (syndromes * locator) mod x^NR
		int tmp = std::min(locator_degree, NR-1);
		int degree = -1;
		for (int i = 0; i <= tmp; ++i) {
			evaluator[i] = syndromes[i] * locator[0];
			for (int j = 1; j <= i; ++j)
				evaluator[i] += syndromes[i-j] * locator[j];
			if (evaluator[i])
				degree = i;
		}
		return degree;
	}
	void compute_magnitudes(ValueType *locator, ValueType *locations, int count, ValueType *evaluator, int evaluator_degree, ValueType *magnitudes)
	{
		for (int i = 0; i < count; ++i) {
			IndexType root(IndexType(locations[i].v) * IndexType(1)), tmp(root);
			ValueType numerator(evaluator[0]);
			for (int j = 1; j <= evaluator_degree; ++j) {
				numerator += evaluator[j] * tmp;
				tmp *= root;
			}
			for (int i = 0; i < FR; ++i)
				numerator *= root;
			ValueType denominator(locator[1]);
			IndexType root2(root * root), tmp2(root2);
			for (int j = 3; j <= count; j += 2) {
				denominator += locator[j] * tmp2;
				tmp2 *= root2;
			}
			denominator *= root;
			magnitudes[i] = numerator / denominator;
		}
	}
	int correct(ValueType *code, ValueType *syndromes)
	{
		ValueType locator[NR+1];
		int errors = Berlekamp_Massey_algorithm(syndromes, locator);
		assert(locator[0] == ValueType(1));
		int locator_degree = NR;
		while (!locator[locator_degree])
			if (--locator_degree < 0)
				return -1;
		ValueType locations[locator_degree];
		int locations_count = Chien_search(locator, locator_degree, locations);
		if (locations_count < locator_degree)
			return -1;
		// Forney algorithm
		ValueType evaluator[NR];
		int evaluator_degree = compute_evaluator(syndromes, locator, locator_degree, evaluator);
		ValueType magnitudes[locations_count];
		compute_magnitudes(locator, locations, locations_count, evaluator, evaluator_degree, magnitudes);
		for (int i = 0; i < locations_count; ++i)
			code[locations[i].v] += magnitudes[i];
#ifndef NDEBUG
		static int init;
		if (!init) {
			init = 1;
			std::cout << "s(x) = ";
			for (int i = NR-1; i > 0; --i) {
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
			for (int i = 0; i < locations_count; ++i)
				std::cout << " " << (int)locations[i].v;
			std::cout << std::endl;
			std::cout << "w(z) = ";
			for (int i = evaluator_degree; i > 0; --i) {
				if (!evaluator[i].v)
					continue;
				if (evaluator[i].v != 1)
					std::cout << (int)evaluator[i].v << "*";
				std::cout << "x";
				if (i != 1)
					std::cout << "^" << i;
				if (i != 1 || evaluator[0].v)
					std::cout << " + ";
			}
			if (evaluator[0].v)
				std::cout << (int)evaluator[0].v;
			std::cout << std::endl;
			std::cout << "magnitudes =";
			for (int i = 0; i < locations_count; ++i)
				std::cout << " " << (int)magnitudes[i].v;
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
