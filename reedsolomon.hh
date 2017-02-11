/*
reedsolomon - Reed-Solomon error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef REEDSOLOMON_HH
#define REEDSOLOMON_HH

#include "galoisfield.hh"

template <int NR, int FCR, typename GF>
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
		// $generator = \prod_{i=0}^{NR}(x-pe^{FCR+i})$
		ValueType tmp[NR+1];
		IndexType root(FCR), pe(1);
		for (int i = 0; i < NR; ++i) {
			tmp[i] = ValueType(1);
			for (int j = i; j > 0; --j)
				tmp[j] = fma(root, tmp[j], tmp[j-1]);
			tmp[0] *= root;
			root *= pe;
		}
		tmp[NR] = ValueType(1);
#ifndef NDEBUG
		std::cout << "generator = ";
		for (int i = NR; i > 0; --i) {
			if (!tmp[i])
				continue;
			if (tmp[i] != ValueType(1))
				std::cout << (int)tmp[i] << "*";
			std::cout << "x";
			if (i != 1)
				std::cout << "^" << i;
			std::cout << " + ";
		}
		std::cout << (int)tmp[0] << std::endl;
#endif
		for (int i = 0; i <= NR; ++i)
			generator[i] = index(tmp[i]);
	}
	void encode(ValueType *code)
	{
		// $code = data * x^{NR} + (data * x^{NR}) \mod{generator}$
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
	int Berlekamp_Massey_algorithm(ValueType *s, ValueType *C, IndexType *erasures = 0, int count = 0)
	{
		C[0] = ValueType(1);
		for (int i = 1; i <= NR; ++i)
			C[i] = ValueType(0);
		// $C = \prod_{i=0}^{count}(1-x\,pe^{N-1-erasures_i})$
		if (count)
			C[1] = value(IndexType(N-1) / erasures[0]);
		for (int i = 1; i < count; ++i) {
			IndexType tmp(IndexType(N-1) / erasures[i]);
			for (int j = i; j >= 0; --j)
				C[j+1] += tmp * C[j];
		}
		ValueType B[NR+1];
		for (int i = 0; i <= NR; ++i)
			B[i] = C[i];
		int L = count;
		for (int n = count, m = 1; n < NR; ++n) {
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
				if (2 * L <= n + count) {
					L = n + count + 1 - L;
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
	int Chien_search(ValueType *locator, int locator_degree, IndexType *locations)
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
				locations[count++] = IndexType(i);
		}
		return count;
	}
	int compute_evaluator(ValueType *syndromes, ValueType *locator, int locator_degree, ValueType *evaluator)
	{
		// $evaluator = (syndromes * locator) \bmod{x^{NR}}$
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
	void compute_magnitudes(ValueType *locator, IndexType *locations, int count, ValueType *evaluator, int evaluator_degree, ValueType *magnitudes)
	{
		// $magnitude = root^{FCR-1} * \frac{evaluator(root)}{locator'(root)}$
		for (int i = 0; i < count; ++i) {
			IndexType root(locations[i] * IndexType(1)), tmp(root);
			ValueType eval(evaluator[0]);
			for (int j = 1; j <= evaluator_degree; ++j) {
				eval += evaluator[j] * tmp;
				tmp *= root;
			}
			if (!eval) {
				magnitudes[i] = ValueType(0);
				continue;
			}
			ValueType deriv(locator[1]);
			IndexType root2(root * root), tmp2(root2);
			for (int j = 3; j <= count; j += 2) {
				deriv += locator[j] * tmp2;
				tmp2 *= root2;
			}
			IndexType magnitude(index(eval) / index(deriv));
			if (FCR == 0)
				magnitude /= root;
			if (FCR > 1)
				for (int j = 1; j < FCR; ++j)
					magnitude *= root;
			magnitudes[i] = value(magnitude);
		}
	}
	int Forney_algorithm(ValueType *syndromes, ValueType *locator, IndexType *locations, int count, ValueType *evaluator, ValueType *magnitudes)
	{
		int evaluator_degree = compute_evaluator(syndromes, locator, count, evaluator);
		compute_magnitudes(locator, locations, count, evaluator, evaluator_degree, magnitudes);
		return evaluator_degree;
	}
	int correct(ValueType *code, ValueType *syndromes, IndexType *erasures = 0, int erasures_count = 0)
	{
		ValueType locator[NR+1];
		int locator_degree = Berlekamp_Massey_algorithm(syndromes, locator, erasures, erasures_count);
		assert(locator_degree);
		assert(locator_degree <= NR);
		assert(locator[0] == ValueType(1));
		while (!locator[locator_degree])
			if (--locator_degree < 0)
				return -1;
		IndexType locations[locator_degree];
		int count;
		if (locator_degree == 1) {
			count = 1;
			locations[0] = (index(locator[0]) / index(locator[1])) / IndexType(1);
		} else if (locator_degree == 2) {
			if (!locator[1] || !locator[0])
				return -1;
			ValueType a(locator[2]), b(locator[1]), c(locator[0]);
			ValueType ba(b/a), R(Artin_Schreier_imap(a*c/(b*b)));
			if (!R)
				return -1;
			count = 2;
			locations[0] = index(ba * R) / IndexType(1);
			locations[1] = index(ba * R + ba) / IndexType(1);
		} else {
			count = Chien_search(locator, locator_degree, locations);
			if (count < locator_degree)
				return -1;
		}
		ValueType evaluator[NR];
		ValueType magnitudes[count];
		int evaluator_degree = Forney_algorithm(syndromes, locator, locations, count, evaluator, magnitudes);
		for (int i = 0; i < count; ++i)
			code[(int)locations[i]] += magnitudes[i];
		int corrections_count = 0;
		for (int i = 0; i < count; ++i)
			corrections_count += !!magnitudes[i];
#ifdef NDEBUG
		(void)evaluator_degree;
#else
		static int init;
		if (!init) {
			init = 1;
			std::cout << "syndromes =";
			for (int i = 0; i < NR; ++i)
				std::cout << " " << (int)syndromes[i];
			std::cout << std::endl;
			std::cout << "locator = ";
			for (int i = NR; i > 0; --i) {
				if (!locator[i])
					continue;
				if (locator[i] != ValueType(1))
					std::cout << (int)locator[i] << "*";
				std::cout << "x";
				if (i != 1)
					std::cout << "^" << i;
				std::cout << " + ";
			}
			std::cout << (int)locator[0] << std::endl;
			std::cout << "locations =";
			for (int i = 0; i < count; ++i)
				std::cout << " " << (int)locations[i];
			std::cout << std::endl;
			std::cout << "evaluator = ";
			for (int i = evaluator_degree; i > 0; --i) {
				if (!evaluator[i])
					continue;
				if (evaluator[i] != ValueType(1))
					std::cout << (int)evaluator[i] << "*";
				std::cout << "x";
				if (i != 1)
					std::cout << "^" << i;
				if (i != 1 || evaluator[0])
					std::cout << " + ";
			}
			if (evaluator[0])
				std::cout << (int)evaluator[0];
			std::cout << std::endl;
			std::cout << "magnitudes =";
			for (int i = 0; i < count; ++i)
				std::cout << " " << (int)magnitudes[i];
			std::cout << std::endl;
		}
#endif
		return corrections_count;
	}
	int compute_syndromes(ValueType *code, ValueType *syndromes)
	{
		// $syndromes_i = code(pe^{FCR+i})$
		for (int i = 0; i < NR; ++i)
			syndromes[i] = code[0];
		for (int j = 1; j < N; ++j) {
			IndexType root(FCR), pe(1);
			for (int i = 0; i < NR; ++i) {
				syndromes[i] = fma(root, syndromes[i], code[j]);
				root *= pe;
			}
		}
		int nonzero = 0;
		for (int i = 0; i < NR; ++i)
			nonzero += !!syndromes[i];
		return nonzero;
	}
	int decode(ValueType *code, IndexType *erasures = 0, int erasures_count = 0)
	{
#if 0
		for (int i = 0; i < erasures_count; ++i)
			code[(int)erasures[i]] = ValueType(0);
#endif
		ValueType syndromes[NR];
		if (compute_syndromes(code, syndromes))
			return correct(code, syndromes, erasures, erasures_count);
		return 0;
	}
	void encode(value_type *code)
	{
		encode(reinterpret_cast<ValueType *>(code));
	}
	int decode(value_type *code, value_type *erasures = 0, int erasures_count = 0)
	{
		return decode(reinterpret_cast<ValueType *>(code), reinterpret_cast<IndexType *>(erasures), erasures_count);
	}
	int compute_syndromes(value_type *code, value_type *syndromes)
	{
		return compute_syndromes(reinterpret_cast<ValueType *>(code), reinterpret_cast<ValueType *>(syndromes));
	}
};

#endif
