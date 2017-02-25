/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef BOSE_CHAUDHURI_HOCQUENGHEM_HH
#define BOSE_CHAUDHURI_HOCQUENGHEM_HH

#include <initializer_list>
#include "galois_field.hh"
#include "correction.hh"

template <int NR, int FCR, int K, typename GF>
class BoseChaudhuriHocquenghem
{
public:
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, NP = N - K;
	ValueType generator[NP+1];
	BoseChaudhuriHocquenghem(std::initializer_list<int> minimal_polynomials)
	{
		// $generator(x) = \prod_i(minpoly_i(x))$
		int generator_degree = 1;
		generator[0] = ValueType(1);
		for (int i = 1; i <= NP; ++i)
			generator[i] = ValueType(0);
		for (auto m: minimal_polynomials) {
			assert(0 < m && m < 1<<(GF::M+1));
			int m_degree = GF::M;
			while (!(m>>m_degree))
				--m_degree;
			assert(generator_degree + m_degree <= NP + 1);
			for (int i = generator_degree; i >= 0; --i) {
				if (!generator[i])
					continue;
				generator[i] = ValueType(m&1);
				for (int j = 1; j <= m_degree; ++j)
					generator[i+j] += ValueType((m>>j)&1);
			}
			generator_degree += m_degree;
		}
		assert(generator_degree == NP + 1);
#ifndef NDEBUG
		IndexType root(FCR), pe(1);
		for (int i = 0; i < NR; ++i) {
			ValueType tmp(generator[NP]);
			for (int j = 1; j <= NP; ++j)
				tmp = fma(root, tmp, generator[NP-j]);
			assert(!tmp);
			root *= pe;
		}
		std::cout << "generator =";
		for (int i = 0; i <= NP; ++i)
			std::cout << " " << (int)generator[i];
		std::cout << std::endl;
#endif
	}
	void encode(ValueType *code)
	{
		// $code = data * x^{NP} + (data * x^{NP}) \mod{generator}$
		for (int i = 0; i < NP; ++i)
			code[K+i] = ValueType(0);
		for (int i = 0; i < K; ++i) {
			if (code[i] != code[K]) {
				for (int j = 1; j < NP; ++j)
					code[K+j-1] = generator[NP-j] + code[K+j];
				code[N-1] = generator[0];
			} else {
				for (int j = 1; j < NP; ++j)
					code[K+j-1] = code[K+j];
				code[N-1] = ValueType(0);
			}
		}
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
		assert(0 <= erasures_count && erasures_count <= NR);
#if 0
		for (int i = 0; i < erasures_count; ++i)
			code[(int)erasures[i]] = ValueType(0);
#endif
		ValueType syndromes[NR];
		if (!compute_syndromes(code, syndromes))
			return 0;
		IndexType locations[NR];
		ValueType magnitudes[NR];
		int count = Correction<NR, FCR, GF>::algorithm(syndromes, locations, magnitudes, erasures, erasures_count);
		if (count <= 0)
			return count;
		for (int i = 0; i < count; ++i)
			if (1 < (int)magnitudes[i])
				return -1;
		for (int i = 0; i < count; ++i)
			code[(int)locations[i]] += magnitudes[i];
		int corrections_count = 0;
		for (int i = 0; i < count; ++i)
			corrections_count += !!magnitudes[i];
		return corrections_count;
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
