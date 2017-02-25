/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef REED_SOLOMON_HH
#define REED_SOLOMON_HH

#include "galois_field.hh"
#include "correction.hh"

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
