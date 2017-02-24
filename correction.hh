/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef CORRECTION_HH
#define CORRECTION_HH

#include "galois_field.hh"
#include "berlekamp_massey.hh"
#include "chien.hh"
#include "forney.hh"

template <int NR, int FCR, typename GF>
struct Correction
{
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, K = N - NR;
	static int algorithm(ValueType *code, ValueType *syndromes, IndexType *erasures = 0, int erasures_count = 0)
	{
		assert(0 <= erasures_count && erasures_count <= NR);
		ValueType locator[NR+1];
		locator[0] = ValueType(1);
		for (int i = 1; i <= NR; ++i)
			locator[i] = ValueType(0);
		// $locator = \prod_{i=0}^{count}(1-x\,pe^{N-1-erasures_i})$
		if (erasures_count)
			locator[1] = value(IndexType(N-1) / erasures[0]);
		for (int i = 1; i < erasures_count; ++i) {
			IndexType tmp(IndexType(N-1) / erasures[i]);
			for (int j = i; j >= 0; --j)
				locator[j+1] += tmp * locator[j];
		}
		int locator_degree = BerlekampMassey<NR, GF>::algorithm(syndromes, locator, erasures_count);
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
			count = Chien<NR, GF>::search(locator, locator_degree, locations);
			if (count < locator_degree)
				return -1;
		}
		ValueType evaluator[NR];
		ValueType magnitudes[count];
		int evaluator_degree = Forney<NR, FCR, GF>::algorithm(syndromes, locator, locations, count, evaluator, magnitudes);
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
};

#endif
