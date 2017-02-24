/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef FORNEY_HH
#define FORNEY_HH

#include "galois_field.hh"

template <int NR, int FCR, typename GF>
struct Forney
{
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, K = N - NR;
	static int compute_evaluator(ValueType *syndromes, ValueType *locator, int locator_degree, ValueType *evaluator)
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
	static void compute_magnitudes(ValueType *locator, IndexType *locations, int count, ValueType *evaluator, int evaluator_degree, ValueType *magnitudes)
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
	static int algorithm(ValueType *syndromes, ValueType *locator, IndexType *locations, int count, ValueType *evaluator, ValueType *magnitudes)
	{
		int evaluator_degree = compute_evaluator(syndromes, locator, count, evaluator);
		compute_magnitudes(locator, locations, count, evaluator, evaluator_degree, magnitudes);
		return evaluator_degree;
	}
};

#endif
