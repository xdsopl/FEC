/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef CHIEN_HH
#define CHIEN_HH

#include "galois_field.hh"

template <int NR, typename GF>
struct Chien
{
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, K = N - NR;
	static int search(ValueType *locator, int locator_degree, IndexType *locations)
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
};

#endif
