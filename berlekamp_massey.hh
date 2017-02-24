/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef BERLEKAMP_MASSEY_HH
#define BERLEKAMP_MASSEY_HH

#include "galois_field.hh"

template <int NR, typename GF>
struct BerlekampMassey
{
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, K = N - NR;
	static int algorithm(ValueType *s, ValueType *C, int count = 0)
	{
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
};

#endif
