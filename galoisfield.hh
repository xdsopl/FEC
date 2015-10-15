/*
reedsolomon - Reed-Solomon error correction
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef GALOISFIELD_HH
#define GALOISFIELD_HH

#include <cassert>

template <int M, int POLY, typename TYPE>
class GaloisField
{
public:
	typedef TYPE value_type;
	static const int Q = 1 << M, N = Q - 1;
	value_type log_table[Q];
	value_type exp_table[Q];
	GaloisField()
	{
		assert(M <= 8 * sizeof(TYPE));
		assert(Q == (POLY & ~N));
		log_table[0] = N;
		exp_table[N] = 0;
		int a = 1;
		for (int i = 0; i < N; ++i) {
			log_table[a] = i;
			exp_table[i] = a;
			a <<= 1;
			a = a & Q ? a ^ POLY : a;
		}
		assert(1 == a);
	}
	TYPE log(TYPE a) { return log_table[a]; }
	TYPE exp(TYPE a) { return exp_table[a]; }
	TYPE add(TYPE a, TYPE b) { return a ^ b; }
	TYPE mul(TYPE a, TYPE b)
	{
		int tmp = (int)log_table[a] + (int)log_table[b];
		return (!a || !b) ? 0 : exp_table[tmp < N ? tmp : tmp - N];

	}
	TYPE fma(TYPE a, TYPE b, TYPE c)
	{
		int tmp = (int)log_table[a] + (int)log_table[b];
		return (!a || !b) ? c : c ^ exp_table[tmp < N ? tmp : tmp - N];

	}
	TYPE div(TYPE a, TYPE b)
	{
		int tmp = (int)log_table[a] - (int)log_table[b];
		return !b ? N : !a ? 0 : exp_table[tmp >= 0 ? tmp : tmp + N];

	}
	TYPE rcp(TYPE a) { return !a ? N : exp_table[N - log_table[a]]; }
};

#endif
