/*
reedsolomon - Reed-Solomon error correction
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef GALOISFIELD_HH
#define GALOISFIELD_HH

template <int M, int POLY, typename TYPE>
class GaloisField
{
public:
	typedef TYPE value_type;
	static const int Q = 1 << M, N = Q - 1;
	static_assert(M <= 8 * sizeof(TYPE), "TYPE not wide enough");
	static_assert(Q == (POLY & ~N), "POLY not of degree Q");
	GaloisField()
	{
		if (log(0))
			return;
		log(0) = N;
		exp(N) = 0;
		TYPE a = 1;
		for (int i = 0; i < N; ++i, a = next(a))
			log(exp(i) = a) = i;
		//assert(1 == a);
	}
	TYPE next(TYPE a)
	{
		return a & (TYPE)(Q >> 1) ? (a << 1) ^ (TYPE)POLY : a << 1;
	}
	static TYPE &log(TYPE a)
	{
		static value_type table[Q];
		return table[a];
	}
	static TYPE &exp(TYPE a)
	{
		static value_type table[Q];
		return table[a];
	}
	TYPE add(TYPE a, TYPE b) { return a ^ b; }
	TYPE mul(TYPE a, TYPE b)
	{
		TYPE loga = log(a);
		TYPE logb = log(b);
		TYPE tmp = loga + logb;
		tmp = (TYPE)N - loga <= logb ? tmp - (TYPE)N : tmp;
		return (!a || !b) ? (TYPE)0 : exp(tmp);

	}
	TYPE fma(TYPE a, TYPE b, TYPE c)
	{
		TYPE loga = log(a);
		TYPE logb = log(b);
		TYPE tmp = loga + logb;
		tmp = (TYPE)N - loga <= logb ? tmp - (TYPE)N : tmp;
		return (!a || !b) ? c : c ^ exp(tmp);

	}
	TYPE div(TYPE a, TYPE b)
	{
		TYPE loga = log(a);
		TYPE logb = log(b);
		TYPE tmp = loga - logb;
		tmp = loga < logb ? tmp + (TYPE)N : tmp;
		return !b ? (TYPE)N : !a ? (TYPE)0 : exp(tmp);

	}
	TYPE rcp(TYPE a) { return !a ? (TYPE)N : exp((TYPE)N - log(a)); }
};

#endif
