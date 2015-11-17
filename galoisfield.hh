/*
reedsolomon - Reed-Solomon error correction
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef GALOISFIELD_HH
#define GALOISFIELD_HH

#include "galoisfieldtables.hh"

template <int M, int POLY, typename TYPE>
class GaloisFieldIndex;

template <int M, int POLY, typename TYPE>
class GaloisFieldValue
{
public:
	static const int Q = 1 << M, N = Q - 1;
	static_assert(M <= 8 * sizeof(TYPE), "TYPE not wide enough");
	static_assert(Q == (POLY & ~N), "POLY not of degree Q");
	TYPE v;
	GaloisFieldValue() {}
	explicit GaloisFieldValue(TYPE v) : v(v) {}
	GaloisFieldValue<M, POLY, TYPE> operator *= (GaloisFieldIndex<M, POLY, TYPE> a)
	{
		return *this = *this * a;
	}
	GaloisFieldValue<M, POLY, TYPE> operator *= (GaloisFieldValue<M, POLY, TYPE> a)
	{
		return *this = *this * a;
	}
	static const GaloisFieldValue<M, POLY, TYPE> inf()
	{
		return GaloisFieldValue<M, POLY, TYPE>(N);
	}
	static const GaloisFieldValue<M, POLY, TYPE> zero()
	{
		return GaloisFieldValue<M, POLY, TYPE>(0);
	}
};

template <int M, int POLY, typename TYPE>
class GaloisFieldIndex
{
public:
	static const int Q = 1 << M, N = Q - 1;
	static_assert(M <= 8 * sizeof(TYPE), "TYPE not wide enough");
	static_assert(Q == (POLY & ~N), "POLY not of degree Q");
	TYPE i;
	GaloisFieldIndex() {}
	explicit GaloisFieldIndex(TYPE i) : i(i) {}
	GaloisFieldIndex<M, POLY, TYPE> operator *= (GaloisFieldIndex<M, POLY, TYPE> a)
	{
		return *this = *this * a;
	}
	static const TYPE modulus()
	{
		return N;
	}
};

template <int M, int POLY, typename TYPE>
struct GaloisField
{
	static const int Q = 1 << M, N = Q - 1;
	typedef TYPE value_type;
	typedef GaloisFieldValue<M, POLY, TYPE> ValueType;
	typedef GaloisFieldIndex<M, POLY, TYPE> IndexType;
};

template <int M, int POLY, typename TYPE>
GaloisFieldIndex<M, POLY, TYPE> index(GaloisFieldValue<M, POLY, TYPE> a)
{
	return GaloisFieldIndex<M, POLY, TYPE>(GaloisFieldTables<M, POLY, TYPE>::log(a.v));
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> value(GaloisFieldIndex<M, POLY, TYPE> a) {
	return GaloisFieldValue<M, POLY, TYPE>(GaloisFieldTables<M, POLY, TYPE>::exp(a.i));
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> operator + (GaloisFieldValue<M, POLY, TYPE> a, GaloisFieldValue<M, POLY, TYPE> b)
{
	return GaloisFieldValue<M, POLY, TYPE>(a.v ^ b.v);
}

template <int M, int POLY, typename TYPE>
GaloisFieldIndex<M, POLY, TYPE> operator * (GaloisFieldIndex<M, POLY, TYPE> a, GaloisFieldIndex<M, POLY, TYPE> b)
{
	TYPE tmp = a.i + b.i;
	return GaloisFieldIndex<M, POLY, TYPE>(a.modulus() - a.i <= b.i ? tmp - a.modulus() : tmp);
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> operator * (GaloisFieldValue<M, POLY, TYPE> a, GaloisFieldValue<M, POLY, TYPE> b)
{
	return (!a.v || !b.v) ? a.zero() : value(index(a) * index(b));
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> rcp(GaloisFieldValue<M, POLY, TYPE> a)
{
	return !a.v ? a.inf() : value(GaloisFieldIndex<M, POLY, TYPE>(a.modulus() - index(a).i));
}

template <int M, int POLY, typename TYPE>
GaloisFieldIndex<M, POLY, TYPE> operator / (GaloisFieldIndex<M, POLY, TYPE> a, GaloisFieldIndex<M, POLY, TYPE> b)
{
	TYPE tmp = a.i - b.i;
	return GaloisFieldIndex<M, POLY, TYPE>(a.i < b.i ? tmp + a.modulus() : tmp);
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> operator / (GaloisFieldValue<M, POLY, TYPE> a, GaloisFieldValue<M, POLY, TYPE> b)
{
	return !b.v ? a.inf() : !a.v ? a.zero() : value(index(a) / index(b));
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> operator * (GaloisFieldIndex<M, POLY, TYPE> a, GaloisFieldValue<M, POLY, TYPE> b)
{
	return !b.v ? a.zero() : value(a * index(b));
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> operator * (GaloisFieldValue<M, POLY, TYPE> a, GaloisFieldIndex<M, POLY, TYPE> b)
{
	return !a.v ? a.zero() : value(index(a) * b);
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> fma(GaloisFieldIndex<M, POLY, TYPE> a, GaloisFieldIndex<M, POLY, TYPE> b, GaloisFieldValue<M, POLY, TYPE> c)
{
	return a * b + c;
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> fma(GaloisFieldIndex<M, POLY, TYPE> a, GaloisFieldValue<M, POLY, TYPE> b, GaloisFieldValue<M, POLY, TYPE> c)
{
	return !b.v ? c : (value(a * index(b)) + c);
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> fma(GaloisFieldValue<M, POLY, TYPE> a, GaloisFieldIndex<M, POLY, TYPE> b, GaloisFieldValue<M, POLY, TYPE> c)
{
	return !a.v ? c : (value(index(a) * b) + c);
}

template <int M, int POLY, typename TYPE>
GaloisFieldValue<M, POLY, TYPE> fma(GaloisFieldValue<M, POLY, TYPE> a, GaloisFieldValue<M, POLY, TYPE> b, GaloisFieldValue<M, POLY, TYPE> c)
{
	return (!a.v || !b.v) ? c : (value(index(a) * index(b)) + c);
}

#endif
