/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef GALOIS_FIELD_HH
#define GALOIS_FIELD_HH

namespace GF {

#include "galois_field_tables.hh"

template <int M, int POLY, typename TYPE>
struct Index;

template <int M, int POLY, typename TYPE>
struct Value
{
	static const int Q = 1 << M, N = Q - 1;
	static_assert(M <= 8 * sizeof(TYPE), "TYPE not wide enough");
	static_assert(Q == (POLY & ~N), "POLY not of degree Q");
	TYPE v;
	Value() {}
	explicit Value(TYPE v) : v(v)
	{
		assert(v <= N);
	}
	explicit operator bool () const { return v; }
	explicit operator int () const { return v; }
	Value<M, POLY, TYPE> operator *= (Index<M, POLY, TYPE> a)
	{
		assert(a.i < a.modulus());
		return *this = *this * a;
	}
	Value<M, POLY, TYPE> operator *= (Value<M, POLY, TYPE> a)
	{
		assert(a.v <= a.N);
		return *this = *this * a;
	}
	Value<M, POLY, TYPE> operator += (Value<M, POLY, TYPE> a)
	{
		assert(a.v <= a.N);
		return *this = *this + a;
	}
	static const Value<M, POLY, TYPE> zero()
	{
		return Value<M, POLY, TYPE>(0);
	}
};

template <int M, int POLY, typename TYPE>
struct Index
{
	static const int Q = 1 << M, N = Q - 1;
	static_assert(M <= 8 * sizeof(TYPE), "TYPE not wide enough");
	static_assert(Q == (POLY & ~N), "POLY not of degree Q");
	TYPE i;
	Index() {}
	explicit Index(TYPE i) : i(i)
	{
		assert(i < modulus());
	}
	explicit operator int () const { return i; }
	Index<M, POLY, TYPE> operator *= (Index<M, POLY, TYPE> a)
	{
		assert(a.i < a.modulus());
		assert(i < modulus());
		return *this = *this * a;
	}
	Index<M, POLY, TYPE> operator /= (Index<M, POLY, TYPE> a)
	{
		assert(a.i < a.modulus());
		assert(i < modulus());
		return *this = *this / a;
	}
	static const TYPE modulus()
	{
		return N;
	}
};

template <int WIDTH, int POLY, typename TYPE>
struct Types
{
	static const int M = WIDTH, Q = 1 << M, N = Q - 1;
	typedef TYPE value_type;
	typedef Value<M, POLY, TYPE> ValueType;
	typedef Index<M, POLY, TYPE> IndexType;
};

template <int M, int POLY, typename TYPE>
Index<M, POLY, TYPE> index(Value<M, POLY, TYPE> a)
{
	assert(a.v <= a.N);
	assert(a.v);
	return Index<M, POLY, TYPE>(Tables<M, POLY, TYPE>::log(a.v));
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> value(Index<M, POLY, TYPE> a) {
	assert(a.i < a.modulus());
	return Value<M, POLY, TYPE>(Tables<M, POLY, TYPE>::exp(a.i));
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> Artin_Schreier_imap(Value<M, POLY, TYPE> a) {
	assert(a.v <= a.N);
	assert(a.v);
	return Value<M, POLY, TYPE>(Tables<M, POLY, TYPE>::Artin_Schreier_imap(a.v));
}

template <int M, int POLY, typename TYPE>
bool operator == (Value<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.v <= b.N);
	return a.v == b.v;
}

template <int M, int POLY, typename TYPE>
bool operator != (Value<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.v <= b.N);
	return a.v != b.v;
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator + (Value<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.v <= b.N);
	return Value<M, POLY, TYPE>(a.v ^ b.v);
}

template <int M, int POLY, typename TYPE>
Index<M, POLY, TYPE> operator * (Index<M, POLY, TYPE> a, Index<M, POLY, TYPE> b)
{
	assert(a.i < a.modulus());
	assert(b.i < b.modulus());
	TYPE tmp = a.i + b.i;
	return Index<M, POLY, TYPE>(a.modulus() - a.i <= b.i ? tmp - a.modulus() : tmp);
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator * (Value<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.v <= b.N);
	return (!a.v || !b.v) ? a.zero() : value(index(a) * index(b));
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> rcp(Value<M, POLY, TYPE> a)
{
	assert(a.v <= a.N);
	assert(a.v);
	return value(Index<M, POLY, TYPE>(a.modulus() - index(a).i));
}

template <int M, int POLY, typename TYPE>
Index<M, POLY, TYPE> operator / (Index<M, POLY, TYPE> a, Index<M, POLY, TYPE> b)
{
	assert(a.i < a.modulus());
	assert(b.i < b.modulus());
	TYPE tmp = a.i - b.i;
	return Index<M, POLY, TYPE>(a.i < b.i ? tmp + a.modulus() : tmp);
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator / (Value<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.v <= b.N);
	assert(b.v);
	return !a.v ? a.zero() : value(index(a) / index(b));
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator / (Index<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.i < a.modulus());
	assert(b.v <= b.N);
	assert(b.v);
	return value(a / index(b));
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator / (Value<M, POLY, TYPE> a, Index<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.i < b.modulus());
	return !a.v ? a.zero() : value(index(a) / b);
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator * (Index<M, POLY, TYPE> a, Value<M, POLY, TYPE> b)
{
	assert(a.i < a.modulus());
	assert(b.v <= b.N);
	return !b.v ? b.zero() : value(a * index(b));
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> operator * (Value<M, POLY, TYPE> a, Index<M, POLY, TYPE> b)
{
	assert(a.v <= a.N);
	assert(b.i < b.modulus());
	return !a.v ? a.zero() : value(index(a) * b);
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> fma(Index<M, POLY, TYPE> a, Index<M, POLY, TYPE> b, Value<M, POLY, TYPE> c)
{
	assert(a.i < a.modulus());
	assert(b.i < b.modulus());
	assert(c.v <= c.N);
	return value(a * b) + c;
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> fma(Index<M, POLY, TYPE> a, Value<M, POLY, TYPE> b, Value<M, POLY, TYPE> c)
{
	assert(a.i < a.modulus());
	assert(b.v <= b.N);
	assert(c.v <= c.N);
	return !b.v ? c : (value(a * index(b)) + c);
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> fma(Value<M, POLY, TYPE> a, Index<M, POLY, TYPE> b, Value<M, POLY, TYPE> c)
{
	assert(a.v <= a.N);
	assert(b.i < b.modulus());
	assert(c.v <= c.N);
	return !a.v ? c : (value(index(a) * b) + c);
}

template <int M, int POLY, typename TYPE>
Value<M, POLY, TYPE> fma(Value<M, POLY, TYPE> a, Value<M, POLY, TYPE> b, Value<M, POLY, TYPE> c)
{
	assert(a.v <= a.N);
	assert(b.v <= b.N);
	assert(c.v <= c.N);
	return (!a.v || !b.v) ? c : (value(index(a) * index(b)) + c);
}

}
#endif
