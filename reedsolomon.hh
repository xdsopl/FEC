/*
reedsolomon - Reed-Solomon error correction
Written in 2015 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef REEDSOLOMON_HH
#define REEDSOLOMON_HH

#include "galoisfield.hh"

template <int NR, int FR, typename GF>
class ReedSolomon
{
public:
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static const int N = GF::N, K = N - NR;
	ValueType G[NR+1];
	ReedSolomon()
	{
		IndexType root(FR), pe(1);
		for (int i = 0; i < NR; ++i) {
			G[i] = ValueType(1);
			for (int j = i; j > 0; --j)
				G[j] = fma(root, G[j], G[j-1]);
			G[0] *= root;
			root *= pe;
		}
		G[NR] = ValueType(1);
	}
	void encode(ValueType *parity, const ValueType *data)
	{
		for (int i = 0; i < NR; ++i)
			parity[i] = ValueType(0);
		for (int i = 0; i < K; ++i) {
			ValueType feedback = data[i] + parity[0];
			for (int j = 1; j < NR; ++j)
				parity[j-1] = fma(feedback, G[NR-j], parity[j]);
			parity[NR-1] = G[0] * feedback;
		}
	}
	void encode(value_type *parity, const value_type *data)
	{
		encode(reinterpret_cast<ValueType *>(parity), reinterpret_cast<const ValueType *>(data));
	}
};

#endif
