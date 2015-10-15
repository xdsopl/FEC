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
	static const int N = GF::N, K = N - NR;
	GF gf;
	value_type G[NR+1];
	ReedSolomon()
	{
		value_type root = gf.exp(FR);
		value_type pe = gf.exp(1);
		for (int i = 0; i < NR; ++i) {
			G[i] = 1;
			for (int j = i; j > 0; --j)
				G[j] = gf.fma(root, G[j], G[j-1]);
			G[0] = gf.mul(G[0], root);
			root = gf.mul(root, pe);
		}
		G[NR] = 1;
	}
	void encode(value_type *parity, const value_type *data)
	{
		for (int i = 0; i < NR; ++i)
			parity[i] = 0;
		for (int i = 0; i < K; ++i) {
			value_type feedback = gf.add(data[i], parity[0]);
			for (int j = 1; j < NR; ++j)
				parity[j-1] = gf.fma(feedback, G[NR-j], parity[j]);
			parity[NR-1] = gf.mul(G[0], feedback);
		}
	}
};

#endif
