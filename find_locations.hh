/*
FEC - Forward error correction
Written in 2017 by <Ahmet Inan> <xdsopl@gmail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#ifndef FIND_LOCATIONS_HH
#define FIND_LOCATIONS_HH

#include "galois_field.hh"
#include "chien.hh"

template <int NR, typename GF>
struct FindLocations
{
	typedef typename GF::value_type value_type;
	typedef typename GF::ValueType ValueType;
	typedef typename GF::IndexType IndexType;
	static int search(ValueType *locator, int locator_degree, IndexType *locations)
	{
		if (locator_degree == 1) {
			locations[0] = (index(locator[0]) / index(locator[1])) / IndexType(1);
			return 1;
		}
		if (locator_degree == 2) {
			if (!locator[1] || !locator[0])
				return 0;
			ValueType a(locator[2]), b(locator[1]), c(locator[0]);
			ValueType ba(b/a), R(Artin_Schreier_imap(a*c/(b*b)));
			if (!R)
				return 0;
			locations[0] = index(ba * R) / IndexType(1);
			locations[1] = index(ba * R + ba) / IndexType(1);
			return 2;
		}
		return Chien<NR, GF>::search(locator, locator_degree, locations);
	}
};

#endif
