/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __METEOCONST_H__
#define __METEOCONST_H__

namespace mio {

namespace Cst {
	const double gravity = 9.80665; // (m s-2)
	const double gas_constant_air = 287.; // (J kg-1 K-1)
	const double std_press = 101325.; // (Pa)
	const double p_water_triple_pt = 610.78; // (Pa)

	const double earth_R0 = 6356766.0; // (m)

	const double solcon = 1366.1; // (W/m^2)
} //end CST namespace

} //end namespace

#endif
