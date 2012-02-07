/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __PROCUNDERCATCH_H__
#define __PROCUNDERCATCH_H__

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcUndercatch
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2012-02-06
 * @brief Correct precipitation for undercatch in winter conditions.
 * @details
 * This implements the standard methods for precipitation correction as described in
 * <i>"WMO Solid Precipitation Measurement Intercomparison"</i>, B. Goodison, P. Louie and D. Yang, <b>872</b>, 1998 as well as
 * the overview given by <i>"Literature Study on the Correction of Precipitation Measurements"</i>, Annette Wagner, 2009.
 * These methods process pure snow and mixed precipitation differently, with the following thresholds:
 * - pure snow below -2 C
 * - mixed precipitation between -2 and +2 C
 * - pure rain above 2 C
 *
 * They also depend on the usage of a shield around the gauge as well as the type of rain gauge that does the measurements,
 * therefore this type must be specified as an argument. The coefficients are not always available both for shielded and
 * unshielded gauges, so most of the time only one variation will be available and is specified below.
 * The following methods can be specified as argument (only one can be specified):
 * - cst {factor for snow} {factor for mixed precipitation} - this applies a constant factor to the precipitation
 * - Nipher - National standard rain gauge in Canada, shielded
 * - Tretyakov - Designed in USSR in the 1950s, deployed by some national networks in ex-USSR territories, shielded
 * - US8sh - US national 8\" rain gauge, shielded
 * - US8unsh - US national 8\" rain gauge, unshielded
 * - Hellmann - the most widely used rain gauge in the world, with some country specific variations, unshielded
 *
 * @code
 * HNW::filter1	= undercatch
 * HNW::arg1	= cst 1.3 1.1
 * @endcode
 */

class ProcUndercatch : public ProcessingBlock {
	public:
		ProcUndercatch(const std::vector<std::string>& vec_args);

		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum SENSOR_TYPE {
			cst,
			nipher,
			tretyakov,
			us8sh,
			us8unsh,
			hellmann
		} sensor_type;

		typedef enum PRECIP_TYPE {
			rain,
			mixed,
			snow
		} precip_type;

		void parse_args(std::vector<std::string> filter_args);

		sensor_type type;
		double factor_snow, factor_mixed;
		static const double Tsnow, Train;
};

} //end namespace

#endif
