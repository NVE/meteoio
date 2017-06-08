/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PPHASEGENERATOR_H
#define PPHASEGENERATOR_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

namespace mio {

/**
 * @class PPhaseGenerator
 * @ingroup parametrizations
 * @brief Generate precipitation splitting according to the selected method
 * @details
 * It takes the following arguments:
 *  - TYPE: the splitting method to use, any of the following:
 *     - THRESH: a provided fixed air temperature threshold splits precipitation as either fully solid or fully liquid
 *     - RANGE: two air temperature thresholds provide the lower and upper range for fully solid / fully liquid precipitation.
 *                 Within the provided range, a linear transition is assumed.
 *  - SNOW: when using a fixed air temperature threshold, this gives the snow/rain threshold (in K). When using two air temperatures
 * thresholds, this provides the temperature below which only solid precipitation is found (in K);
 *  - RAIN: when using two air temperatures thresholds, this provides the temperature above which only liquid precipitation is found (in K);
 *
 * @code
 * PSUM_PH::generators     = PPHASE
 * PSUM_PH::PPHASE::type   = THRESH
 * PSUM_PH::PPHASE::snow   = 274.35
 * @endcode
 */
class PPhaseGenerator : public GeneratorAlgorithm {
	public:
		PPhaseGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo)
			: GeneratorAlgorithm(vecArgs, i_algo), model(THRESH), fixed_thresh(IOUtils::nodata),
			range_start(IOUtils::nodata), range_norm(IOUtils::nodata) { parse_args(vecArgs); }

		bool generate(const size_t& param, MeteoData& md);
		bool create(const size_t& param, std::vector<MeteoData>& vecMeteo);

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);

		typedef enum PARAMETRIZATION {
			THRESH,
			RANGE
		} parametrization;
		parametrization model;
		double fixed_thresh, range_start, range_norm;
};

} //end namespace mio

#endif
