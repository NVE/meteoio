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
#ifndef __FILTERPASSIVET_H__
#define __FILTERPASSIVET_H__

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterPassiveT
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2012-04-23
 * @brief Filters and correct temperatures from unventillated sensor.
 * This filter can ONLY be applied to air temperatures.
 *
 * @code
 * TA::filter2	= passive_T
 * TA::arg2	= nakamura
 * @endcode
 */

class FilterPassiveT : public FilterBlock {
	public:
		FilterPassiveT(const std::vector<std::string>& vec_args);

		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum TA_CORRECTION {
			none,
			nakamura
		} ta_correction;
		
		void parse_args(std::vector<std::string> vec_args);
		
		ta_correction type;
};

} //end namespace

#endif
