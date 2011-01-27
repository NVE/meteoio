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
#ifndef __RATEFILTER_H__
#define __RATEFILTER_H__

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  RateFilter
 * @author Thomas Egger
 * @date   2011-01-26
 * @brief Rate of change filter.
 * Calculate the change rate (ie: slope) between two points, if it is above a user given value, reject the point. Remarks:
 * - the maximum permissible rate of change (per seconds) has to be provided as an argument
 *
 * @code
 * TA::filter1	= rate
 * TA::arg1	= 0.01
 * @endcode
 */
class RateFilter : public FilterBlock {
	public:
		RateFilter(const std::vector<std::string>& vec_args);

		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
						 std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
		double max_rate_of_change;
};

} //end namespace

#endif
