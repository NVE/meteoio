/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __FilterFlatline_H__
#define __FilterFlatline_H__

#include <meteoio/meteoFilters/WindowedFilter.h> //use this one for filters relying on a data window, for example std_dev

#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterFlatline
 * @ingroup processing
 * @brief Searching for periods of constant values at set them invalid
 * The filter searches for time periods in which the value of the certain variable stagnates / doesn't change 
 * This is done with calculating the variance. 
 * References/Literature:  
 * Remarks:
 * - nodata values are excluded from the mean ?????? => checken  
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *     - minimal number of points in window
 *     - minimal time interval spanning the window (in seconds)
 * - only window position "center" possible
 * - keyword "soft" not allowed 
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position ?????? => checken
 * - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present ?????? => checken
 *
 * @code
 * Valid examples for the io.ini file:
 *          HS::filter1 = flat_line
 *          HS::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 *          TA::filter1 = flat_line
 *          TA::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 * 
 * @author Anna-Maria Tilg
 * @date   2015-12-04
 */

class FilterFlatline : public WindowedFilter {
	public:
		FilterFlatline(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
};

} //end namespace

#endif
