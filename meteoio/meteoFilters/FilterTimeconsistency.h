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
#ifndef __FilterTimeconsistency_H__
#define __FilterTimeconsistency_H__

#include <meteoio/meteoFilters/WindowedFilter.h> //use this one for filters relying on a data window, for example std_dev
//#include <meteoio/meteoFilters/FilterBlock.h> //use this one for all others

#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterTimeconsistency
 * @ingroup processing
 * @author Anna-Maria Tilg
 * @date   2015-12-08
 * @brief Compare sum of differences between snow depth (HS) value and HS value before and HS value afterwards, respectively, 
 *  with 4 times of the standard deviation of HS in a defined time period (something like MAD)
 * First, the standard deviation of HS is calculated for a certain time period. Afterwards, the difference between the value and 
 * the value before and the difference between the value and the following value are calculated. Then, the sum of the two 
 * differences is calculated and compared with 4 times of the standard deviation. Is the sum lower than the standard deviation, 
 * the HS value is accepted. Otherwise the HS value gets invalid. 
 * References/Literature:  Zahumensky, Igor, 2004: Guidelines on Quality Control Procedures for Data from Automatic Weather Stations, World Meteorological Organisation 
 * Remarks:
 * - nodata values are excluded from the calculation of the standard deviation ?????? => checken  
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - only window position "center" possible
 * - keyword "soft" not allowed 
	* - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position ?????? => checken
	* - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present ?????? => checken
 *
 * @code
 * Valid examples for the io.ini file:
 *          HS::filter1 = time_consistency
 *          HS::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 *          TA::filter1 = time_consistency
 *          TA::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */

//class FilterTimeconsistency : public FilterBlock { //use this one for simple filter that only look at one data point at a time, for example min_max
//class FilterTimeconsistency : public ProcessingBlock { //use this one for data corrections
class FilterTimeconsistency : public WindowedFilter { //use this one for filters relying on a data window, for example std_dev
	public:
		FilterTimeconsistency(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
};

} //end namespace

#endif
