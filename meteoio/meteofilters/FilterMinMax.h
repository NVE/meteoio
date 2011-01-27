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
#ifndef __FILTERMINMAX_H__
#define __FILTERMINMAX_H__

//#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) 

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterMinMax
 * @brief  
 * @author Thomas Egger
 * @date   2011-01-02
 */

class FilterMinMax : public FilterBlock {
	public:
		FilterMinMax(const std::vector<std::string>& vec_args);

		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
						 std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);

		bool is_soft;
		double min_val, max_val;
		double min_soft, max_soft;
};

} //end namespace

#endif
