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
#ifndef __PROCADD_H__
#define __PROCADD_H__

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcAdd
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2012-02-06
 * @brief Add an offset to the values.
 * This simply adds to all values a given constant.
 * @code
 * TSS::filter1	= add
 * TSS::arg1	= 0.5
 * @endcode
 */

class ProcAdd : public ProcessingBlock {
	public:
		ProcAdd(const std::vector<std::string>& vec_args);

		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector<std::string>& vec_args);

		double offset;
};

} //end namespace

#endif
