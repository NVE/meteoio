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
#ifndef __FILTERBLOCK_H__
#define __FILTERBLOCK_H__

#include <meteoio/MeteoData.h>
#include <meteoio/meteofilters/ProcessingBlock.h>
#include <vector>
#include <string>
#include <set>

namespace mio {

/**
 * @class  FilterBlock
 * @brief  An abstract class 
 * @author Thomas Egger
 * @date   2011-01-02
 */
class FilterBlock : public ProcessingBlock {
	public:
		virtual ~FilterBlock();
		
		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
						 std::vector<MeteoData>& ovec) = 0;

		static bool is_soft(std::vector<std::string>& vec_args);

		void convert_args(const unsigned int& min_nargs, const unsigned int& max_nargs,
					   const std::vector<std::string>& vec_args, std::vector<double>& dbl_args);

	protected:
		FilterBlock(const std::string& filter_name); ///< protected constructor only to be called by children
};

}
#endif
