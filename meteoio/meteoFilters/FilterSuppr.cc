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
#include <meteoio/meteoFilters/FilterSuppr.h>

using namespace std;

namespace mio {

FilterSuppr::FilterSuppr(const std::vector<std::string>& vec_args, const std::string& name)
          : FilterBlock(name)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void FilterSuppr::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	for (size_t ii=0; ii<ovec.size(); ii++){
		ovec[ii](param) = IOUtils::nodata;
	}
}


void FilterSuppr::parse_args(std::vector<std::string> vec_args) {

	if ( !vec_args.empty() )
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);
}

} //end namespace
