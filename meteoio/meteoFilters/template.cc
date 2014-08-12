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
#include <meteoio/meteoFilters/TEMPLATE.h>

using namespace std;

namespace mio {

TEMPLATE::TEMPLATE(const std::vector<std::string>& vec_args, const std::string& name)
          : FilterBlock(name) //this has to match the class you are inheriting from! ie FilterBlock or ProcessingBlock or WindowedFilter
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void TEMPLATE::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	for (size_t ii=0; ii<ovec.size(); ii++){
		//here, implement what has to be done on each data point
		//for example:
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		if (tmp < 0.){ //delete all values less than zero
			tmp = IOUtils::nodata;
		}
	}
}


void TEMPLATE::parse_args(std::vector<std::string> vec_args) {

	//implement the parsing of the arguments
	//you can inspire yourself with the min filter

	if ( !vec_args.empty() ) //ie if there are arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);
}

} //end namespace
