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
#ifndef __FILTERPROPERTIES_H__
#define __FILTERPROPERTIES_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>

#include <string>
#include <vector>

namespace mio {

typedef void(*funcptr)(const std::vector<MeteoData>& vecM,
                       const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                       std::vector<MeteoData>& vecWindowM);

class FilterProperties {
	public:
		bool checkonly;
		funcptr filterfunc;
		
 		FilterProperties() : checkonly(false), filterfunc(NULL){}
 		FilterProperties(const bool& _co, const funcptr& _ptr ) : checkonly(_co), filterfunc(_ptr){}
};

} //end namespace mio

#endif
