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

#include "MeteoData.h"
#include "StationData.h"
#include <string>
#include <vector>

namespace mio {

typedef bool(*funcptr)(const std::vector<MeteoData>&, const std::vector<StationData>&, 
				   const unsigned int&, const Date&, const std::vector<std::string>&,
				   const unsigned int&, std::vector<MeteoData>&, std::vector<StationData>&);

class FilterProperties {
	public:
		bool checkonly;
		unsigned int minNbPoints;
		Date deltatime;
		funcptr filterfunc;
		
 		FilterProperties() : checkonly(false), minNbPoints(0), deltatime(0.0), filterfunc(NULL){}
 		FilterProperties(const bool& _co, const unsigned int& _points, const Date& _date, const funcptr& _ptr ) 
			: checkonly(_co), minNbPoints(_points), deltatime(_date), filterfunc(_ptr){}
};

} //end namespace mio

#endif
