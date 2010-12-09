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
#ifndef __RESAMPLINGALGORITHMS_H__
#define __RESAMPLINGALGORITHMS_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>
#include <meteoio/libinterpol1D.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <map>

namespace mio {

typedef void(*resamplingptr)(const unsigned int& position, const unsigned int& paramindex,
					    const std::vector<std::string>& taskargs, std::vector<MeteoData>& vecM);

/**
 * @class ResamplingAlgorithms
 * @brief 
 * @author Thomas Egger
 * @date   2010-06-25
 */
class ResamplingAlgorithms {
	public:

		static const resamplingptr& getAlgorithm(const std::string& algorithmname);

		//Available algorithms
		static void LinearResampling(const unsigned int& position, const unsigned int& paramindex,
		                             const std::vector<std::string>& taskargs, std::vector<MeteoData>& vecM);
		static void NearestNeighbour(const unsigned int& position, const unsigned int& paramindex, 
		                             const std::vector<std::string>& taskargs, std::vector<MeteoData>& vecM);

		static void Accumulate(const unsigned int& pos, const unsigned int& paramindex,
		                       const std::vector<std::string>& taskargs, std::vector<MeteoData>& vecM);
		
 	private:
		static double funcval(const std::vector<MeteoData>& vecM, const unsigned int& index, 
						  const Date& date, const unsigned int& paramindex);

		static std::map<std::string, resamplingptr> algorithmMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map algorithmMap
};
} //end namespace

#endif
