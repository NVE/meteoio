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
#ifndef __METEO1DINTERPOLATOR_H__
#define __METEO1DINTERPOLATOR_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>
#include <meteoio/Config.h>
#include <meteoio/ResamplingAlgorithms.h>
#include <meteoio/ProcessingBlock.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>

namespace mio {

/**
 * @class Meteo1DInterpolator
 * @brief A class that can resample MeteoData objects
 *
 * @ingroup stats
 * @author Thomas Egger
 * @date   2010-06-24
 */

class Meteo1DInterpolator {
	public:

		/**
		* @brief 	The default constructor
		* Set up the interpolation algorithm for each parameter
		* Init tasklist: a vector that holds one std::string for each parameter,
		*                representing the interpolation algorithm that will be executed
		*                for the respective parameter
		*                e.g. tasklist for TA: linear
		* taskargs:      a vector that holds the respective arguments for the algorithms
		*                as a std::vector<std::string>, so there can be multiple arguments
		*
		* @param[in] _cfg Config object that holds the MeteoFilter configuration in the [Filters] section
		*/
		Meteo1DInterpolator(const Config& _cfg);

		/**
		 * @brief A function that executes all the resampling algorithms that have been setup in the constructor
		 * @param[in] date The requested date for a MeteoData object (to be resampled if not present)
		 * @param[in] vecM A vector of MeteoData where the new object will be inserted if not present
		 * @return    The position of the newly constructed MeteoData/StationData pair within vecM & vecS
		 */
		unsigned int resampleData(const Date& date, std::vector<MeteoData>& vecM);

		void getWindowSize(ProcessingProperties& o_properties);

		friend std::ostream& operator<<(std::ostream& os, const Meteo1DInterpolator& Interpolator);

 	private:
		std::string getInterpolationForParameter(const std::string& parname, std::vector<std::string>& vecArguments);

		Config cfg;
		std::vector<std::string> tasklist;
		std::vector< std::vector< std::string > > taskargs;
		std::map<	std::string, std::pair < std::string, std::vector < std::string > > >	extended_tasklist;
};
} //end namespace

#endif
