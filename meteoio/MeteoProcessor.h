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
#ifndef __METEOPROCESSOR_H__
#define __METEOPROCESSOR_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>
#include <meteoio/ConfigReader.h>
#include <meteoio/MeteoFilter.h>
#include <meteoio/Meteo1DInterpolator.h>

#include <vector>

namespace mio {

/**
 * @class MeteoProcessor
 * @brief A facade class that invokes the processing of the filters and the resampling
 * @author Thomas Egger
 * @date   2010-06-25
 */

class MeteoProcessor {
	public:
		MeteoProcessor(const ConfigReader& _cfg);

		/**
		 * @brief A function that executes all the filters that have been setup in the constructor
		 * @param[in] date The requested date for a MeteoData object
		 * @param[in] vecM The raw sequence of MeteoData objects for a given station
		 * @param[in] vecS The meta data for the MeteoData objects in vecM
		 * @param[out] md The MeteoData object to be returned
		 * @param[out] sd The associated StationData object for md
		 */
		void processData(const Date& date, const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
		                 MeteoData& md, StationData& sd);

 	private:
		ConfigReader cfg;
		MeteoFilter mf;
		Meteo1DInterpolator mi1d;
};
} //end namespace

#endif
