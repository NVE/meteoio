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
#ifndef __METEOFILTER_H__
#define __METEOFILTER_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>
#include <meteoio/ConfigReader.h>
#include <meteoio/libinterpol1D.h>
#include <meteoio/FilterProperties.h>
#include <meteoio/FilterAlgorithms.h>

#include <iostream>
#include <string>
#include <vector>

namespace mio {

/**
 * @class MeteoFilter
 * @brief A class that can filter (i.e. clean, resample, ...) MeteoData objects
 * @author Thomas Egger
 * @date   2009-11-01
 */

class MeteoFilter {
	public:

		/**
		* @brief 	The default constructor
		* Set up all the filters for each parameter
		* Init tasklist: a vector that holds one vector\<string\> for each parameter,
		*                representing the sequence of filters that will be executed
		*                for the respective parameter
		*                e.g. tasklist for TA: min_max, resample, min_max
		* taskargs:      a vector that holds the respective arguments for the filters 
		*                listed in tasklist
		*
		* Important: the filtering occurs in two passes:
		* Pass 1   : all filters specified are executed, resampling (if required) is the last filter applied
		* Pass 2   : all filters, that only perform checks are reapplied to the resampled values
		* @param[in] _cfg ConfigReader object that holds the MeteoFilter configuration in the [Filters] section
		*/
		MeteoFilter(const ConfigReader& _cfg);

		/**
		 * @brief A function that executes all the filters that have been setup in the constructor
		 * @param[in] vecM The raw sequence of MeteoData objects for a given station
		 * @param[in] vecS The meta data for the MeteoData objects in vecM
		 * @param[out] vecWindowM The sequence of MeteoData objects that shall be filtered
		 * @param[out] vecWindowS The meta data for the MeteoData objects in vecWindowM
		 * @param[in] checkonly indicates whether to execute all filters or only "checkonly" filters
		 */
		void filterData(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
                          std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS, 
                          const bool& checkonly = false);


 	private:
		unsigned int getFiltersForParameter(const std::string& parname, std::vector<std::string>& vecFilters);
		unsigned int getArgumentsForFilter(const std::string& keyname, std::vector<std::string>& vecArguments);		
		
		ConfigReader cfg;
		std::vector< std::vector<std::string> > tasklist;
		std::vector< std::vector< std::vector<std::string> > > taskargs;
};
} //end namespace

#endif
