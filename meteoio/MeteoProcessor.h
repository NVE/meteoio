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
#include <meteoio/Config.h>
#include <meteoio/Meteo1DInterpolator.h>
#include <meteoio/meteofilters/ProcessingStack.h>

#include <vector>
#include <set>

namespace mio {

/**
 * @class MeteoProcessor
 * @brief A facade class that invokes the processing of the filters and the resampling
 * @author Thomas Egger
 * @date   2010-06-25
 */

class MeteoProcessor {
	public:
		MeteoProcessor(const Config& cfg);
		~MeteoProcessor();

		/**
		 * @brief A function that executes all the filters for all meteo parameters 
		 *        configuered by the user
		 * @param[in] ivec A dataset of MeteoData
		 * @param[in] ovec The filtered output of MeteoData
		 * @param[in] second_pass Whether this is the second pass (check only filters)
		 */
		void process(const std::vector< std::vector<MeteoData> >& ivec, 
				   std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass=false);

		unsigned int resample(const Date& date, std::vector<MeteoData>& ivec);

		void getWindowSize(ProcessingProperties& o_properties);
		friend std::ostream& operator<<(std::ostream& os, const MeteoProcessor& data);

 	private:
		unsigned int get_parameters(const Config& cfg, std::set<std::string>& set_parameters);
		void compareProperties(const ProcessingProperties& newprop, ProcessingProperties& current);

		std::map<std::string, ProcessingStack*> processing_stack;

		Meteo1DInterpolator mi1d;
};
} //end namespace

#endif
