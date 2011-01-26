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
#ifndef __FILTERALGORITHMS_H__
#define __FILTERALGORITHMS_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>
#include <meteoio/Config.h>
#include <meteoio/libinterpol1D.h>
#include <meteoio/FilterProperties.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <map>

namespace mio {

/**
 * @class FilterAlgorithms
 * @brief A class to filter time series of data
 *
 * @author Thomas Egger
 * @date   2009-11-03
 */
class FilterAlgorithms {
	public:

		static const FilterProperties& filterProperties(const std::string& filtername);

		//Available filters
		static void RateFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                       const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                         const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                           const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                           const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void MedianAbsoluteDeviationFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                                          const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void StandardDeviationFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                                    const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void Tukey53HFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                           const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void AccumulateProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                              const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void MedianAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                             const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void MeanAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                           const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void WindAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                           const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void ExpSmoothingProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                                const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);
		static void WMASmoothingProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
		                                const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM);

 	private:
		static bool compareMeteoData (const MeteoData& m1, const MeteoData& m2);
		static void parseFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
		                                 const unsigned int& minArgs, const unsigned int& maxArgs,
		                                 bool& isSoft, std::vector<double>& vecArgs_out);
		static void parseWindowFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
		                                       const unsigned int& minArgs, const unsigned int& maxArgs,
		                                       bool& isSoft, std::string& windowposition, std::vector<double>& vecArgs_out);
		static unsigned int getWindowData(const std::string& filtername, const std::vector<MeteoData>& vecM,
                                                  const unsigned int& pos, 
                                                  const std::vector<std::string>& _vecArgs, std::vector<MeteoData>& vecResult);
		static bool getWindowData(const std::string& filtername, const std::vector<MeteoData>& vecM,
		                          const unsigned int& pos,
		                          const Date& date, const std::vector<std::string>& _vecArgs,
		                          const unsigned int& paramindex, std::vector<double>& vecWindow,
		                          std::vector<Date> *vecDate = NULL);

		static double ExpSmoothingAlgorithm(const std::vector<MeteoData>& vecMeteo, 
		                                    const unsigned int& paramindex, const double& alpha);
		static double WMASmoothingAlgorithm(const std::vector<MeteoData>& vecMyMeteo, const unsigned int& paramindex);

		static std::map<std::string, FilterProperties> filterMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map filterMap
};
} //end namespace

#endif
