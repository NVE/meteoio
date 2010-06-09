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

#include "MeteoData.h"
#include "StationData.h"
#include "ConfigReader.h"
#include "libinterpol1D.h"
#include "FilterProperties.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <map>

namespace mio {

/**
 * @class FilterAlgorithms
 * @brief 
 * @author Thomas Egger
 * @date   2009-11-03
 */
class FilterAlgorithms {
	public:

		static const FilterProperties& filterProperties(const std::string& filtername);

		//Available filters
		static bool RateFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MedianAbsoluteDeviationFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool AccumulateProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool LinResamplingProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool NearestNeighbourResamplingProcess(const std::vector<MeteoData>& vecM, 
                         const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MedianAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MeanAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool WindAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
					const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);

 	private:
		static void parseFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
							const unsigned int& minArgs, const unsigned int& maxArgs,
							bool& isSoft, std::vector<double>& vecArgs_out);
		static void parseWindowFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
								const unsigned int& minArgs, const unsigned int& maxArgs,
								bool& isSoft, std::string& windowposition, std::vector<double>& vecArgs_out);

		static bool getWindowData(const std::string& filtername, const std::vector<MeteoData>& vecM,
						const unsigned int& pos,
						const Date& date, const std::vector<std::string>& _vecArgs,
						const unsigned int& paramindex, std::vector<double>& vecWindow,
						std::vector<Date> *vecDate = NULL);
		static std::map<std::string, FilterProperties> filterMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map filterMap
};
} //end namespace

#endif
