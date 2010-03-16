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
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool ResamplingFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		static bool WindAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		/**
		 * @details The median average filter returns the median value of all values within a certain window.
		 *          The size of the window is defined by an argument describing the minimal number of points 
		 *          and one describing the time span (minutes) of the window. Whether the window is centered, 
		 *          left shifted or right shifted may be configured with the keywords "left", "right" and "center" 
		 *          ("center" is the default). Finally the keyword "soft" indicated whether the filter shall strictly 
		 *          enforce the user settings with regards to the centricity of the window. If the filter is 
		 *          configured as "soft" and the centricity of the window is right then in the case, that there
		 *          is not enough data on the right, data will be added to the window from the left.
		 *
		 * @code
		 * Valid examples for the io.ini file:
		 *          TA::filter1 = median_avg
		 *          TA::arg1    = soft left 1 300 (300 minutes time span for the left leaning window)
		 *          RH::filter1 = median_avg
		 *          RH::arg1    = 10 100          (strictly centered window spanning 100 minutes and at least 10 points)
		 * @endcode
		 */
		static bool MedianAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);
		/**
		 * @details The mean average filter returns the mean value of all values within a certain window.
		 *          The size of the window is defined by an argument describing the minimal number of points 
		 *          and one describing the time span (minutes) of the window. Whether the window is centered, 
		 *          left shifted or right shifted may be configured with the keywords "left", "right" and "center" 
		 *          ("center" is the default). Finally the keyword "soft" indicated whether the filter shall strictly 
		 *          enforce the user settings with regards to the centricity of the window. If the filter is 
		 *          configured as "soft" and the centricity of the window is right then in the case, that there
		 *          is not enough data on the right, data will be added to the window from the left.
		 *
		 * @code
		 * Valid examples for the io.ini file:
		 *          TA::filter1 = mean_avg
		 *          TA::arg1    = soft left 1 300 (300 minutes time span for the left leaning window)
		 *          RH::filter1 = mean_avg
		 *          RH::arg1    = 10 100          (strictly centered window spanning 100 minutes and at least 10 points)
		 * @endcode
		 */
		static bool MeanAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
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
						const Date_IO& date, const std::vector<std::string>& _vecArgs,
						const unsigned int& paramindex, std::vector<double>& vecWindow,
						std::vector<Date_IO> *vecDate = NULL);
		static std::map<std::string, FilterProperties> filterMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map filterMap
};

#endif
