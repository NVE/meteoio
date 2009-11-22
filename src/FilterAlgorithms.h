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

		static void parseFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
					const unsigned int& minArgs, const unsigned int& maxArgs, 
					bool& isSoft, std::vector<double>& vecArgs_out);

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
		static bool MedianAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
					const unsigned int& paramindex, std::vector<MeteoData>& vecFilteredM, 
					std::vector<StationData>& vecFilteredS);

 	private:
		static std::map<std::string, FilterProperties> filterMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map filterMap
};

#endif
