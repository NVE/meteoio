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

		static void parseFilterArguments(const string& filtername, const vector<string>& vecArgs_in,
								   const unsigned int& minArgs, const unsigned int& maxArgs, 
								   bool& isSoft, vector<double>& vecArgs_out);
		static bool RateFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
						   const unsigned int& paramindex,
						   vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS);
		static bool ResamplingFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
							    const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
							    const unsigned int& paramindex,
							    vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS);
		static bool MinMaxFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
							const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
							const unsigned int& paramindex,
							vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS);
		static bool MinValueFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
							  const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
							  const unsigned int& paramindex,
							  vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS);
		static bool MaxValueFilter(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
							  const unsigned int& pos, const Date_IO& date, const vector<string>& _vecArgs,
							  const unsigned int& paramindex,
							  vector<MeteoData>& vecFilteredM, vector<StationData>& vecFilteredS);

 	private:
		static std::map<std::string, FilterProperties> filterMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map filterMap
};

#endif
