#ifndef __METEOFILTER_H__
#define __METEOFILTER_H__

#include "MeteoData.h"
#include "StationData.h"
#include "ConfigReader.h"
#include "libinterpol1D.h"
#include "FilterProperties.h"
#include "FilterAlgorithms.h"
#include <iostream>
#include <string>
#include <vector>

/**
 * @class MeteoFilter
 * @brief 
 * @author Thomas Egger
 * @date   2009-11-01
 */
class MeteoFilter {
	public:

		/**
		* @brief The default constructor 
		*/
		MeteoFilter(const ConfigReader& _cfg);

		bool filterData(const std::vector<MeteoData>& vecM, 
					 const std::vector<StationData>& vecS, 
					 const unsigned int& pos, const Date_IO& date,
					 MeteoData& md, StationData& sd);


 	private:
		unsigned int getFiltersForParameter(const std::string& parname, std::vector<std::string>& vecFilters);
		unsigned int getArgumentsForFilter(const std::string& keyname, std::vector<std::string>& vecArguments);		

		ConfigReader cfg;
		std::vector< std::vector<std::string> > tasklist;
		std::vector< std::vector< std::vector<std::string> > > taskargs;
};

#endif
