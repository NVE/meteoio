#ifndef __FILTERPROPERTIES_H__
#define __FILTERPROPERTIES_H__

#include "MeteoData.h"
#include "StationData.h"
#include <string>
#include <vector>

typedef bool(*funcptr)(const std::vector<MeteoData>&, const std::vector<StationData>&, 
				   const unsigned int&, const Date_IO&, const std::vector<std::string>&,
				   const unsigned int&, std::vector<MeteoData>&, std::vector<StationData>&);

class FilterProperties {
	public:
		bool checkonly;
		unsigned int minNbPoints;
		Date_IO deltatime;
		funcptr filterfunc;
		
 		FilterProperties() : checkonly(false), minNbPoints(0), deltatime(0.0), filterfunc(NULL){}
 		FilterProperties(const bool& _co, const unsigned int& _points, const Date_IO& _date, const funcptr& _ptr ) 
			: checkonly(_co), minNbPoints(_points), deltatime(_date), filterfunc(_ptr){}
};

#endif
