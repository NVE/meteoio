#ifndef FILTERBASEMSTN_H_INCLUDED
#define FILTERBASEMSTN_H_INCLUDED

#include "FilterBase.h"

using namespace std;

/**
 * @class FilterBase1Stn
 * @brief Base class (interface) for data filters working on a multiple stations.
 * @author Florian Hof
 * @date   2009-03-13
 */
class FilterBaseMStn : public FilterBase {

	public:

		// check handling

		/**
		* Check the meteo data at a single time with the filters.
		* Before calling doCheck, the filter's parameters have to be filled and the prepareCheck method called. 
		* @param stations   [in] The list of stations to which the meteo data corresponds.
		* @param meteoBuffers   [in] For each station, the buffer (over the past) of meteo data, used by filters that needs a time frame. 
		* @param meteoDatas   [in out] For each station, the meteo data to filter. 
		*/
		virtual void doCheck(vector<StationData>& stations, vector<MeteoBuffer>& filteredMeteoBuffers, vector<MeteoData>& meteoDatas) = 0;

};

/**
 * @brief Type for functions that create concrete multi-stations filters (constructor-like).
 */
typedef FilterBaseMStn*(*FilterMStnCreator)(); 

#endif 
