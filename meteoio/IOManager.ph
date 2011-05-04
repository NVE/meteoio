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
#ifndef __IOMANAGER_H__
#define __IOMANAGER_H__

#include <meteoio/BufferedIOHandler.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/MeteoProcessor.h>
#include <meteoio/marshal_meteoio.h>
#include <meteoio/IOPlugin.h>

using namespace mio; //HACK for POPC

namespace mio {

typedef std::map<std::string, IOPlugin>::iterator PLUGIN_ITERATOR; //HACK for POPC
//This one line above does absolutely nothing, but if removed, popc does not compile the file....


/**
* @file IOManager.ph
* The is the parclass implementing the user accessible interface.
*/

parclass IOManager {
		classuid(1004);
	public:
		IOManager([in]const Config& i_cfg);

		//Legacy support to support functionality of the IOInterface superclass:
		void read2DGrid([out]Grid2DObject& grid_out, [in]const std::string& parameter="");
		void readDEM([out]DEMObject& dem_out);
		void readAssimilationData([in]const Date& date_in, [out]Grid2DObject& da_out);
		void readLanduse([out]Grid2DObject& landuse_out);
		void readSpecialPoints([out]std::vector<Coords>& pts);
		void write2DGrid([in]const Grid2DObject& grid_in, [in]const std::string& options="");
		//end legacy support

		unsigned int getStationData([in]const Date& date, [out, proc=marshal_STATION_TIMESERIE]STATION_TIMESERIE& vecStation);

		unsigned int getMeteoData([in]const Date& dateStart, [in]const Date& dateEnd,
		                          [out, proc=marshal_vector_METEO_TIMESERIE]std::vector< METEO_TIMESERIE >& vecMeteo);

		unsigned int getMeteoData([in]const Date& i_date, [out, proc=marshal_METEO_TIMESERIE]METEO_TIMESERIE& vecMeteo);

		void push_meteo_data([in] const unsigned int& level, [in] const Date& date_start, [in] const Date& date_end,
		                     [in] const std::vector< METEO_TIMESERIE >& vecMeteo);

		void interpolate([in] Date& date, [in] DEMObject& dem, [in, proc=marshal_MeteoParameters] MeteoData::Parameters& meteoparam,
		                 [out]Grid2DObject& result, [out]std::string& info_string);

		void interpolate([in] Date& date, [in] DEMObject& dem, [in, proc=marshal_MeteoParameters] MeteoData::Parameters& meteoparam,
		                 [out]Grid2DObject& result); //HACK popc

		void setProcessingLevel([in]const unsigned int& i_level);

		double getAvgSamplingRate();

		void writeMeteoData([in ,proc=marshal_vector_METEO_TIMESERIE] std::vector< METEO_TIMESERIE >& vecMeteo, [in] std::string name=""); //HACK popc

		std::string toString() /*const*/; //HACK popc
		//friend std::ostream& operator<<(std::ostream& os, const IOManager& io);

	private:
		void add_to_cache(const Date& i_date, const METEO_TIMESERIE& vecMeteo);
		void fill_filtered_cache();
		bool read_filtered_cache(const Date& start_date, const Date& end_date,
		                         std::vector< METEO_TIMESERIE >& vec_meteo);

		//const Config& cfg; //HACK popc
		Config cfg; //HACK popc
		IOHandler rawio;
		BufferedIOHandler bufferedio;
		MeteoProcessor meteoprocessor;
		ProcessingProperties proc_properties;

		std::map<Date, METEO_TIMESERIE > resampled_cache;  ///< stores already resampled data points
		std::vector< METEO_TIMESERIE > filtered_cache; ///< stores already filtered data intervals
		Date fcache_start, fcache_end; ///< store the beginning and the end date of the filtered_cache
		unsigned int processing_level;

		unsigned int raw, filtered, resampled, num_of_levels; //HACK POPC
};

} //end namespace
#endif
