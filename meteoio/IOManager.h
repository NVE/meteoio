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

namespace mio {

class IOManager {

	public:
		enum ProcessingLevel {
			raw           = 0,
			filtered      = 1 << 0,
			resampled     = 1 << 1,
			num_of_levels = 1 << 2
		};

		IOManager(const Config& i_cfg);

		//Legacy support to support functionality of the IOInterface superclass:
		void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		void readDEM(DEMObject& dem_out);
		void readAssimilationData(const Date& date_in, Grid2DObject& da_out);
		void readLanduse(Grid2DObject& landuse_out);
		void readSpecialPoints(std::vector<Coords>& pts);
		void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");
		//end legacy support

		unsigned int getStationData(const Date& date, std::vector<StationData>& vecStation);

		//for an intervall of data: decide whether data should be filtered or raw
		unsigned int getMeteoData(const Date& dateStart, const Date& dateEnd,
		                          std::vector< std::vector<MeteoData> >& vecMeteo);

		//data can be raw or processed (filtered, resampled)
		unsigned int getMeteoData(const Date& i_date, std::vector<MeteoData>& vecMeteo);
		
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, 
		                 Grid2DObject& result, std::string& info_string);
		
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, 
		                 Grid2DObject& result);
		
		void setProcessingLevel(const unsigned int& i_level);

		double getAvgSamplingRate();

		void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string& name="");

		friend std::ostream& operator<<(std::ostream& os, const IOManager& io);

	private:
		void add_to_cache(const Date& i_date, const std::vector<MeteoData>& vecMeteo);

		Config cfg;
		IOHandler rawio;
		BufferedIOHandler bufferedio;
		MeteoProcessor meteoprocessor;

		std::map<Date, std::vector<MeteoData> > meteo_cache; ///< stores already fetched data points
		unsigned int processing_level;
};
} //end namespace
#endif
