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

		/**
		 * @brief Fill vector<MeteoData> object with multiple datasets
		 * corresponding to the time indicated by the Date object.
		 * Matching rule: Find first data set for every station which has an event time (measurement time)
		 * that is greater (newer) or equal to the time represented by the Date object parameter. The
		 * vector<StationData> object holds multiple StationData objects representing meta information
		 * about the meteo stations that recorded the meteo data.
		 *
		 * NOTE:
		 * - vecMeteo will contain nodata objects if an exact time match is impossible
		 *   and resampling is turned off. If resampling is turned on a resampled value is returned
		 *   if resampling is possible (enough measurements), otherwise nodata objects will be returned
		 * - is there is absolutely no data to be found, and hence not even station data, vecMeteo and vecStation
		 *   will be filled with only one nodata obejct of MeteoData and StationData respectively
		 *
		 * Example Usage:
		 * @code
		 * vector<MeteoData> vecMeteo;      //empty vector
		 * IOManager iomanager(Config("io.ini"));
		 * iomanager.getMeteoData(Date(2008,06,21,11,00), vecMeteo); //21.6.2008 11:00
		 * @endcode
		 * @param i_date      A Date object representing the date/time for the sought MeteoData objects
		 * @param vecMeteo    A vector of MeteoData objects to be filled with data
		 */
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
