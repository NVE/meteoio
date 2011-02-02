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
			raw           = 1,
			filtered      = 1 << 1,
			resampled     = 1 << 2,
			num_of_levels = 1 << 3
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

		/**
		* @brief Fill vecMeteo with a time series of objects
		* corresponding to the interval indicated by dateStart and dateEnd.
		* Depending on the ProcessingLevel for the instance of the IOManager
		* the data returned will be either raw (read directly from the IOHandler)
		* or processed (read from an BufferedIOHandler and filtered through the 
		* MeteoProcessor
		*
		* vecMeteo will be empty if no datasets were retrieved in the interval defined
		* by dateStart and dateEnd
		*    
		* Example Usage:
		* @code
		* vector< vector<MeteoData> > vecMeteo;      //empty vector
		* Date d1(2008,06,21,11,00);       //21.6.2008 11:00
		* Date d2(2008,07,21,11,00);       //21.7.2008 11:00
		* IOManager iom(Config("io.ini"));
		* unsigned int nstations = iom.getMeteoData(d1, d2, vecMeteo);
		* @endcode
		* @param dateStart   A Date object representing the beginning of an interval (inclusive)
		* @param dateEnd     A Date object representing the end of an interval (inclusive)
		* @param vecMeteo    A vector of vector<MeteoData> objects to be filled with data
		* @return            Number of stations for which data has been found in the interval
		*/
		unsigned int getMeteoData(const Date& dateStart, const Date& dateEnd,
		                          std::vector< std::vector<MeteoData> >& vecMeteo);

		/**
		 * @brief Fill vector<MeteoData> object with multiple instances of MeteoData
		 * corresponding to the instant indicated by a Date object. Each MeteoData
		 * instance within the vector represents the data for one station at the given 
		 * instant. Depending on the ProcessingLevel configured data will be either 
		 * raw (read directly from the IOHandler)
		 *
		 * NOTE:
		 * - vecMeteo will be empty if there is no data found for any station
		 *
		 * Example Usage:
		 * @code
		 * vector<MeteoData> vecMeteo;      //empty vector
		 * IOManager iomanager(Config("io.ini"));
		 * iomanager.getMeteoData(Date(2008,06,21,11,00), vecMeteo); //21.6.2008 11:00
		 * @endcode
		 * @param i_date      A Date object representing the date/time for the sought MeteoData objects
		 * @param vecMeteo    A vector of MeteoData objects to be filled with data
		 * @return            Number of stations for which data has been found in the interval
		 */
		unsigned int getMeteoData(const Date& i_date, std::vector<MeteoData>& vecMeteo);
		
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, 
		                 Grid2DObject& result, std::string& info_string);
		
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, 
		                 Grid2DObject& result);
		
		/**
		 * @brief Set the desired ProcessingLevel of the IOManager instance
		 *        The processing level affects the way meteo data is read and processed
		 *        Three values are possible:
		 *        - IOManager::raw data shall be read directly from the buffer
		 *        - IOManager::filtered data shall be filtered before returned to the user
		 *        - IOManager::resampled data shall be resampled before returned to the user
		 *          this only affects the function getMeteoData(const Date&, METEO_DATASET&);
		 *
		 *        The three values can be combined: e.g. IOManager::filtered | IOManager:resampled
		 * @param i_level The ProcessingLevel values that shall be used to process data
		 */
		void setProcessingLevel(const unsigned int& i_level);

		double getAvgSamplingRate();

		void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string& name="");

		friend std::ostream& operator<<(std::ostream& os, const IOManager& io);

	private:
		void add_to_cache(const Date& i_date, const std::vector<MeteoData>& vecMeteo);
		void fill_filtered_cache();
		bool read_filtered_cache(const Date& start_date, const Date& end_date,
							std::vector<METEO_DATASET>& vec_meteo);

		const Config& cfg;
		IOHandler rawio;
		BufferedIOHandler bufferedio;
		MeteoProcessor meteoprocessor;
		ProcessingProperties proc_properties;

		std::map<Date, std::vector<MeteoData> > resampled_cache;  ///< stores already resampled data points
		std::vector< std::vector<MeteoData> > filtered_cache; ///< stores already filtered data intervals
		Date fcache_start, fcache_end;
		unsigned int processing_level;
};
} //end namespace
#endif
