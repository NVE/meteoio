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
#ifndef __BUFFEREDIOHANDLER_H__
#define __BUFFEREDIOHANDLER_H__

#ifdef _POPC_
#include <meteoio/IOHandler.ph>
#else
#include <meteoio/IOHandler.h>
#endif

#include <meteoio/Config.h>

#include <map>
#include <vector>
#include <string>

namespace mio {

/**
 * @class BufferedIOHandler
 * @brief This class is the class to use for buffered I/O operations. It is responsible for transparently loading the plugins
 * and transparently buffering the data. It follows the interface defined by the IOInterface class with the addition of
 * a few convenience methods.
 *
 * @author Thomas Egger
 * @date   2009-07-25
 */

class MeteoFilter;

#ifdef _POPC_
class BufferedIOHandler {
#else
class BufferedIOHandler : public IOInterface {
#endif
	public:

		/**
		 * @brief The constructor accepts an already initialized child of IOInterface (e.g. A3DIO, BormaIO, ImisIO)
		 *        and a Config object
		 *
		 * Example Usage:
		 * @code
		 * IOHandler *io1;
		 * Config cfg("io.ini");
		 * io1 = new A3DIO(cfg);
		 * BufferedIOHandler bio(*io1, cfg);
		 * @endcode
		 */
		BufferedIOHandler(IOHandler& _iohandler, const Config& _cfg);
	#ifdef _POPC_
		virtual ~BufferedIOHandler();
	#else
		virtual ~BufferedIOHandler() throw();
	#endif

		///Keywords for slope computation algorithm
		typedef enum BUFFER_POLICY {
			KEEP_NODATA, ///< when a data point is nodata in the buffer, return the buffered value
			RECHECK_NODATA ///< when a data point is nodata in the buffer, refresh the buffer to see if a value could be found
		} buffer_policy;

		/**
		 * @brief The function returns the next MeteoData object for each station with a
		 *        date >= to the parameter _date. vecMeteo and vecStation will be empty if there
		 *        is no MeteoData to be found.
		 *        NOTE: only the real measured data is looked at: no resampled values are taken into account
		 *
		 * @param _date start date of the data search for each station
		 * @param vecMeteo   A vector of MeteoData objects to be filled with data
		 */
		//void getNextMeteoData(const Date& _date, std::vector<MeteoData>& vecMeteo);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);

		/**
		 * @brief Clear all buffers in BufferedIOHandler and hence force rebuffering
		 */
		void clearBuffer();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void readDEM(DEMObject& dem_out);
		virtual void readAssimilationData(const Date& date_in, Grid2DObject& da_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const unsigned int& stationindex=IOUtils::npos);
#ifdef _POPC_
		virtual void writeMeteoData(std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");
#else
		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");
#endif
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");

		/**
		 * @brief Manually tune the buffering policy
		 * @param policy flag to define how to handle nodata (see BufferedIOHandler::buffer_policy)
		 */
		void setBufferPolicy(const buffer_policy& policy);

		/**
		 * @brief Manually tune the buffer
		 * @param _bufferbefore start date of the buffer (in offset to initially requested date)
		 * @param _bufferafter end date of the buffer (in offset to initially requested date)
		 * @code
		 * setBufferDuration(Date(2.0), Date(20.0)); //to get 2 days before the requested date to 20 days after
		 * @endcode
		 */
		void setBufferDuration(const Date& _bufferbefore, const Date& _bufferafter);

		double getAvgSamplingRate();

		friend std::ostream& operator<<(std::ostream& os, const BufferedIOHandler& data);

		friend class IOManager;

	private:
		const std::vector<METEO_DATASET>& get_complete_buffer(Date& start, Date& end);

		void setDfltBufferProperties();
		void bufferAllData(const Date& date_start, const Date& date_end);

		IOHandler& iohandler;
		Config cfg;

		bool always_rebuffer;
		Date buffer_start, buffer_end, default_chunk_size;

		std::vector< std::vector<MeteoData> > vec_buffer_meteo;
		std::map<std::string, Grid2DObject> mapBufferedGrids;
};
} //end namespace
#endif
