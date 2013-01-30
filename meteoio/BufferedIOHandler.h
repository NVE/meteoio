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

#include <meteoio/IOHandler.h>
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
 * @section buffereiohandler_keywords Keywords
 * This module uses the following keywords to customize the buffering (all in the [General] section):
 * - BUFF_CHUNKS: how many chunks of data to buffer; (NOT YET USED)
 * - BUFF_CHUNK_SIZE: size in days of a chunk of data to read at once;
 * - BUFF_CENTERING: centering of the buffer. When rebuffering, the new date will be located BUFF_CENTERING % from the
 *                   begining of the buffer (therefore, it takes a value between 0 and 1); Optional, 10% by default.
 * - BUFF_BEFORE: alternate way of buffer centering: When rebuffering, the new date will be located BUFF_BEFORE days from the
 *                beginning of the buffer (therefore, it takes a value in days); Optional, only one of
 *                two centering option can be used.
 * - BUFF_GRIDS: how many grids to keep in the buffer. If more grids have to be read, the oldest ones will be removed from
 *               the buffer. (10 by default, 0 means no buffering for grids)
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
		BufferedIOHandler(IOHandler& in_iohandler, const Config& in_cfg);
	#ifdef _POPC_
		virtual ~BufferedIOHandler();
	#else
		virtual ~BufferedIOHandler() throw();
	#endif

		BufferedIOHandler& operator=(const BufferedIOHandler&); ///<Assignement operator

		///Keywords for slope computation algorithm
		typedef enum BUFFER_POLICY {
			KEEP_NODATA, ///< when a data point is nodata in the buffer, return the buffered value
			RECHECK_NODATA ///< when a data point is nodata in the buffer, refresh the buffer to see if a value could be found
		} buffer_policy;

		/**
		 * @brief Read the metadata for a given date.
		 * @param date date for which to read the metadata
		 * @param vecStation vector of metadata
		 */
		virtual void readStationData(const Date& date, STATION_TIMESERIE& vecStation);

		/**
		 * @brief Clear all buffers in BufferedIOHandler and hence force rebuffering
		 */
		void clearBuffer();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readAssimilationData(const Date& date_in, Grid2DObject& da_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< METEO_TIMESERIE >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);
#ifdef _POPC_
		virtual void writeMeteoData(std::vector< METEO_TIMESERIE >& vecMeteo,
		                            const std::string& name="");
#else
		virtual void writeMeteoData(const std::vector< METEO_TIMESERIE >& vecMeteo,
		                            const std::string& name="");
#endif
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

		/**
		 * @brief Returns the average sampling rate in the data.
		 * This computes the average sampling rate of the data that is contained in the buffer. This is a quick
		 * estimate, centered on how often a station measures "something" (ie, how many timestamps do we have
		 * for this station in the buffer). if the station measures TA at h+0 and h+30 and
		 * RH at h+15 and h+45, it would return 4 measurements per hour. If the station measures TA and RH at h+0 and h+30,
		 * it would return 2 measurements per hour.
		 * @return average sampling rate in Hz, nodata if the buffer is empty
		 */
		double getAvgSamplingRate();

		friend std::ostream& operator<<(std::ostream& os, const BufferedIOHandler& data);

		friend class IOManager;

		/**
		 * @brief Set buffer window properties requirements as known to the application itself.
		 * This will compare these requirements with the ones expressed by the end user and keep the max between them.
		 * The method can be called several times, it will NOT reset the calculated buffer's requirements but keep
		 * on merging with new submissions. Any parameter given as IOUtils::nodata will be ignored.
		 * @param i_chunk_size buffer size in days
		 * @param i_buff_before buffer centering in days
		 */
		void setMinBufferRequirements(const double& i_chunk_size, const double& i_buff_before);

	private:
		//private methods
		const std::vector<METEO_TIMESERIE>& get_complete_buffer(Date& start, Date& end);

		void push_meteo_data(const Date& date_start, const Date& date_end,
		                     const std::vector< METEO_TIMESERIE >& vecMeteo);

		void setDfltBufferProperties();
		void bufferData(const Date& date_start, const Date& date_end, std::vector< METEO_TIMESERIE >& vecvecMeteo);
		void bufferGrid(const Grid2DObject& in_grid2Dobj, const std::string& in_filename);

		//private members
		IOHandler& iohandler;
		const Config& cfg;

		std::vector< METEO_TIMESERIE > vec_buffer_meteo;
		std::map<std::string, Grid2DObject> mapBufferedGrids;
		std::vector<std::string> IndexBufferedGrids;

		Date buffer_start, buffer_end;
		Duration chunk_size; ///< How much data to read at once
		Duration buff_before; ///< How much data to read before the requested date in buffer
		size_t chunks; ///< How many chuncks to buffer
		size_t max_grids; ///< How many grids to buffer (grids, dems, landuse and assimilation grids together)
};

} //end namespace
#endif
