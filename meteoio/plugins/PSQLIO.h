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
#ifndef __PSQLIO_H__
#define __PSQLIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>

#include <libpq-fe.h>
#include <string>
#include <algorithm>

namespace mio {

/**
 * @class PSQLIO
 * @brief 
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2014-01-28
 */
class PSQLIO : public IOInterface {
	public:
		PSQLIO(const std::string& configfile);
		PSQLIO(const PSQLIO&);
		PSQLIO(const Config& cfgreader);
		~PSQLIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readPOI(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void getParameters();
		void open_connection();
		PGresult* get_data(const std::string& sqlcommand);
		bool replace(std::string& str, const std::string& from, const std::string& to);
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		void map_parameters(PGresult* result, MeteoData& md, std::vector<size_t>& index);
		void parse_row(PGresult* result, const int& row, const int& cols, 
		               MeteoData& md, std::vector<size_t>& index, std::vector<mio::MeteoData>& vecMeteo);
		void close_connection(PGconn *conn);
		void convertUnits(MeteoData& meteo) const;

		const Config cfg;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::string endpoint, port, dbname, userid, passwd; ///< Variables for endpoint configuration
		PGconn *psql; ///<holds the current connection
		double default_timezone;
		std::vector<StationData> vecMeta;
		std::map<size_t, double> multiplier, offset;
		std::vector<std::string> vecFixedStationID, vecMobileStationID;
		std::string sql_meta, sql_data;
};

} //namespace
#endif
