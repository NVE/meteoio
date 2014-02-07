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
#include "PSQLIO.h"

using namespace std;

namespace mio {
/**
 * @page psqlio PSQLIO
 * @section psql_format Format
 *
 *
 * @section psql_units Units
 *
 *
 * @section psql_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - PSQL_URL: The URL or IP of the database server
 * - PSQL_DB: The name of the database to access
 * - PSQL_USER: The username to access the server
 * - PSQL_PASS: The password to authenticate the PSQL_USER
 */

const double PSQLIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

PSQLIO::PSQLIO(const std::string& configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), endpoint(), port(), 
									   dbname(), userid(), passwd()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getParameters();
}

PSQLIO::PSQLIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), endpoint(), port(), 
								  dbname(), userid(), passwd()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getParameters();
}

PSQLIO::~PSQLIO() throw()
{

}

void PSQLIO::getParameters()
{
	port = "5432"; //The default PostgreSQL port

	cfg.getValue("PSQL_URL", "Input", endpoint);
	cfg.getValue("PSQL_PORT", "Input", port, IOUtils::nothrow);
	cfg.getValue("PSQL_DB", "Input", dbname);
	cfg.getValue("PSQL_USER", "Input", userid);
	cfg.getValue("PSQL_PASS", "Input", passwd);

	cfg.getValues("STATION", "Input", vecFixedStationID);
	cfg.getValues("MOBILE", "Input", vecMobileStationID);

	cfg.getValue("SQL_META", "Input", sql_meta);
}

void PSQLIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*name_in*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::read2DGrid(Grid2DObject& /*grid_out*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	string station_list;

	if (vecFixedStationID.empty() && vecMobileStationID.empty()) {
		return; //nothing to do
	} else {
		for (vector<string>::const_iterator it = vecFixedStationID.begin(); it != vecFixedStationID.end(); ++it) {
			if (it != vecFixedStationID.begin()) {
				station_list += ", ";
			}
			station_list += "'" + *it + "'";
		}
	}

	PGresult *result = get_data(sql_meta + " (" + station_list + ") ORDER BY id;");
	if (result) {
		int rows = PQntuples(result);
		int columns = PQnfields(result);

		int col_id = PQfnumber(result, "id");
		int col_name = PQfnumber(result, "name");
		int col_x = PQfnumber(result, "x");
		int col_y = PQfnumber(result, "y");
		int col_alt = PQfnumber(result, "altitude");
		int col_epsg = PQfnumber(result, "epsg");

		if ((col_id * col_name * col_x * col_y * col_alt * col_epsg) < 0) { //missing column
			throw IOException("Result set does not have all necessary columns", AT);
		}

		for (int ii=0; ii<rows; ii++) {
			int epsg;
			double easting, northing, altitude;

			IOUtils::convertString(epsg, PQgetvalue(result, ii, col_epsg));
			IOUtils::convertString(easting, PQgetvalue(result, ii, col_x));
			IOUtils::convertString(northing, PQgetvalue(result, ii, col_y));
			IOUtils::convertString(altitude, PQgetvalue(result, ii, col_alt));

			Coords point;
			point.setEPSG(epsg);
			point.setXY(easting, northing, altitude);

			StationData sd(point, PQgetvalue(result, ii, col_id), PQgetvalue(result, ii, col_name));
			vecStation.push_back(sd); //this is ordered ascending by id
		}

		PQclear(result);
	}
}

void PSQLIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                             const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                              const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::write2DGrid(const Grid2DObject& /*grid_in*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::cleanup() throw()
{

}

void PSQLIO::open_connection()
{
	string connect = "hostaddr = '" + endpoint +
		"' port = '" + port +
		"' dbname = '" + dbname +
		"' user = '" + userid +
		"' password = '" + passwd +
		"' connect_timeout = '10'";

	psql = PQconnectdb(connect.c_str());

	if (!psql) {
		throw IOException("PSQLIO connection error: PQconnectdb returned NULL", AT);
	}
	if (PQstatus(psql) != CONNECTION_OK) {
		throw IOException("PSQLIO connection error: PQstatus(psql) != CONNECTION_OK", AT);
	}

	//cout << "Connection established" << endl;
}

PGresult *PSQLIO::get_data(const string& sql_command)
{
	open_connection();

	PGresult *result = PQexec(psql, sql_command.c_str());
	ExecStatusType status = PQresultStatus(result);
	if (status == PGRES_TUPLES_OK) { //Successful completion of a SELECT data request
		// cout << "Select executed normally... " << endl;

		// PQprintOpt        options = {0};
		// options.header    = 1;    /* Ask for column headers            */
		// options.align     = 1;    /* Pad short columns for alignment   */
		// options.fieldSep  = "|";  /* Use a pipe as the field separator */
		// PQprint(stdout, result, &options);

	} else {
		//cout << "BAD SELECT: " << PQresStatus(status) << endl;
		PQclear(result);
		return NULL;
	}

	close_connection(psql);
	return result;
}

void PSQLIO::close_connection(PGconn *conn)
{
    PQfinish(conn);
    //cout << "Connection closed" << endl;
}

} //namespace
