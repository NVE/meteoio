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
#include <set>
#include <algorithm>
#include "PSQLIO.h"

using namespace std;

namespace mio {
/**
 * @page psqlio PSQLIO
 * @section psql_format Format
 * This plugin connects to a <i>generic</i> <A HREF="www.postgresql.org/">PostgreSQL</A> server to retrieve its meteorological data. The server
 * parameters must be provided as well as the queries to retrieve the stations' data and metadata. In order to compile this plugin,
 * the development package of libpq is required (this is the PostgreSQL c client library).
 *
 * @subsection psql_meta_query Metadata query
 * This query is used to retrieve the stations' metadata. This SQL query string should retrieve the following columns as result set (in this very order):
 *
 *      id (int), name (string), x (easting as double), y (northing as double), altitude (height above sea level as double), epsg (int)
 *
 * The user is allowed to select stations with the STATIONS keyword (see below). That is why the SQL query has to end with a 'WHERE id_column_name IN' clause, for example:
 *
 *      SELECT id, station_name AS name, x_coord AS x, y_coord AS y, z AS altitude, epsg from all_stations WHERE id IN
 *
 * @subsection psql_data_query Data query
 * This query is used to retrieve the data for the user selected stations within a given time interval.
 * The SQL query may retrieve the following columns as result set (any order, only date is mandatory):
 *
 *      date (mandatory, as date), ta (double), rh (double), p (double), vw (double), dw (double), iprec (the HNW value, double), iswr (double)
 *
 * The SQL query must retrieve the data for one station only, which has to be specified as \a STATIONID (this will be dynamically replaced by the plugin).
 * To set the upper and lower bounds for the date the SQL query has to contain \a DATE_START and \a DATE_END. These keywords will be dynamically replaced by
 * the plugin with the correct date. Furthermore the resultset should be ordered by date ascending. An example for a correct SQL data query string is therefore:
 *
 *      SELECT * FROM all_measurements WHERE id = ''STATIONID'' AND date>=''DATE_START'' AND date<=''DATE_END'' ORDER BY date
 *
 * @subsection psql_exclude_file Exclude file
 * It is possible to exclude specific parameters from specific stations. This is done by listing in a CSV file for each station id, which parameters should be excluded.
 * An example of an exclude file is given below:
 * @code
 *       # Example of an exclude file, comment line
 *       ; another comment line
 *       ; the parameters to exclude are specified as comma separated values:
 *       # stationid,parameter1,parameter2
 *       1,p,RH,Iprec
 *       230,RH
 *       231,RH
 * @endcode
 *
 * @section psql_units Units
 * Units are assumed to be pure SI, except:
 *  - temperatures in °C
 *  - relative humidity in %
 *  - snow height in cm
 *  - pressure in mbar
 *
 * @section psql_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] section
 * - database connection keywords; [Input] section:
 *      - PSQL_URL: The URL or IP of the database server
 *      - PSQL_PORT: the port to use to connect
 *      - PSQL_DB: The name of the database to access
 *      - PSQL_USER: The username to access the server
 *      - PSQL_PASS: The password to authenticate the PSQL_USER
 * - database structure keywords; [Input] section
 *      - SQL_META: SQL query to use to get the stations' metadata.
 *      - SQL_DATA: SQL query to use to get the stations' data.
 * - STATIONS: comma separated list of station ids that the user is interested in; [Input] section
 * - EXCLUDE: File containing a list of parameters to exclude listed per station id (optional; [Input] section)
 *
 */

const double PSQLIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

PSQLIO::PSQLIO(const std::string& configfile) : coordin(), coordinparam(), coordout(), coordoutparam(), endpoint(), port(),
                                                dbname(), userid(), passwd(), psql(NULL), default_timezone(1.), vecMeta(),
                                                vecFixedStationID(), vecMobileStationID(), sql_meta(), sql_data(), shadowed_parameters()
{
	Config cfg(configfile);
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getParameters(cfg);
}

PSQLIO::PSQLIO(const Config& cfg) : coordin(), coordinparam(), coordout(), coordoutparam(), endpoint(), port(),
                                          dbname(), userid(), passwd(), psql(NULL), default_timezone(1.), vecMeta(),
                                          vecFixedStationID(), vecMobileStationID(), sql_meta(), sql_data(), shadowed_parameters()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getParameters(cfg);
}

PSQLIO::PSQLIO(const PSQLIO& in) : coordin(in.coordin), coordinparam(in.coordinparam), coordout(in.coordout),
                                   coordoutparam(in.coordoutparam), endpoint(in.endpoint), port(in.port), dbname(in.dbname), userid(in.userid),
                                   passwd(in.passwd), psql(NULL), default_timezone(1.), vecMeta(in.vecMeta),
                                   vecFixedStationID(in.vecFixedStationID), vecMobileStationID(in.vecMobileStationID),
                                   sql_meta(in.sql_meta), sql_data(in.sql_data), shadowed_parameters(in.shadowed_parameters) {}

PSQLIO& PSQLIO::operator=(const PSQLIO& in)
{
	PSQLIO tmp(in);

	 swap(coordin, tmp.coordin);
	 swap(coordinparam, tmp.coordinparam);
	 swap(coordout, tmp.coordout);
	 swap(coordoutparam, tmp.coordoutparam);
	 swap(endpoint, tmp.endpoint);
	 swap(port, tmp.port);
	 swap(dbname, tmp.dbname);
	 swap(userid, tmp.userid);
	 swap(passwd, tmp.passwd);
	 swap(psql, tmp.psql);
	 swap(default_timezone, tmp.default_timezone);
	 swap(vecMeta, tmp.vecMeta);
	 swap(vecFixedStationID, tmp.vecFixedStationID);
	 swap(vecMobileStationID, tmp.vecMobileStationID);
	 swap(sql_meta, tmp.sql_meta);
	 swap(sql_data, tmp.sql_data);
	 swap(shadowed_parameters, tmp.shadowed_parameters);

      return *this;
}

PSQLIO::~PSQLIO() throw() {}

void PSQLIO::getParameters(const Config& cfg)
{
	port = "5432"; //The default PostgreSQL port

	cfg.getValue("PSQL_URL", "Input", endpoint);
	cfg.getValue("PSQL_PORT", "Input", port, IOUtils::nothrow);
	cfg.getValue("PSQL_DB", "Input", dbname);
	cfg.getValue("PSQL_USER", "Input", userid);
	cfg.getValue("PSQL_PASS", "Input", passwd);

	string stations;
	cfg.getValue("STATIONS", "Input", stations);
	IOUtils::readLineToVec(stations, vecFixedStationID, ',');

	string exclude_file;
	cfg.getValue("EXCLUDE", "Input", exclude_file, IOUtils::nothrow);
	if (!exclude_file.empty() && IOUtils::fileExists(exclude_file)) {
		//if this is a relative path, prefix the path with the current path
		const std::string prefix = ( IOUtils::isAbsolutePath(exclude_file) )? "" : cfg.getConfigRootDir()+"/";
		const std::string path = IOUtils::getPath(prefix+exclude_file, true);  //clean & resolve path
		const std::string filename = path + "/" + IOUtils::getFilename(exclude_file);

		create_shadow_map(exclude_file);
	}

	cfg.getValue("SQL_META", "Input", sql_meta);
	cfg.getValue("SQL_DATA", "Input", sql_data);

	cfg.getValue("TIME_ZONE", "Input", default_timezone, IOUtils::nothrow);
}

void PSQLIO::create_shadow_map(const std::string& exclude_file)
{
	std::ifstream fin; //Input file streams
	fin.open(exclude_file.c_str(), std::ifstream::in);
	if (fin.fail()) throw FileAccessException(exclude_file, AT);

	try {
		const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		vector<string> tmpvec;
		string line;

		while (!fin.eof()) { //Go through file
			getline(fin, line, eoln); //read complete line meta information
			IOUtils::stripComments(line);
			const size_t ncols = IOUtils::readLineToVec(line, tmpvec, ',');

			if (ncols > 1) {
				for(vector<string>::iterator it = tmpvec.begin()+1; it != tmpvec.end(); ++it) {
					IOUtils::toUpper(*it);
				}

				set<string> tmpset(tmpvec.begin()+1, tmpvec.end());
				shadowed_parameters[ tmpvec[0] ] = tmpset;
			}
		}
	} catch (const std::exception&) {
		fin.close();
		throw;
	}

	fin.close();
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
	if (!vecMeta.empty()) {
		vecStation = vecMeta;
		return;
	}

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

		int col_id = PQfnumber(result, "id");
		int col_name = PQfnumber(result, "name");
		int col_x = PQfnumber(result, "x");
		int col_y = PQfnumber(result, "y");
		int col_alt = PQfnumber(result, "altitude");
		int col_epsg = PQfnumber(result, "epsg");

		if ((col_id * col_name * col_x * col_y * col_alt * col_epsg) < 0) { //missing column
			throw IOException("Result set does not have all necessary columns", AT);
		}

		vector<StationData> tmp_station;
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
			tmp_station.push_back(sd); //this is ordered ascending by id
		}

		//order according to station numbers in io.ini, PGresult is not ordered
		for (vector<string>::const_iterator it = vecFixedStationID.begin(); it != vecFixedStationID.end(); ++it) {
			station_list += "'" + *it + "'";

			for (vector<StationData>::const_iterator station_it = tmp_station.begin(); station_it != tmp_station.end(); ++station_it) {
				if ((*station_it).stationID == *it) {
					vecStation.push_back(*station_it);
				}
			}
		}

		PQclear(result);
	}
}

void PSQLIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                           std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& stationindex)
{
	if (vecMeta.empty()) readStationData(dateStart, vecMeta);
	if (vecMeta.empty()) return; //if there are no stations -> return

	size_t indexStart=0, indexEnd=vecMeta.size();

	//The following part decides whether all the stations are rebuffered or just one station
	if (stationindex == IOUtils::npos){
		vecMeteo.clear();
		vecMeteo.insert(vecMeteo.begin(), vecMeta.size(), vector<MeteoData>());
	} else {
		if (stationindex < vecMeteo.size()){
			indexStart = stationindex;
			indexEnd   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("You tried to access a stationindex in readMeteoData that is out of bounds", AT);
		}
	}

	for (size_t ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		readData(dateStart, dateEnd, vecMeteo[ii], ii);
	}
}

bool PSQLIO::replace(std::string& str, const std::string& from, const std::string& to)
{
    const size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void PSQLIO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	string sql_query(sql_data);

	string id = vecFixedStationID.at(stationindex);
	string date_start = dateStart.toString(Date::ISO);
	string date_end = dateEnd.toString(Date::ISO);
	std::replace(date_start.begin(), date_start.end(), 'T', ' ');
	std::replace(date_end.begin(), date_end.end(), 'T', ' ');

	replace(sql_query, "STATIONID", vecMeta.at(stationindex).stationID);
	replace(sql_query, "DATE_START", date_start);
	replace(sql_query, "DATE_END", date_end);

	PGresult *result = get_data(sql_query);
	if (result) {
		int rows = PQntuples(result);
		int columns = PQnfields(result);

		vector<size_t> index;
		MeteoData tmpmeteo;
		tmpmeteo.meta = vecMeta.at(stationindex);

		map_parameters(result, tmpmeteo, index);

		for (int ii=0; ii<rows; ii++) {
			parse_row(result, ii, columns, tmpmeteo, index, vecMeteo);
		}

		PQclear(result);
	}
}

void PSQLIO::parse_row(PGresult* result, const int& row, const int& cols, MeteoData& md, std::vector<size_t>& index, std::vector<mio::MeteoData>& vecMeteo)
{
	MeteoData tmp(md);
	IOUtils::convertString(md.date, PQgetvalue(result, row, 0), 0.0);

	for (int ii=1; ii<cols; ii++) {
		if (index[ii] != IOUtils::npos) {
			string val( PQgetvalue(result, row, ii) );
			if (!val.empty()) IOUtils::convertString(tmp(index[ii]), val);
		}
	}

	convertUnits(tmp);
	vecMeteo.push_back(tmp);
}

void PSQLIO::map_parameters(PGresult* result, MeteoData& md, std::vector<size_t>& index)
{
	const int columns = PQnfields(result);

	set<string> shadowed;
	map< string, set<string> >::iterator it = shadowed_parameters.find(md.meta.stationID);
	if (it != shadowed_parameters.end()) shadowed = it->second;

	for (int ii=0; ii<columns; ii++) {
		const string field_name( IOUtils::strToUpper(PQfname(result, ii)) );
		const bool is_in = shadowed.find(field_name) != shadowed.end();
		if (is_in) { // Certain parameters may be shadowed
			index.push_back(IOUtils::npos);
			continue;
		}

		if (field_name == "RH") {
			index.push_back(MeteoData::RH);
		} else if (field_name == "TA") {
			index.push_back(MeteoData::TA);
		} else if (field_name == "DW") {
			index.push_back(MeteoData::DW);
		} else if (field_name == "VW") {
			index.push_back(MeteoData::VW);
		} else if (field_name == "ISWR") {
			index.push_back(MeteoData::ISWR);
		} else if (field_name == "RSWR") {
			index.push_back(MeteoData::RSWR);
		} else if (field_name == "HS") {
			index.push_back(MeteoData::HS);
		} else if (field_name == "IPREC") {
			index.push_back(MeteoData::HNW);
		} else if (field_name == "TSS") {
			index.push_back(MeteoData::TSS);
		} else if (field_name == "TSG") {
			index.push_back(MeteoData::TSG);
		} else if (field_name == "P") {
			index.push_back(MeteoData::P);
		} else { //this is an extra parameter
			md.addParameter(field_name);
			const size_t parindex = md.getParameterIndex(field_name);
			index.push_back(parindex);
		}
	}
}

/**
* This function checks whether all the MeteoData elements in vecMeteo are consistent
* regarding their meta data (position information, station name). If they are consistent
* true is returned, otherwise false
*/
bool PSQLIO::checkConsistency(const std::vector<MeteoData>& vecMeteo, StationData& sd)
{
	if (!vecMeteo.empty()) //to get the station data even when in bug 87 conditions
		sd = vecMeteo[0].meta;

	for (size_t ii=1; ii<vecMeteo.size(); ii++){
		const Coords& p1 = vecMeteo[ii-1].meta.position;
		const Coords& p2 = vecMeteo[ii].meta.position;
		if (p1 != p2) {
			//we don't mind if p1==nodata or p2==nodata
			if(p1.isNodata()==false && p2.isNodata()==false) return false;
		}
	}

	return true;
}

void PSQLIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);

	//Loop through all stations
	/*for (size_t ii=0; ii<vecMeteo.size(); ii++){
		//1. check consistency of station data position -> write location in header or data section
		StationData sd;
		sd.position.setProj(coordout, coordoutparam);
		const bool isConsistent = checkConsistency(vecMeteo.at(ii), sd); // sd will hold valid meta info

		if (isConsistent) { //static station

		} else { //mobile station

		}
	}*/
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

void PSQLIO::convertUnits(MeteoData& meteo)
{
	//converts °C to Kelvin, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	if (ta != IOUtils::nodata)
		ta = C_TO_K(ta);

	double& tsg = meteo(MeteoData::TSG);
	if (tsg != IOUtils::nodata)
		tsg = C_TO_K(tsg);

	double& tss = meteo(MeteoData::TSS);
	if (tss != IOUtils::nodata)
		tss = C_TO_K(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
		rh /= 100.;

	double& hs = meteo(MeteoData::HS); //is in cm
	if (hs != IOUtils::nodata)
		hs /= 100.;

	double& p = meteo(MeteoData::P); //is in mbar
	if (p != IOUtils::nodata)
		p *= 100.;
}

void PSQLIO::open_connection()
{
	const string connect = "hostaddr = '" + endpoint +
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
		cerr << "ERROR" << PQstatus(psql) << endl;
		throw IOException("PSQLIO connection error: PQstatus(psql) != CONNECTION_OK", AT);
	}
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
}

} //namespace
