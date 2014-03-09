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
									   dbname(), userid(), passwd(), psql(NULL), default_timezone(1.), vecMeta(), multiplier(), offset(),
                                                vecFixedStationID(), vecMobileStationID(), sql_meta(), sql_data()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getParameters();
}

PSQLIO::PSQLIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), endpoint(), port(), 
								  dbname(), userid(), passwd(), psql(NULL), default_timezone(1.), vecMeta(), multiplier(), offset(),
                                          vecFixedStationID(), vecMobileStationID(), sql_meta(), sql_data()
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

	string stations("");
	cfg.getValue("STATIONS", "Input", stations);
	IOUtils::readLineToVec(stations, vecFixedStationID, ',');

	string exclude_file("");
	cfg.getValue("EXCLUDE", "Input", exclude_file, IOUtils::nothrow);
	if (IOUtils::fileExists(exclude_file)) {
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
		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
		
		vector<string> tmpvec;
		string line("");
		
		while (!fin.eof()) { //Go through file
			getline(fin, line, eoln); //read complete line meta information
			IOUtils::stripComments(line);
			const size_t ncols = IOUtils::readLineToVec(line, tmpvec, ',');

			if (ncols > 1) {
				set<string> tmpset(tmpvec.begin()+1, tmpvec.end());
				shadowed_parameters[tmpvec[0]] = tmpset;
			}
		}
	} catch (const std::exception&) {
		fin.close();
		throw;
	}
	/*
	map< string, set<string> >::iterator it;
	set<string>::iterator setit;
	for (it = shadowed_parameters.begin(); it != shadowed_parameters.end(); ++it) {
		cout << "Shadowed for station " << it->first << "   ";
		for (setit = it->second.begin(); setit != it->second.end(); ++setit) {
			cout << *setit << ";";
		}
		cout << endl;
	}
	*/
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
    size_t start_pos = str.find(from);
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

	// cout << sql_query << endl;

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
			string val(PQgetvalue(result, row, ii));
			if (!val.empty()) IOUtils::convertString(tmp(index[ii]), val);
		}
	}

	convertUnits(tmp);	
	vecMeteo.push_back(tmp);
}

void PSQLIO::map_parameters(PGresult* result, MeteoData& md, std::vector<size_t>& index)
{
	multiplier.clear();
	offset.clear();

	int columns = PQnfields(result);

	set<string> shadowed;
	map< string, set<string> >::iterator it = shadowed_parameters.find(md.meta.stationID);
	if (it != shadowed_parameters.end()) shadowed = it->second;

	for (int ii=0; ii<columns; ii++) {
		const string field_name(IOUtils::strToUpper(PQfname(result, ii)));
		//cout << "field(" << ii << "): " << field_name << endl;

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

bool PSQLIO::checkConsistency(const std::vector<MeteoData>& vecMeteo, StationData& sd)
{
	/**
	 * This function checks whether all the MeteoData elements in vecMeteo are consistent
	 * regarding their meta data (position information, station name). If they are consistent
	 * true is returned, otherwise false
	 */

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

void PSQLIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	//Loop through all stations
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		//1. check consistency of station data position -> write location in header or data section
		StationData sd;
		sd.position.setProj(coordout, coordoutparam);
		const bool isConsistent = checkConsistency(vecMeteo.at(ii), sd); // sd will hold valid meta info

		if (isConsistent) { //static station
			
		} else { //mobile station

		}
	}
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

void PSQLIO::convertUnits(MeteoData& meteo) const
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

	// For all parameters that have either an offset or an multiplier to bring to MKSA
	/*
	map<size_t, double>::const_iterator it;
	for (it = multiplier.begin(); it != multiplier.end(); it++) {
		double& tmp = meteo(it->first);
		if (tmp != IOUtils::nodata) tmp *= it->second;
	}

	for (it = offset.begin(); it != offset.end(); it++) {
		double& tmp = meteo(it->first);
		if (tmp != IOUtils::nodata) tmp += it->second;
	}
	*/
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
		cerr << "ERROR" << PQstatus(psql) << endl;
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
