// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/WWCSIO.h>
#include <meteoio/plugins/libMysqlWrapper.h>

#ifdef _WIN32
	#include <winsock.h>
#endif // _WIN32

#include <mysql.h>
#include <stdio.h>
//#include <cstring>
#include <algorithm>

using namespace std;

namespace mio {
/**
* @page wwcs WWCSIO
* @section WWCS_format Format
* This is the plugin required to get meteorological data from the WWCS MySQL database. The connection 
* to the database is encrypted (with SSL) and compressed (with zlib) if supported by the server.
*
* @section WWCS_units Units
* The units are assumed to be the following:
* - __temperatures__ in celsius
* - __relative humidity__ in %
* - __pressures__ in Pa
*
* @section WWCS_keywords Keywords
* This plugin uses the following keywords:
* - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
* - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
* - TIME_ZONE: For [Input] and [Output] sections (Input::TIME_ZONE should describe the timezone 
* of the data in the database while the resulting data will be converted to Output::TIME_ZONE )
* - WWCS_HOST: MySQL Host Name (e.g. localhost or 191.168.145.20); [Input] section
* - WWCS_DB: MySQL Database (e.g. wwcs); [Input] section
* - WWCS_USER: MySQL User Name (e.g. wwcs); [Input] section
* - WWCS_PASS: MySQL password; [Input] section
* - STATION#: station code for the given number #; [Input] section
* 
* Example configuration:
* @code
* [Input]
* COORDSYS  = CH1903
* TIME_ZONE = 1
* 
* METEO     = WWCS
* WWCS_HOST = nesthorn.slf.ch
* WWCS_DB   = wwcs
* WWCS_USER = wwcs
* WWCS_PASS = XXX
* 
* STATION1  = SLF01
* STATION2  = WWCS_BAL01
* STATION3  = WWCS_LUC
* @endcode
* 
* 
*/

const string WWCSIO::MySQLQueryStationMetaData = "SELECT stationName, latitude, longitude, altitude, slope, azimuth FROM sites WHERE StationID=?";
const string WWCSIO::MySQLQueryMeteoData = "SELECT timestamp, ta, rh, p FROM meteoseries WHERE loggerID=? and timestamp>=? AND timestamp<=? ORDER BY timestamp ASC";

WWCSIO::WWCSIO(const std::string& configfile)
        : cfg(configfile), vecStationIDs(), vecStationMetaData(),
          mysqlhost(), mysqldb(), mysqluser(), mysqlpass(),
          coordin(), coordinparam(), coordout(), coordoutparam(),
          in_dflt_TZ(1.), out_dflt_TZ(1.), mysql_options(mysql_wrp::COMPRESSION | mysql_wrp::ENCRYPTION)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	readConfig();
}

WWCSIO::WWCSIO(const Config& cfgreader)
        : cfg(cfgreader), vecStationIDs(), vecStationMetaData(),
          mysqlhost(), mysqldb(), mysqluser(), mysqlpass(),
          coordin(), coordinparam(), coordout(), coordoutparam(),
          in_dflt_TZ(1.), out_dflt_TZ(1.), mysql_options(mysql_wrp::COMPRESSION | mysql_wrp::ENCRYPTION)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	readConfig();
}

void WWCSIO::readConfig()
{
	cfg.getValue("TIME_ZONE","Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE","Output", out_dflt_TZ, IOUtils::nothrow);
	
	cfg.getValue("WWCS_HOST", "Input", mysqlhost);
	cfg.getValue("WWCS_DB", "Input", mysqldb);
	cfg.getValue("WWCS_USER", "Input", mysqluser);
	cfg.getValue("WWCS_PASS", "Input", mysqlpass);
}

/**
* @brief Build the list of user provided station IDs to read
* @return list of station IDs
*/
std::vector<std::string> WWCSIO::readStationIDs() const
{
	std::vector<std::string> vecStationID;
	cfg.getValues("STATION", "INPUT", vecStationID);

	if (vecStationID.empty())
		cerr << "\tNo stations specified for WWCSIO... is this what you want?\n";
	
	return vecStationID;
}

/**
* @brief Retrieve the stations' metadata from the database.
* This metadata can then be used either to return only the metadata as a vector of StationData or
* later when reading the meteo data, in order to provide the necessary metadata.
* The results are stored in the vecStationMetaData member.
*/
void WWCSIO::readStationMetaData()
{
	vecStationMetaData.clear();
	const std::vector<std::string> vecStationID( readStationIDs() );
	
	MYSQL *mysql = mysql_wrp::initMysql(mysqlhost, mysqluser, mysqlpass, mysqldb, mysql_options);
	MYSQL_STMT *stmt = mysql_wrp::initStmt(&mysql, MySQLQueryStationMetaData, 1);
	std::vector<mysql_wrp::fType> result_fields{ mysql_wrp::fType(MYSQL_TYPE_STRING), mysql_wrp::fType(MYSQL_TYPE_DOUBLE), mysql_wrp::fType(MYSQL_TYPE_DOUBLE), mysql_wrp::fType(MYSQL_TYPE_DOUBLE), mysql_wrp::fType(MYSQL_TYPE_DOUBLE), mysql_wrp::fType(MYSQL_TYPE_DOUBLE) };
	
	for (size_t ii=0; ii<vecStationID.size(); ii++) {
		const std::string stationID( vecStationID[ii] );
		std::vector<mysql_wrp::fType> params_fields{ mysql_wrp::fType(stationID)};
		mysql_wrp::bindParams(&stmt, params_fields);
		
		if (mysql_stmt_execute(stmt)) {
			throw IOException("Error executing statement", AT);
		} else { //retrieve results
			mysql_wrp::bindResults(&stmt, result_fields);
			if (mysql_stmt_num_rows(stmt)!=1) throw IOException("stationID is not unique in the database!", AT);
			mysql_stmt_fetch(stmt); //we only have one result, see check above
			
			const std::string station_name( result_fields[0].str );
			Coords location(coordin,coordinparam);
			location.setLatLon(result_fields[1].val, result_fields[2].val, result_fields[3].val);
			StationData sd(location, stationID, station_name);
			sd.setSlope(result_fields[4].val, result_fields[5].val);
			vecStationMetaData.push_back( sd );
		}
	}
	
	if (mysql_stmt_close(stmt)) {
		throw IOException("Failed closing Mysql connection: "+std::string(mysql_error(mysql)), AT);
	}
	
	mysql_close(mysql);
}

//standard call to get the metadata as defined in IOInterface.h
void WWCSIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	readStationMetaData(); //reads all the station meta data into the vecStationMetaData (member vector)
	vecStation = vecStationMetaData; //vecStationMetaData is a global vector holding all meta data
}

//standard call to get the meteo data as defined in IOInterface.h
void WWCSIO::readMeteoData(const Date& dateStart , const Date& dateEnd,
                            std::vector<std::vector<MeteoData> >& vecMeteo)
{
	if (vecStationMetaData.empty()) readStationMetaData();

	vecMeteo.clear();
	vecMeteo.insert(vecMeteo.begin(), vecStationMetaData.size(), vector<MeteoData>());

	for (size_t ii=0; ii<vecStationMetaData.size(); ii++) { //loop through relevant stations
		readData(dateStart, dateEnd, vecMeteo, ii);
	}
}

/**
* @brief Retrieve the meteo data for one given station from the database
* @param[in] dateStart start date of the data to retrieve
* @param[in] dateEnd end date of the data to retrieve
* @param[out] vecMeteo vector to store the data for all stations
* @param[in] stationindex index of the current station of interest in both vecStationMetaData and vecMeteo
*/
void WWCSIO::readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
                       const size_t& stationindex) const
{
	vecMeteo.at(stationindex).clear();

	Date dateS(dateStart), dateE(dateEnd);
	dateS.setTimeZone(in_dflt_TZ);
	dateE.setTimeZone(in_dflt_TZ);
	
	MYSQL *mysql = mysql_wrp::initMysql(mysqlhost, mysqluser, mysqlpass, mysqldb, mysql_options);
	MYSQL_STMT *stmt = mysql_wrp::initStmt(&mysql, MySQLQueryMeteoData, 3);
	std::vector<mysql_wrp::fType> result_fields{ mysql_wrp::fType(MYSQL_TYPE_DATETIME), mysql_wrp::fType(MYSQL_TYPE_DOUBLE), mysql_wrp::fType(MYSQL_TYPE_DOUBLE), mysql_wrp::fType(MYSQL_TYPE_DOUBLE) };
	
	const StationData sd( vecStationMetaData[stationindex] );
	const std::string stationID( sd.getStationID() );
	//const std::string stationID( "98:f4:ab:39:96:90" ); //HACK for debuging since we don't have the view yet
	std::vector<mysql_wrp::fType> params_fields{ mysql_wrp::fType(stationID), mysql_wrp::fType(dateStart), mysql_wrp::fType(dateEnd)};
	mysql_wrp::bindParams(&stmt, params_fields);
	
	if (mysql_stmt_execute(stmt)) {
		throw IOException("Error executing statement", AT);
	} else { //retrieve results
		mysql_wrp::bindResults(&stmt, result_fields);
		
		do {
			const int status = mysql_stmt_fetch(stmt);
			if (status==1 || status==MYSQL_NO_DATA)
				break;
			
			if (result_fields[0].is_null==1) continue; //this should not happen, but better safe than sorry!
			
			MeteoData md( result_fields[0].getDate(in_dflt_TZ), sd); //get from Mysql without TZ, set it to input TZ
			md.date.setTimeZone(out_dflt_TZ); //set to requested TZ
			
			md("TA") = retrieveData(result_fields[1], mysql_wrp::C_TO_K);
			md("RH") = retrieveData(result_fields[2], mysql_wrp::NORMALIZE_PC);
			md("P") = retrieveData(result_fields[3]);
			
			vecMeteo[stationindex].push_back( md );
		} while (true);
	}
	
	if (mysql_stmt_close(stmt)) {
		throw IOException("Failed closing Mysql connection: "+std::string(mysql_error(mysql)), AT);
	}
	
	mysql_close(mysql);
}

} //namespace
