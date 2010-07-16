/***********************************************************************************/
/*  Copyright 2009, 2010 WSL Institute for Snow and Avalanche Research   SLF-DAVOS */
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
#include "ImisIO.h"

using namespace std;
using namespace oracle;
using namespace oracle::occi;

namespace mio {
/**
 * @page imis IMIS
 * @section imis_format Format
 * This plugin reads data directly from the IMIS network database (Oracle database). It retrieves standard IMIS data as well as ENETZ and ANETZ data.
 *
 * @section imis_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 *
 * @section imis_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - DBNAME: name of the database to use (exemple: sdbo)
 * - DBUSER: user name to use when connecting to the database
 * - DBPASS: password to use when connecting to the database
 * - NROFSTATIONS: total number of stations listed for use
 * - STATION#: station code for the given number #
 */

const double ImisIO::plugin_nodata = -999.0; //plugin specific nodata value

const string ImisIO::sqlQueryMeteoData = "SELECT to_char(datum, 'YYYY-MM-DD HH24:MI') as datum, avg(ta) as ta, avg(iswr) as iswr, avg(vw) as vw, avg(dw) as dw, avg(rh) as rh, avg(ilwr) as ilwr, avg(hnw) as hnw, avg(tsg) as tsg, avg(tss) as tss, avg(hs) as hs, avg(rswr) as rswr from ams.v_ams_raw WHERE stat_abk=:1 and stao_nr=:2 and datum>=:3 and datum<=:4 group by datum order by datum asc";

const string ImisIO::sqlQueryStationData = "SELECT stao_name,stao_x,stao_y,stao_h from station2.standort WHERE stat_abk like :1 AND stao_nr=:2";

void ImisIO::getDBParameters()
{
	cfg.getValue("DBNAME", "Input", oracleDBName_in);
	cfg.getValue("DBUSER", "Input", oracleUserName_in);
	cfg.getValue("DBPASS", "Input", oraclePassword_in);
	/*cfg.getValue("DBNAME", "Output", oracleDBName_out);
	cfg.getValue("DBUSER", "Output", oracleUserName_out);
	cfg.getValue("DBPASS", "Output", oraclePassword_out);*/
}

ImisIO::ImisIO(void (*delObj)(void*), const std::string& filename) : IOInterface(delObj), cfg(filename)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getDBParameters();
}

ImisIO::ImisIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getDBParameters();
}

ImisIO::ImisIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getDBParameters();
}

ImisIO::~ImisIO() throw()
{
	cleanup();
}

void ImisIO::read2DGrid(Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readAssimilationData(const Date&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&,
                            const std::vector< std::vector<StationData> >&,
                            const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	if (vecMyStation.size() == 0)
		readStationMetaData(); //reads all the station meta data into the vecMyStation

	vecStation = vecMyStation;
}

/**
 * @brief A meta function that extracts all station names from the ConfigReader,
 *        parses them and retrieves all meta data from the IMIS database
 */
void ImisIO::readStationMetaData()
{
	vector<string> vecStationName;
	readStationNames(vecStationName);

	for (unsigned int ii=0; ii<vecStationName.size(); ii++){

		const string& stationName = vecStationName.at(ii);
		string stName = "", stationNumber = "";
		vector<string> resultset;

		//the stationName consists of the STAT_ABK and the STAO_NR, e.g. "KLO2" consists of "KLO" and "2"
		parseStationName(stationName, stName, stationNumber);

		//Now connect to the database and retrieve the meta data - this only needs to be done once per instance
		getStationData(stName, stationNumber, resultset);

		if (resultset.size() < 4)
			throw IOException("Could not read enough meta data", AT);

		double east, north, alt;
		if ((!IOUtils::convertString(east, resultset.at(1), std::dec))
		    || (!IOUtils::convertString(north, resultset.at(2), std::dec))
		    || (!IOUtils::convertString(alt, resultset.at(3), std::dec)))
			throw ConversionFailedException("Error while converting station coordinate from Imis DB", AT);

		Coords myCoord(coordin, coordinparam);
		myCoord.setXY(east, north, alt);
		vecMyStation.push_back(StationData(myCoord, stationName, resultset.at(0)));
	}
}

/**
 * @brief This function breaks up the station name into two components (a string and a number e.g. KLO2 -> "KLO","2")
 * @param stationName The full name of the station (e.g. "KLO2")
 * @param stName      The string part of the name  (e.g. "KLO")
 * @param stNumber    The integer part of the name (e.g. "2")
 */
void ImisIO::parseStationName(const std::string& stationName, std::string& stName, std::string& stNumber)
{
	stName    = stationName.substr(0, stationName.length()-1); //The station name: e.g. KLO
	stNumber  = stationName.substr(stationName.length()-1, 1); //The station number: e.g. 2
}

/**
 * @brief This function extracts all info about the stations that are to be used from global ConfigReader object
 * @param vecStationName A vector that will hold all relevant stations as std::strings
 */
void ImisIO::readStationNames(std::vector<std::string>& vecStationName)
{
	vecStationName.clear();

	//Read in the StationNames
	string xmlpath="", str_stations="";
	unsigned int stations=0;

	cfg.getValue("NROFSTATIONS", "Input", str_stations);

	if (str_stations == "")
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);

	if (!IOUtils::convertString(stations, str_stations, std::dec))
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);

	for (unsigned int ii=0; ii<stations; ii++) {
		stringstream tmp_stream;
		string stationname="", tmp_file="";

		tmp_stream << (ii+1); //needed to construct key name
		cfg.getValue(string("STATION"+tmp_stream.str()), "Input", stationname);
		std::cout << "\tRead io.ini stationname: '" << stationname << "'" << std::endl;
		vecStationName.push_back(stationname);
	}
}


void ImisIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                           std::vector< std::vector<MeteoData> >& vecMeteo,
                           std::vector< std::vector<StationData> >& vecStation,
                           const unsigned int& stationindex)
{
	if (vecMyStation.size() == 0)
		readStationMetaData(); //reads all the station meta data into the vecMyStation

	if (vecMyStation.size() == 0) //if there are no stations -> return
		return;

	unsigned int indexStart=0, indexEnd=vecMyStation.size();

	//The following part decides whether all the stations are rebuffered or just one station
	if (stationindex == IOUtils::npos){
		vecMeteo.clear();
		vecStation.clear();

		vecMeteo.insert(vecMeteo.begin(), vecMyStation.size(), vector<MeteoData>());
		vecStation.insert(vecStation.begin(), vecMyStation.size(), vector<StationData>());
	} else {
		if ((stationindex < vecMeteo.size()) && (stationindex < vecStation.size())){
			indexStart = stationindex;
			indexEnd   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("You tried to access a stationindex in readMeteoData that is out of bounds", AT);
		}
	}

	for (unsigned int ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		readData(dateStart, dateEnd, vecMeteo, vecStation, ii);
	}
}

/**
 * @brief A meta function to read meteo data for one specific station (specified by the stationindex)
 * @param dateStart    The beginning of the interval to retrieve data for
 * @param dateEnd      The end of the interval to retrieve data for
 * @param vecMeteo     The vector that will hold all MeteoData for each station
 * @param vecStation   The vector that will hold all StationData for each station
 * @param stationindex The index of the station as specified in the ConfigReader
 */
void ImisIO::readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
                      std::vector< std::vector<StationData> >& vecStation, const unsigned int& stationindex)
{
	vecMeteo.at(stationindex).clear();
	vecStation.at(stationindex).clear();

	string stationName="", stationNumber="";
	vector< vector<string> > vecResult;
	vector<int> datestart = vector<int>(5);
	vector<int> dateend   = vector<int>(5);

	parseStationName(vecMyStation.at(stationindex).getStationID(), stationName, stationNumber);

	dateStart.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	dateEnd.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);

	//Oracle can't deal with an integer for the hour of 24, hence the following workaround
	if (datestart[3] == 24){
		Date tmpDate = dateStart + Date(3.0/(60*60*24)); //add three seconds to omit 24 for 00 
		tmpDate.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	}
	if (dateend[3] == 24){
		Date tmpDate = dateEnd + Date(3.0/(60*60*24)); //add three seconds to omit 24 for 00 
		tmpDate.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);
	}

	getImisData(stationName, stationNumber, datestart, dateend, vecResult);

	MeteoData tmpmd;
	for (unsigned int ii=0; ii<vecResult.size(); ii++){
		parseDataSet(vecResult[ii], tmpmd);
		convertUnits(tmpmd);

		//Now insert tmpmd and a StationData object
		vecMeteo.at(stationindex).push_back(tmpmd);
		vecStation.at(stationindex).push_back(vecMyStation.at(stationindex));
	}
}

/**
 * @brief Puts the data that has been retrieved from the database into a MeteoData object
 * @param _meteo a row of meteo data from the database (note: order important, matches SQL query)
 * @param md     the object to copy the data to
 */
void ImisIO::parseDataSet(const std::vector<std::string>& _meteo, MeteoData& md)
{
	IOUtils::convertString(md.date, _meteo.at(0), dec);
	IOUtils::convertString(md.param(MeteoData::TA),   _meteo.at(1),  dec);
	IOUtils::convertString(md.param(MeteoData::ISWR), _meteo.at(2),  dec);
	IOUtils::convertString(md.param(MeteoData::VW),   _meteo.at(3),  dec);
	IOUtils::convertString(md.param(MeteoData::DW),   _meteo.at(4),  dec);
	IOUtils::convertString(md.param(MeteoData::RH),   _meteo.at(5),  dec);
	IOUtils::convertString(md.param(MeteoData::ILWR), _meteo.at(6),  dec);
	IOUtils::convertString(md.param(MeteoData::HNW),  _meteo.at(7),  dec);
	IOUtils::convertString(md.param(MeteoData::TSG),  _meteo.at(8),  dec);
	IOUtils::convertString(md.param(MeteoData::TSS),  _meteo.at(9),  dec);
	IOUtils::convertString(md.param(MeteoData::HS),   _meteo.at(10), dec);
	IOUtils::convertString(md.param(MeteoData::RSWR), _meteo.at(11), dec);
}

/**
 * @brief This function gets back data from table station2 and fills vector with station data
 * @param stat_abk :      a string key of table station2
 * @param stao_nr :       a string key of table station2
 * @param vecStationData: string vector in which data will be filled
 */
void ImisIO::getStationData(const std::string& stat_abk, const std::string& stao_nr, std::vector<std::string>& vecStationData)
{
	Environment *env = NULL;
	vecStationData.clear();

	try {
		Connection *conn = NULL;
		Statement *stmt = NULL;
		ResultSet *rs = NULL;

		env = Environment::createEnvironment();// static OCCI function
		conn = env->createConnection(oracleUserName_in, oraclePassword_in, oracleDBName_in);

		stmt = conn->createStatement(sqlQueryStationData);
		stmt->setString(1, stat_abk); // set 1st variable's value
		stmt->setString(2, stao_nr);  // set 2nd variable's value
		rs = stmt->executeQuery();    // execute the statement stmt

		while (rs->next() == true) {
			for (unsigned int ii=0; ii<4; ii++) {
				vecStationData.push_back(rs->getString(ii+1));
			}
		}

		stmt->closeResultSet(rs);
		conn->terminateStatement(stmt);
		env->terminateConnection(conn);

		Environment::terminateEnvironment(env); // static OCCI function
	} catch (exception& e){
		Environment::terminateEnvironment(env); // static OCCI function
		throw IOException("Oracle Error: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

/**
 * @brief This is a private function. It gets back data from ams.v_imis which is a table of the database
 * and fill them in a vector of vector of string. Each record returned is a string vector.
 * @param stat_abk :     a string key of ams.v_imis
 * @param stao_nr :      a string key of ams.v_imis
 * @param datestart :    a vector of five(5) integer corresponding to the recording date
 * @param dateend :      a vector of five(5) integer corresponding to the recording date
 * @param vecMeteoData : a vector of vector of string in which data will be filled
 */
void ImisIO::getImisData (const std::string& stat_abk, const std::string& stao_nr,
                          const std::vector<int>& datestart, const std::vector<int>& dateend,
                          std::vector< std::vector<std::string> >& vecMeteoData)
{
	Environment *env = NULL;
	vecMeteoData.clear();

	try {
		env = Environment::createEnvironment();// static OCCI function
		Connection *conn = NULL;
		Statement *stmt = NULL;
		ResultSet *rs = NULL;

		conn = env->createConnection(oracleUserName_in, oraclePassword_in, oracleDBName_in);
		stmt = conn->createStatement(sqlQueryMeteoData);

		// construct the oracle specific Date object: year, month, day, hour, minutes
		occi::Date begindate(env, datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
		occi::Date enddate(env, dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);
		stmt->setString(1, stat_abk); // set 1st variable's value (station name)
		stmt->setString(2, stao_nr);  // set 2nd variable's value (station number)
		stmt->setDate(3, begindate);  // set 3rd variable's value (begin date)
		stmt->setDate(4, enddate);    // set 4th variable's value (enddate)

		rs = stmt->executeQuery(); // execute the statement stmt

		rs->setMaxColumnSize(7,22);
		vector<string> vecTmpMeteoData;
		while (rs->next() == true) {
			vecTmpMeteoData.clear();
			for (unsigned int ii=1; ii<=12; ii++) { // 12 columns
				vecTmpMeteoData.push_back(rs->getString(ii));
			}
			vecMeteoData.push_back(vecTmpMeteoData);
		}

		stmt->closeResultSet(rs);
		conn->terminateStatement(stmt);
		env->terminateConnection(conn);
		Environment::terminateEnvironment(env); // static OCCI function
	} catch (exception& e){
		Environment::terminateEnvironment(env); // static OCCI function
		throw IOException("Oracle Error: " + string(e.what()), AT); //Translation of OCCI exception to IOException
	}
}

void ImisIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts ilwr to ea, converts RH to [0,1]
	if(meteo.ta!=IOUtils::nodata) {
		meteo.ta=C_TO_K(meteo.ta);
	}

	if(meteo.tsg!=IOUtils::nodata) {
		meteo.tsg=C_TO_K(meteo.tsg);
	}

	if(meteo.tss!=IOUtils::nodata) {
		meteo.tss=C_TO_K(meteo.tss);
	}

	if(meteo.rh!=IOUtils::nodata) {
		meteo.rh /= 100.;
	}
}

void ImisIO::cleanup() throw()
{
}

#ifndef _METEOIO_JNI
extern "C"
{
	//using namespace MeteoIO;
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const string& classname, const string& filename) {
		if(classname == "ImisIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new ImisIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
