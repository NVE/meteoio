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
#include "ImisIO.h"

const double ImisIO::plugin_nodata = -999.0; //plugin specific nodata value

using namespace oracle;
using namespace oracle::occi;

/**
 * @page imis IMIS
 * @section imis_format Format
 * This plugin reads data directly from the IMIS network database (Oracle database). It retrieves standard IMIS data as well as ENETZ data for the precipitations.
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
 * - COORDSYS: input coordinate system (see Coords)
 * - COORDPARAM: extra input coordinates parameters (see Coords)
 * - NROFSTATIONS: total number of stations listed for use
 * - STATION#: station code for the given number #
 */

/**
 * @class ImisIO 
 * @brief The class with-in the data from the database are treated. The MeteoData and the StationData will be set in.
 * This class also herited to IOInterface class which is abstract.
 * @author Moustapha Mbengue
 * @date 2009-05-12
 */

void ImisIO::getProjectionParameters() {
	//get projection parameters
	try {
		cfg.getValue("COORDSYS", "Input", coordsys);
		cfg.getValue("COORDPARAM", "Input", coordparam, ConfigReader::nothrow);
	} catch(std::exception& e){
		//problems while reading values for COORDSYS or COORDPARAM
		std::cerr << "[E] " << AT << ": reading configuration file: " << "\t" << e.what() << std::endl;
		throw;
	}
}

ImisIO::ImisIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename)
{
	getProjectionParameters();
}

ImisIO::ImisIO(const string& configfile) : IOInterface(NULL), cfg(configfile)
{
	getProjectionParameters();
}

ImisIO::ImisIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	getProjectionParameters();
}

ImisIO::~ImisIO() throw()
{
	cleanup();
}

void ImisIO::read2DGrid(Grid2DObject&, const string&)
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

void ImisIO::readAssimilationData(const Date_IO&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject&, const string&)
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

void ImisIO::readStationData(const Date_IO&, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	if (vecMyStation.size() == 0)
		readStationMetaData(); //reads all the station meta data into the vecMyStation

	vecStation = vecMyStation;
}

void ImisIO::readStationMetaData()
{
	vector<string> vecStationName;
	readStationNames(vecStationName);

	for (unsigned int ii=0; ii<vecStationName.size(); ii++){

		const string& stationName = vecStationName.at(ii);
		string stName = "";
		unsigned int stationNumber = 0;
		vector<string> resultset;

		//the stationName consists of the STAT_ABK and the STAO_NR
		parseStationName(stationName, stName, stationNumber);

		//Now connect to the database and retrieve the meta data - this only needs to be done once per instance
		getStation2Data(stName, stationNumber, resultset);

		if (resultset.size() < 4)
			throw IOException("Could not read enough meta data", AT);

		double east, north, alt;
		if ((!convertString(east, resultset.at(1), std::dec))
		    || (!convertString(north, resultset.at(2), std::dec))
		    || (!convertString(alt, resultset.at(3), std::dec)))
			throw ConversionFailedException("Error while converting station coordinate from Imis DB", AT);

		Coords myCoord(coordsys, coordparam);
		myCoord.setXY(east, north, alt);
		vecMyStation.push_back(StationData(myCoord, stationName));
	}
}

void ImisIO::parseStationName(const string& stationName, string& stName, unsigned int& stNumber)
{		
	stName    = stationName.substr(0, stationName.length()-1);
	string stNum  = stationName.substr(stationName.length()-1, 1);

	if (!convertString(stNumber, stNum))
		throw ConversionFailedException("Error while converting station number", AT);
}

void ImisIO::readStationNames(vector<string>& vecStationName)
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
		Date_IO tmp_date(0.0);
		
		tmp_stream << (ii+1); //needed to construct key name
		cfg.getValue(string("STATION"+tmp_stream.str()), "Input", stationname);
		std::cout << "\tRead io.ini stationname: '" << stationname << "'" << std::endl;
		vecStationName.push_back(stationname);
	}    
}


void ImisIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
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

void ImisIO::readData(const Date_IO& dateStart, const Date_IO& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, 
				 std::vector< std::vector<StationData> >& vecStation, const unsigned int& stationindex)
{
	vecMeteo.at(stationindex).clear();
	vecStation.at(stationindex).clear();

	unsigned int stationNumber;
	string stationName;
	vector< vector<string> > vecResult;
	vector<int> datestart = vector<int>(5);
	vector<int> dateend   = vector<int>(5);

	parseStationName(vecMyStation.at(stationindex).getStationName(), stationName, stationNumber);

	dateStart.getDate(datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]);
	dateEnd.getDate(dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]);

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
* @brief Puts the data that has been retrieved from the database into a MeteoBuffer
* which contains the meteo data and the station data of each single station in the configfile.
* @param meteo_in (vector \<vector \<string\>\>&) meteo data from the database.
* @param station_in (vector \<string\>&) station data from the database.
* @param mb (MeteoBuffer&) variable in which stationdata and meteodata are filled.
*/
void ImisIO::parseDataSet(const vector<string>& meteo_in, MeteoData& md)
{
	Date_IO tmpDate;
	double ta, iswr, vw, dw, rh, lwr, hnw, tsg, tss, hs, rswr;

	convertString(tmpDate, meteo_in.at(0), dec);
	convertString(ta,      meteo_in.at(1), dec);
	convertString(iswr,    meteo_in.at(2), dec);
	convertString(vw,      meteo_in.at(3), dec);
	convertString(dw,      meteo_in.at(4), dec);
	convertString(rh,      meteo_in.at(5), dec);
	convertString(lwr,     meteo_in.at(6), dec);
	convertString(hnw,     meteo_in.at(7), dec);
	convertString(tsg,     meteo_in.at(8), dec);
	convertString(tss,     meteo_in.at(9), dec);
	convertString(hs,      meteo_in.at(10), dec);
	convertString(rswr,    meteo_in.at(11), dec);
	
	md.setMeteoData(tmpDate, ta, iswr, vw, dw, rh, lwr, hnw, tsg, tss, hs, rswr);
}

/**
* @brief This is a private function. It gets back data from station2 which is a table of the database and fill them in a string vector
* @param stat_abk : a string key of station2 
* @param stao_nr : an integer key of station2
* @param data2S : string vector in which data will be filled
*/
void ImisIO::getStation2Data(const std::string stat_abk, unsigned int stao_nr, std::vector<std::string>& data2S)
{
	const string userName = "slf";
	const string password = "SDB+4u";
	const string dbName = "sdbo";
	unsigned int timeOut = 0, seconds = 60;

	Environment *env = Environment::createEnvironment();// static OCCI function
	{
		Connection *conn;
		Statement *stmt;
		ResultSet *rs;
		while (timeOut != 3) {
			timeOut = 0;
			try {
				conn = env->createConnection(userName, password, dbName);
				timeOut++;
			} catch (SQLException &connex) {
				cout <<"getStation2Data : Connection failed, please verify if userName, password and dbName are correct........"<< endl;
				cout << connex.getMessage();
				exit(1);
			}
			try {
				stmt = conn->createStatement("select stao_name,stao_x,stao_y,stao_h from station2.standort                              								where STAT_ABK =: 1 AND STAO_NR =: 2");
				stmt->setString(1, stat_abk); // set 1st variable's value
				stmt->setInt(2, stao_nr); // set 2nd variable's value 		
				rs = stmt->executeQuery(); // execute the statement stmt
				timeOut++;
			} catch (SQLException &stmtex) {
				cout <<"getStation2Data : Statement failed, please verify if it is correctly written............"<< endl;
				cout << stmtex.getMessage();
				exit(1);
			}
			try {			
				while (rs->next() == true) {
					for (int i=0; i<4; i++) {
						data2S.push_back(rs->getString(i+1));
					}
				}
				timeOut++;
			} catch (SQLException &rsex) {
				cout <<"getStation2Data : ResultSet manipulation failed, please verify if there is no mistake............."<< endl;
				cout << rsex.getMessage();
				exit(1);
			}catch (exception &cppex) { // C++ exception
				cout<< "Error "<< cppex.what()<<endl;
			}
			if (timeOut != 3 && seconds <= 27*60) {
				sleep(seconds);
				seconds *= 3;
			} else if (seconds > 27*60) {
				break;
			}	
		}   	   
		stmt->closeResultSet(rs);
		conn->terminateStatement(stmt);
		env->terminateConnection(conn);
	}
	Environment::terminateEnvironment(env); // static OCCI function
	
}

/**
* @brief This is a private function. It gets back data from ams.v_imis which is a table of the database
* and fill them in a vector of vector of string. It seems that each record is a string vector
* @param stat_abk : a string key of ams.v_imis
* @param stao_nr : an integer key of ams.v_imis
* @param date_in : a vector of five(5) integer corresponding to the recording date
* @param datatImis : a vector of vector of string in which data will be filled
*/
void ImisIO::getImisData (const string &stat_abk, const unsigned int &stao_nr, 
					 const vector<int>& datestart, const vector<int>& dateend, vector< vector<string> >& dataImis)
{
	const string userName = "slf";
	const string password = "SDB+4u";
	const string dbName = "sdbo";
	vector<string> vec;
	unsigned int timeOut = 0, seconds = 60;

	Environment *env = Environment::createEnvironment();// static OCCI function
	{
		Connection *conn;
		Statement *stmt;
		ResultSet *rs;
		while (timeOut != 3) {
			timeOut = 0;
			try {
				conn = env->createConnection(userName, password, dbName);
				timeOut++;
			} catch (SQLException &connex) {
				cout <<"getImisData : Connection failed, please verify if userName, password and dbName are correct........."<< endl;
				cout << connex.getMessage();
				exit(1);
			}
			try {
				stmt = conn->createStatement("select to_char(datum, 'YYYY-MM-DD HH24:MI') as datum,ta,iswr,vw,dw,rh,lwr,nswc,tsg,tss,hs,rswr from ams.v_amsio where STAT_ABK =: 1 AND STAO_NR =: 2 and DATUM >=: 3 and DATUM <=: 4 and rownum<=4800");
				// construct the oracle specific Date object: year, month, day, hour, minutes
				Date begindate(env, datestart[0], datestart[1], datestart[2], datestart[3], datestart[4]); 
				Date enddate(env, dateend[0], dateend[1], dateend[2], dateend[3], dateend[4]); 
				stmt->setString(1, stat_abk); // set 1st variable's value
				stmt->setInt(2, stao_nr); // set 2nd variable's value
				stmt->setDate(3, begindate); // set 3rd variable's value
				stmt->setDate(4, enddate); // set 4th variable's value

				rs = stmt->executeQuery(); // execute the statement stmt
				timeOut++;
			} catch (SQLException &stmtex) {
				cout <<"getImisData : Statement failed, please verify if it is correctly written............"<< endl;
				cout << stmtex.getMessage();
				exit(1);
			}
			try {		
				rs->setMaxColumnSize(7,22);
				while (rs->next() == true) {
					vec.clear();
					for (int i=1; i<=12; i++) { // 12 columns 
						vec.push_back(rs->getString(i));
					}
					dataImis.push_back(vec);
				}
				timeOut++;
			} catch (SQLException &rsex) {
				cout <<"getImisData : ResultSet manipulation failed, please verify if there is no mistake............."<< endl;
				cout << rsex.getMessage();
				exit(1);
			} catch (exception &cppex) { // C++ exception
				cout<< "Error "<< cppex.what()<<endl;
			}
			if (timeOut != 3 && seconds <= 27*60) {
				sleep(seconds);
				seconds *= 3;
			} else if (seconds > 27*60) {
				break;
			}
		}   	   
		stmt->closeResultSet(rs);
		conn->terminateStatement(stmt);
		env->terminateConnection(conn);
	}
	Environment::terminateEnvironment(env); // static OCCI function
}

void ImisIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts lwr to ea, converts RH to [0,1]
	if(meteo.ta!=IOUtils::nodata) {
		meteo.ta=C_TO_K(meteo.ta);
	}
	
	if(meteo.tsg!=IOUtils::nodata) {
		meteo.tsg=C_TO_K(meteo.tss);
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
