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

ImisIO::ImisIO(const string& configfile) : IOInterface(NULL), cfg(configfile)
{
	mbImis.clear();
	vecStationName.clear();
}

ImisIO::ImisIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

ImisIO::~ImisIO() throw()
{
	cleanup();
}

void ImisIO::read2DGrid(Grid2DObject& grid_out, const string& parameter)
{
	//Nothing so far
	(void)parameter;
	(void)grid_out;
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readDEM(DEMObject& dem_out)
{
	//Nothing so far
	(void)dem_out;
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readLanduse(Grid2DObject& landuse_out)
{
	//Nothing so far
	(void)landuse_out;
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	//Nothing so far
	(void)date_in;
	(void)da_out;
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::readSpecialPoints(CSpecialPTSArray& pts)
{
	//Nothing so far
	(void)pts;
	throw IOException("Nothing implemented here", AT);
}

void ImisIO::write2DGrid(const Grid2DObject& grid_in, const string& name)
{
	//Nothing so far
	(void)grid_in;
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

/**
* @brief refer to void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
*/
void ImisIO::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo)
{
	vector<StationData> vecStation;
	readMeteoData(date_in, vecMeteo, vecStation);
}

/**
* @brief this method is used to get back meteo and station data for a given day
* into the respective vectors MeteoData and StationData
* @param date_in (Date_IO&) recording date
* @param vecMeteo (vector \<MeteoData\>&) vector in which meteo data will be filled
* @param vecStation (vector \<StationData\>&) vector in which station data will be filled
*/
void ImisIO::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{
	if (vecStationName.size() == 0)
		readStationNames(); //reads station names into vector<string> vecStationName

	if (vecStationName.size() == 0) //if there are no stations -> return
		return;

	unsigned int indexStart=0, indexEnd=vecStationName.size();

	//The following part decides whether all the stations are rebuffered or just one station
	if (stationindex == IOUtils::npos){
		vecMeteo.clear();
		vecStation.clear();

		vecMeteo.insert(vecMeteo.begin(), vecStationName.size(), vector<MeteoData>());
		vecStation.insert(vecStation.begin(), vecStationName.size(), vector<StationData>());
	} else {
		if ((stationindex < vecMeteo.size()) && (stationindex < vecStation.size())){
			indexStart = stationindex;
			indexEnd   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("You tried to access a stationindex in readMeteoData that is out of bounds", AT);
		}
	}

	for (unsigned int ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		StationData sd;
		readStationMetaData(sd, ii);

		readData(dateStart, dateEnd, vecMeteo, vecStation, sd, ii);
	}
/*	vecMeteo.clear();
	vecStation.clear();
	
	if (mbImis.size() == 0) {
		getStationName();
		createBuffer();
	}
	
	unsigned int size = mbImis.size();
	
	for(unsigned int ii=0; ii<size; ii++) {
		MeteoData md;
		md.date = date_in;
		StationData sd;
		if (mbImis[ii].seek(date_in) == MeteoBuffer::npos) {
			setMbImis(date_in, vecStationName[ii], mbImis[ii]);
			resampleMbImis(md, sd, date_in, mbImis[ii]);
		} else {
			if (mbImis[ii].getMeteoData(ii).date <= date_in && mbImis[ii].getMeteoData(mbImis[ii].size()-1).date >= date_in) {
				//cerr << "[I] Buffered data found for date: " << date_in.toString() << endl;
				resampleMbImis(md, sd, date_in, mbImis[ii]);
			}
		}
		vecMeteo.push_back(md);
		vecStation.push_back(sd);
	}

	if (vecMeteo.size() == 0) {//No data found
		throw IOException("[E] No data for any station for date " + date_in.toString() + " found", AT);
	}*/

}

/**
* @brief Puts the data that has been retrieved from the database into a MeteoBuffer
* which contains the meteo data and the station data of each single station in the configfile.
* @param meteo_in (vector \<vector \<string\>\>&) meteo data from the database.
* @param station_in (vector \<string\>&) station data from the database.
* @param mb (MeteoBuffer&) variable in which stationdata and meteodata are filled.
*/
void ImisIO::createData(vector< vector<string> >& meteo_in, vector<string>& station_in, MeteoBuffer& mb)
{
	MeteoData md;
	StationData sd;
	Date_IO tmpDate;
	
	double east, north, lat, lon, alt;
	convertString(east, station_in[1], dec);
	convertString(north, station_in[2], dec);
	convertString(alt, station_in[3], dec);
	string sName = "";
	if (station_in[0].size() == 0) {
		sName = vecStationName[mbImis.size()];
	} else {
		sName = station_in[0];
	}
	CH1903_to_WGS84(east, north, lat, lon);
	sd.setStationData(east, north, alt, sName, lat, lon);
	
	double ta, iswr, vw, dw, rh, lwr, hnw, tsg, tss, hs, rswr;
	const unsigned int size = meteo_in.size();
	for (unsigned int i=0; i<size; i++) {
		//loop over the stations
		convertString(tmpDate, meteo_in[i][0], dec);
		convertString(ta, meteo_in[i][1], dec);
		convertString(iswr, meteo_in[i][2], dec);
		convertString(vw, meteo_in[i][3], dec);
		convertString(dw, meteo_in[i][4], dec);
		convertString(rh, meteo_in[i][5], dec);
		convertString(lwr, meteo_in[i][5], dec);
		convertString(hnw, meteo_in[i][7], dec);
		convertString(tsg, meteo_in[i][8], dec);
		convertString(tss, meteo_in[i][9], dec);
		convertString(hs, meteo_in[i][10], dec);
		convertString(rswr, meteo_in[i][11], dec);
		md.setMeteoData(tmpDate, ta, iswr, vw, dw, rh, lwr, hnw, tsg, tss, hs, rswr);
		
		mb.put(md, sd);
	}
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
	const string password = "sdb+4u";
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
void ImisIO::getImisData (const string &stat_abk, const unsigned int &stao_nr, vector<int> date_in, vector< vector<string> >& dataImis)
{
	const string userName = "slf";
	const string password = "sdb+4u";
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
				if (stao_nr != 0) {
					stmt = conn->createStatement("select to_char(datum, 'YYYY-MM-DD HH24:MI') as datum,ta,iswr,vw,dw,rh,lwr,nswc, 									    tsg,tss,hs,rswr from ams.v_amsio where STAT_ABK =: 1 AND STAO_NR =: 2             									    and DATUM >=: 3 and rownum<=4800");
					Date edate(env, date_in[0], date_in[1], date_in[2], date_in[3], date_in[4]); // year, month, day, hour, minutes
					stmt->setString(1, stat_abk); // set 1st variable's value
					stmt->setInt(2, stao_nr); // set 2nd variable's value
					stmt->setDate(3, edate); // set 3rd variable's value
				} else {
					string sql = "select to_char(datum, 'YYYY-MM-DD HH24:MI') as datum,ta,iswr,vw,dw,rh,lwr,nswc,tsg,tss,hs,rswr 								from ams.v_amsio where STAT_ABK=:1 AND STAO_NR is null and DATUM>=:2 and rownum<=4800";
					stmt = conn->createStatement(sql);
					Date edate(env, date_in[0], date_in[1], date_in[2], date_in[3], date_in[4]); // year, month, day, hour, minutes
					stmt->setString(1, stat_abk); // set 1st variable's value
					stmt->setDate(2, edate); // set 2nd variable's value				
				}
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

/**
* @brief Get back station's name from the configfile
*/
void ImisIO::getStationName()
{
	vecStationName.clear();
	
	string str_stations="";
	int stations=0;
	
	cfg.getValue("NROFSTATIONS", str_stations);
	
	if (!convertString(stations, str_stations, dec)) {
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);
	}
	for (int ii=0; ii<stations; ii++) {
		stringstream tmp_stream;
		string stationname="";
		
		tmp_stream << (ii+1); //needed to construct key name
		cfg.getValue(string("STATION"+tmp_stream.str()), stationname);
		
		vecStationName.push_back(stationname);
	}
}
		
void ImisIO::setMbImis(Date_IO date_in, const string& stationName, MeteoBuffer& buffer)
{
	buffer.clear();
	int date[5], stao;
	date_in -= 1./24.; // one hour before
	date_in += 1./(24.*3600.*100.); //Oracle does not want 24:00, so we must make sure we use 00:00 instead
	date_in.getDate(date[0],date[1],date[2],date[3],date[4]);
	vector<int> date_io(date, date+sizeof(date)/sizeof(int));
	string station, name, number;
	
	vector<string>* data2s = new vector<string>;
	vector< vector<string> >* data_imis = new vector< vector<string> >;
	station = stationName;
	number = station.substr(station.length()-1);
	name = station.substr(0, station.length()-1);

	if (!convertString(stao, number, dec)) {
		throw ConversionFailedException("Error while reading station number in readMeteoData(...) ", AT);
	}

	getStation2Data(name, stao, *data2s);
	getImisData(name, stao, date_io, *data_imis);
	createData(*data_imis, *data2s, buffer);
	delete(data2s);
	delete(data_imis);
}

void ImisIO::resampleMbImis(MeteoData& meteo, StationData& station, const Date_IO& date_in, MeteoBuffer& mb)
{
	unsigned int index = mb.seek(date_in);
	if (index != MeteoBuffer::npos) {
		if (mb.getMeteoData(index).date == date_in) {
			meteo = mb.getMeteoData(index);
			station = mb.getStationData(index);
		} else {
			Meteo1DResampler mresampler;
			mresampler.resample(index, date_in, mb);
			if (index != MeteoBuffer::npos) {
				meteo = mb.getMeteoData(index);
				station = mb.getStationData(index);
			} else {
				cout << "[i] Buffering data for Station " << station.stationName << " at date " 
					<< date_in.toString() << " failed" << endl;
			}
		}
	}
}

/**
* @brief ImisIO cleaner function
*/
void ImisIO::cleanup() throw()
{
	mbImis.clear();
	vecStationName.clear();
}

void ImisIO::createBuffer()
{
	//Clear the buffers
	mbImis.clear();
	const unsigned int stations = vecStationName.size();

	//Allocate one MeteoBuffer per station
	for (unsigned int ii=0; ii<stations; ii++) {
		mbImis.push_back(MeteoBuffer(IMIS_BUFF_SIZE));
	}
	cout << "[I] "<<AT<<": Created Buffer for " << stations << " stations" << endl;
}

/**
* @brief Returns the ConfigFile
*/
ConfigReader ImisIO::getCfg()
{
	return cfg;
}

/**
* @brief Returns vecStationName, a string vector in which StationName are saved
*/
vector<string> ImisIO::getVecStationName()
{
	return vecStationName;
}

vector<MeteoBuffer> ImisIO::getMbImis()
{
	return mbImis;
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
