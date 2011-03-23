/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include "SMETIO.h"
#include <meteoio/IOUtils.h>

using namespace std;

namespace mio {
/**
 * @page smetio SMET
 * @section template_format Format
 * The Station METeo data files is a station centered, ascii file format that has been designed with flexibility and ease of use in mind. Please refer to its <a href="../SMET_specifications.pdf">official format specification</a> for more information.
 *
 * @section template_units Units
 * All units are MKSA, the only exception being the precipitations that are in mm/h. It is however possible to use  multipliers and offsets (but they must be specified in the file header).
 *
 * @section template_keywords Keywords
 * This plugin uses the following keywords:
 * - STATION#: input filename (in METEOPATH). As many meteofiles as needed may be specified
 * - METEOPATH: meteo files directory where to read/write the meteofiles; [Input] and [Output] sections
 * - METEOPARAM: output file format options (ASCII or BINARY that might be followed by GZIP)
 *
 * Example:
 * @code
 * [Input]
 * METEOPATH = ./input
 * STATION1 = uppper_station.smet
 * STATION2 = lower_station.smet
 * STATION3 = outlet_station.smet
 * [Output]
 * METEOPATH = ./output
 * METEOPARAM = ASCII GZIP
 * @endcode
 */

//TODO: keep a pointer to last read position, so we can fseek for the next read?

const std::string SMETIO::smet_version = "1.1";
const unsigned int SMETIO::buffer_reserve = 23*24*2; //kind of average size of a buffer for optimizing vectors
map<string, MeteoData::Parameters> SMETIO::mapParameterByName;
const bool SMETIO::__init = SMETIO::initStaticData();

bool SMETIO::initStaticData()
{
	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){
		mapParameterByName[MeteoData::getParameterName(ii)] = MeteoData::Parameters(ii);
	}

	return true;
}

double& SMETIO::getParameter(const std::string& columnName, MeteoData& md)
{//HACK: the whole name mapping is a big hack. Replace with proper mapping!!!
	if(columnName=="OSWR") {
		MeteoData::Parameters paramindex = mapParameterByName["RSWR"];
		return md.param(paramindex);
	}
	if(columnName=="PSUM") {
		MeteoData::Parameters paramindex = mapParameterByName["HNW"];
		return md.param(paramindex);
	}

	MeteoData::Parameters paramindex = mapParameterByName[columnName];
	return md.param(paramindex);
}

void SMETIO::checkColumnNames(const std::vector<std::string>& vecColumns, const bool& locationInHeader)
{
	/*
	 * This function checks whether the sequence of keywords specified in the 
	 * [HEADER] section (key 'fields') is valid
	 */
	vector<unsigned int> paramcounter = vector<unsigned int>(MeteoData::nrOfParameters, 0);

	for (unsigned int ii=0; ii<vecColumns.size(); ii++){
		std::string column = vecColumns[ii];

		//column names mapping
		if(column=="OSWR") column="RSWR";
		if(column=="PSUM") column="HNW";

		if ((column == "timestamp") || (column == "longitude")
		    || (column == "latitude") || (column == "altitude")){
			//everything ok
		} else {
			map<string, MeteoData::Parameters>::iterator it = mapParameterByName.find(column);
			if (it == mapParameterByName.end())
				throw InvalidFormatException("Key 'fields' specified in [HEADER] section contains invalid names", AT);
			
			paramcounter.at(mapParameterByName[column])++;
		}
	}
	
	//Check for multiple usages of parameters
	for (unsigned int ii=0; ii<paramcounter.size(); ii++){
		if (paramcounter[ii] > 1)
			throw InvalidFormatException("In 'fields': Multiple use of " + MeteoData::getParameterName(ii), AT);
	}
	
	//If there is no location information in the [HEADER] section, then
	//location information must be part of fields
	if (!locationInHeader){
		unsigned int latcounter = 0, loncounter=0, altcounter=0;
		for (unsigned int ii=0; ii<vecColumns.size(); ii++){
			if (vecColumns[ii] == "longitude")
				loncounter++;
			else if (vecColumns[ii] == "latitude")
				latcounter++;
			else if (vecColumns[ii] == "altitude")
				altcounter++;
		}

		if ((latcounter != loncounter) && (loncounter != altcounter) && (altcounter != 1))
			throw InvalidFormatException("Key 'fields' must contain 'latitude', 'longitude' and 'altitude'", AT);
	}
}

//END STATIC SECTION


SMETIO::SMETIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	parseInputOutputSection();
	plugin_nodata = IOUtils::nodata;
}

SMETIO::SMETIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	parseInputOutputSection();
	plugin_nodata = IOUtils::nodata;
}

SMETIO::SMETIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	parseInputOutputSection();
	plugin_nodata = IOUtils::nodata;
}

SMETIO::~SMETIO() throw()
{
	cleanup();
}

void SMETIO::cleanup() throw()
{
	//clear ios flags
	fout << resetiosflags(ios_base::fixed | ios_base::left);

	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void SMETIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{//big HACK: this is a barbaric code duplication!! Plus it should support coordinates in the data
//ie: it should use the given date! (and TZ)
	vecStation.clear();
	vecStation.reserve(nr_stations);

	//Now loop through all requested stations, open the respective files and parse them
	for (unsigned int ii=0; ii<nr_stations; ii++){
		bool isAscii = true;
		string filename = vecFiles.at(ii); //filename of current station
		
		if (!IOUtils::fileExists(filename))
			throw FileNotFoundException(filename, AT);

		fin.clear();
		fin.open (filename.c_str(), ios::in);
		if (fin.fail())
			throw FileAccessException(filename, AT);

		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		//Go through file, save key value pairs
		string line="";
		std::vector<std::string> tmpvec, vecDataSequence;
		std::vector<double> vecUnitsOffset, vecUnitsMultiplier;
		double timezone = 0.0;
		bool locationInHeader = false;
		StationData sd;

		try {
			//1. Read signature
			getline(fin, line, eoln); //read complete signature line
			IOUtils::stripComments(line);
			IOUtils::readLineToVec(line, tmpvec);
			checkSignature(tmpvec, filename, isAscii);

			//2. Read Header
			readHeader(eoln, filename, locationInHeader, timezone, sd, vecDataSequence, vecUnitsOffset, vecUnitsMultiplier);
			cleanup();
		} catch(std::exception& e) {
			cleanup();
			throw;
		}
		vecStation.push_back(sd);
	}
}

void SMETIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE","Input",in_dflt_TZ,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_dflt_TZ,Config::nothrow);

	// Parse the [Input] and [Output] sections within Config object cfg
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

	//Parse input section: extract number of files to read and store filenames in vecFiles
	std::string inpath="", in_meteo="";
	cfg.getValue("METEO", "Input", in_meteo, Config::nothrow);
	if (in_meteo == "SMET") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("METEOPATH", "Input", inpath);
		unsigned int counter = 1;
		string filename = "";

		do {
			stringstream ss;
			filename = "";
			
			ss << "STATION" << counter;
			cfg.getValue(ss.str(), "Input", filename, Config::nothrow);
			
			if (filename != ""){
				stringstream file_and_path;
				file_and_path << inpath << "/" << filename;
				if (!IOUtils::validFileName(file_and_path.str())) //Check whether filename is valid
					throw InvalidFileNameException(file_and_path.str(), AT);
				
				vecFiles.push_back(file_and_path.str());
			}
			counter++;
		} while (filename != "");
		
		nr_stations = counter - 1;
	}

	//Parse output section: extract info on whether to write ASCII or BINARY format, gzipped or not
	outpath = "";
	outputIsAscii = true; 
	outputIsGzipped = false;

	vector<string> vecArgs;
	cfg.getValue("METEOPATH", "Output", outpath, Config::nothrow);
	cfg.getValue("METEOPARAM", "Output", vecArgs, Config::nothrow); //"ASCII|BINARY GZIP"

	if (outpath == "")
		return;

	if (vecArgs.size() == 0)
		vecArgs.push_back("ASCII");

	if (vecArgs.size() > 2)
		throw InvalidFormatException("Too many values for key METEOPARAM", AT);

	if (vecArgs[0] == "BINARY")
		outputIsAscii = false;
	else if (vecArgs[0] == "ASCII")
		outputIsAscii = true;
	else 
		throw InvalidFormatException("The first value for key METEOPARAM may only be ASCII or BINARY", AT);

	if (vecArgs.size() == 2){
		if (vecArgs[1] != "GZIP")
			throw InvalidFormatException("The second value for key METEOPARAM may only be GZIP", AT);

		outputIsGzipped = true;
	}

}

void SMETIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                           std::vector< std::vector<MeteoData> >& vecMeteo,
                           const unsigned int& stationindex)
{
	//Make sure that vecMeteo have the correct dimension and stationindex is valid
	unsigned int startindex=0, endindex=vecFiles.size();	
	if (stationindex != IOUtils::npos){
		if ((stationindex < vecFiles.size()) || (stationindex < vecMeteo.size())){
			startindex = stationindex;
			endindex = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("Invalid stationindex", AT);
		} 

		vecMeteo[stationindex].clear();
	} else {
		vecMeteo.clear();
		vecMeteo = vector< vector<MeteoData> >(vecFiles.size());
		vecMeteo.reserve(nr_stations);
	}

	//Now loop through all requested stations, open the respective files and parse them
	for (unsigned int ii=startindex; ii<endindex; ii++){
		bool isAscii = true;
		string filename = vecFiles.at(ii); //filename of current station
		
		if (!IOUtils::fileExists(filename))
			throw FileNotFoundException(filename, AT);

		fin.clear();
		fin.open (filename.c_str(), ios::in);
		if (fin.fail())
			throw FileAccessException(filename, AT);

		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		//Go through file, save key value pairs
		string line="";
		std::vector<std::string> tmpvec, vecDataSequence;
		std::vector<double> vecUnitsOffset, vecUnitsMultiplier;
		double timezone = 0.0;
		bool locationInHeader = false;
		StationData sd;

		try {
			//1. Read signature
			getline(fin, line, eoln); //read complete signature line
			IOUtils::stripComments(line);
			IOUtils::readLineToVec(line, tmpvec);
			checkSignature(tmpvec, filename, isAscii);

			//2. Read Header
			readHeader(eoln, filename, locationInHeader, timezone, sd, vecDataSequence, vecUnitsOffset, vecUnitsMultiplier);
			SMETIO::checkColumnNames(vecDataSequence, locationInHeader);

			//now, timezone MUST contain something valid
			if(timezone==IOUtils::nodata) {
				stringstream ss;
				ss << "No timezone information available for file " << filename;
				ss << " either in file header or io.ini";
				throw InvalidFormatException(ss.str(), AT);
			}

			//3. Read DATA
			if (isAscii){
				readDataAscii(eoln, filename, timezone, sd, vecDataSequence, vecUnitsOffset,
			                      vecUnitsMultiplier, dateStart, dateEnd, vecMeteo[ii]);
			} else {
				streampos currpos = fin.tellg();
				fin.close();
				fin.open (filename.c_str(), ios::in|ios::binary);				
				if (fin.fail()) throw FileAccessException(filename, AT);
				fin.seekg(currpos); //jump to binary data section

				readDataBinary(eoln, filename, timezone, sd, vecDataSequence, vecUnitsOffset,
				               vecUnitsMultiplier, dateStart, dateEnd, vecMeteo[ii]);
			}
			cleanup();
		} catch(std::exception& e) {
			cleanup();
			throw;
		}
	}
}

void SMETIO::readDataBinary(const char&, const std::string&, const double& timezone,
                            const StationData& sd, const std::vector<std::string>& vecDataSequence,
                            const std::vector<double>& vecUnitsOffset, std::vector<double>& vecUnitsMultiplier,
                            const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo)
{
	const unsigned int nrOfColumns = vecDataSequence.size();

	vecMeteo.reserve(buffer_reserve);

	while (!fin.eof()){
		MeteoData md;
		StationData tmpsd = sd;
		double lat=IOUtils::nodata, lon=IOUtils::nodata, alt=IOUtils::nodata;
		unsigned int poscounter = 0;

		for (unsigned int ii=0; ii<nrOfColumns; ii++){
			if (vecDataSequence[ii] == "timestamp"){
				double tmpval;
				fin.read(reinterpret_cast < char * > (&tmpval), sizeof(double));

				md.date.setDate(tmpval, timezone);

				if (md.date < dateStart)
					continue;
				if (md.date > dateEnd)
					return;
			} else {
				float tmpval;
				fin.read(reinterpret_cast < char * > (&tmpval), sizeof(float));			
				double val = (double)tmpval;

				if (vecDataSequence[ii] == "latitude"){
					lat = val;
					poscounter++;
				} else if (vecDataSequence[ii] == "longitude"){
					lon = val;
					poscounter++;
				} else if (vecDataSequence[ii] == "altitude"){
					alt = val;
					poscounter++;
				} else {
					if(val==plugin_nodata)
						SMETIO::getParameter(vecDataSequence[ii], md) = IOUtils::nodata;
					else
						SMETIO::getParameter(vecDataSequence[ii], md) = (val * vecUnitsMultiplier[ii]) + vecUnitsOffset[ii];
				}
			}
		}

		char c;
		fin.read(&c, sizeof(char));
		if (c != '\n')
			throw InvalidFormatException("Corrupted data in section [DATA]", AT);

		if (poscounter == 3) {
			lat=IOUtils::standardizeNodata(lat, plugin_nodata);
			lon=IOUtils::standardizeNodata(lon, plugin_nodata);
			alt=IOUtils::standardizeNodata(alt, plugin_nodata);
			tmpsd.position.setLatLon(lat, lon, alt);
		}

		if (md.date >= dateStart){
			md.meta = tmpsd;
			vecMeteo.push_back(md);
		}
	}	
}

void SMETIO::readDataAscii(const char& eoln, const std::string& filename, const double& timezone,
                           const StationData& sd, const std::vector<std::string>& vecDataSequence,
                           const std::vector<double>& vecUnitsOffset, std::vector<double>& vecUnitsMultiplier,
                           const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo)
{
	string line = "";
	vector<string> tmpvec;
	const unsigned int nrOfColumns = vecDataSequence.size();

	tmpvec.reserve(nrOfColumns);
	vecMeteo.reserve(buffer_reserve);

	while (!fin.eof()){
		//HACK nodata mapping is NOT done!!!!!!
		//something like lat = IOUtils::standardizeNodata(lat, plugin_nodata); should be done
		getline(fin, line, eoln);
		IOUtils::stripComments(line);
		IOUtils::trim(line);
		if (line == "") continue; //Pure comment lines and empty lines are ignored

		unsigned int ncols = IOUtils::readLineToVec(line, tmpvec);

		if (ncols != nrOfColumns)
			throw InvalidFormatException("In "+ filename + ": Invalid amount of data in data line "+line, AT);

		MeteoData md;
		StationData tmpsd = sd;
		double lat=IOUtils::nodata, lon=IOUtils::nodata, alt=IOUtils::nodata;
		unsigned int poscounter = 0;
		for (unsigned int ii=0; ii<nrOfColumns; ii++){
			if (vecDataSequence[ii] == "timestamp"){
				if (!IOUtils::convertString(md.date, tmpvec[ii], timezone, std::dec))
					throw InvalidFormatException("In "+filename+": Timestamp "+tmpvec[ii]+" invalid in data line", AT);
				if (md.date < dateStart)
					continue;
				if (md.date > dateEnd)
					return;

			} else if (vecDataSequence[ii] == "latitude"){
				if (!IOUtils::convertString(lat, tmpvec[ii])) 
					throw InvalidFormatException("In "+filename+": Latitude invalid", AT);
				poscounter++;
			} else if (vecDataSequence[ii] == "longitude"){
				if (!IOUtils::convertString(lon, tmpvec[ii])) 
					throw InvalidFormatException("In "+filename+": Longitude invalid", AT);
				poscounter++;
			} else if (vecDataSequence[ii] == "altitude"){
				if (!IOUtils::convertString(alt, tmpvec[ii])) 
					throw InvalidFormatException("In "+filename+": Altitude invalid", AT);
				poscounter++;
			} else {
				double val;
				if (!IOUtils::convertString(val, tmpvec[ii]))
					throw InvalidFormatException("In "+filename+": Invalid value for param", AT);
				if(val==plugin_nodata)
					SMETIO::getParameter(vecDataSequence[ii], md) = IOUtils::nodata;
				else
					SMETIO::getParameter(vecDataSequence[ii], md) = (val * vecUnitsMultiplier[ii]) + vecUnitsOffset[ii];
			}
		}

		if (poscounter == 3) {
			lat=IOUtils::standardizeNodata(lat, plugin_nodata);
			lon=IOUtils::standardizeNodata(lon, plugin_nodata);
			alt=IOUtils::standardizeNodata(alt, plugin_nodata);
			tmpsd.position.setLatLon(lat, lon, alt);
		}

		if (md.date >= dateStart){
			md.meta = tmpsd;
			vecMeteo.push_back(md);
		}
	}
}

void SMETIO::readHeader(const char& eoln, const std::string& filename, bool& locationInHeader,
                        double& timezone, StationData& sd, std::vector<std::string>& vecDataSequence,
                        std::vector<double>& vecUnitsOffset, std::vector<double>& vecUnitsMultiplier)
{
	string line="";
	while (!fin.eof() && (fin.peek() != '['))
		getline(fin, line, eoln);
	
	getline(fin, line, eoln);
	IOUtils::stripComments(line);
	IOUtils::trim(line);
	IOUtils::toUpper(line);

	if (line != "[HEADER]")
		throw InvalidFormatException("Section " + line + " in "+ filename + " invalid", AT);

	std::map<std::string,std::string> mapHeader;
	while (!fin.eof() && (fin.peek() != '[')){
		getline(fin, line, eoln);

		IOUtils::stripComments(line);
		IOUtils::trim(line);

		if (line != "") {
			if (!IOUtils::readKeyValuePair(line, "=", mapHeader))
				throw InvalidFormatException("Invalid key value pair in section [Header]", AT);
		}
	}

	//Now extract info from mapHeader
	IOUtils::getValueForKey(mapHeader, "nodata", plugin_nodata);
	IOUtils::getValueForKey(mapHeader, "station_id", sd.stationID);
	IOUtils::getValueForKey(mapHeader, "station_name", sd.stationName, IOUtils::nothrow);
	timezone = in_dflt_TZ;
	IOUtils::getValueForKey(mapHeader, "tz", timezone, IOUtils::nothrow);
	if(timezone==plugin_nodata) //if a nodata was given in the header, we replace it with the io.ini timezone
		timezone = in_dflt_TZ;
	sd.position.setProj(coordin, coordinparam); //set the default projection from config file

	//trying to read easting/northing
	double easting=IOUtils::nodata, northing=IOUtils::nodata, alt=IOUtils::nodata;
	short int epsg=IOUtils::snodata;
	IOUtils::getValueForKey(mapHeader, "epsg", epsg, IOUtils::nothrow);
	if(epsg!=IOUtils::snodata) {
		sd.position.setEPSG(epsg);
	}
	
	IOUtils::getValueForKey(mapHeader, "easting", easting, IOUtils::nothrow);
	easting=IOUtils::standardizeNodata(easting, plugin_nodata);
	if (easting != IOUtils::nodata){
		IOUtils::getValueForKey(mapHeader, "northing", northing);
		IOUtils::getValueForKey(mapHeader, "altitude", alt);
		northing=IOUtils::standardizeNodata(northing, plugin_nodata);
		alt=IOUtils::standardizeNodata(alt, plugin_nodata);
		sd.position.setXY(easting, northing, alt);
		locationInHeader = true;
	} else {
		locationInHeader = false;
	}

	//now trying to read lat/long (second, so that it has precedence over east/north coordinates)
	double lat=IOUtils::nodata, lon=IOUtils::nodata;
	IOUtils::getValueForKey(mapHeader, "latitude", lat, IOUtils::nothrow);
	lat=IOUtils::standardizeNodata(lat, plugin_nodata);
	if (lat != IOUtils::nodata){ 
		IOUtils::getValueForKey(mapHeader, "longitude", lon);
		IOUtils::getValueForKey(mapHeader, "altitude", alt);
		lon=IOUtils::standardizeNodata(lon, plugin_nodata);
		alt=IOUtils::standardizeNodata(alt, plugin_nodata);
		sd.position.setLatLon(lat, lon, alt);
		locationInHeader = true;
	}

	IOUtils::getValueForKey(mapHeader, "fields", vecDataSequence);

	//trying to read units offsets
	std::vector<std::string> vecTmp;
	IOUtils::getValueForKey(mapHeader, "units_offset", vecTmp, IOUtils::nothrow); //TODO: implement these!!
	vecUnitsOffset.clear();
	if(vecTmp.size()==0) {
		for(unsigned int i=0; i<vecDataSequence.size(); i++)
			vecUnitsOffset.push_back(0.);
	} else if(vecTmp.size()==vecDataSequence.size()) {
		for(unsigned int i=0; i<vecDataSequence.size(); i++) {
			double tmp;
			IOUtils::convertString(tmp, vecTmp[i]);
			vecUnitsOffset.push_back(tmp);
		}
	} else {
		throw InvalidFormatException("header lines \"units_offset\" and \"fields\" do not match in "+ filename, AT);
	}

	//trying to read units multipliers
	vecTmp.clear();
	IOUtils::getValueForKey(mapHeader, "units_multiplier", vecTmp, IOUtils::nothrow);
	vecUnitsMultiplier.clear();
	if(vecTmp.size()==0) {
		for(unsigned int i=0; i<vecDataSequence.size(); i++)
			vecUnitsMultiplier.push_back(1.);
	} else if(vecTmp.size()==vecDataSequence.size()) {
		for(unsigned int i=0; i<vecDataSequence.size(); i++) {
			double tmp;
			IOUtils::convertString(tmp, vecTmp[i]);
			vecUnitsMultiplier.push_back(tmp);
		}
	} else {
		throw InvalidFormatException("header lines \"units_multipliers\" and \"fields\" do not match in "+ filename, AT);
	}

	//Read [DATA] section tag
	getline(fin, line, eoln);
	IOUtils::stripComments(line);
	IOUtils::trim(line);
	IOUtils::toUpper(line);

	if (line != "[DATA]")
		throw InvalidFormatException("Section " + line + " in "+ filename + " invalid, expected [DATA]", AT);
}

void SMETIO::checkSignature(const std::vector<std::string>& vecSignature, const std::string& filename, bool& isAscii)
{
	if ((vecSignature.size() != 3) || (vecSignature[0] != "SMET"))
		throw InvalidFormatException("The signature of file " + filename + " is invalid", AT);

	std::string version = vecSignature[1];
	if ((version != "0.9") && (version != "0.95") && (version != "0.99") && (version != "1.0") && (version != smet_version))
		throw InvalidFormatException("Unsupported file format version for file " + filename, AT);

	if(version=="0.9" || version=="0.95" || version=="0.99" || version=="1.0") {
		std::cout << "[W] SMET specification 1.1 changes the priorities of units_multiplier and units_offset. Please check/update your files and bring them to 1.1!!\n";
	}

	const std::string type = vecSignature[2];
	if (type == "ASCII")
		isAscii = true;
	else if (type == "BINARY")
		isAscii = false;
	else 
		throw InvalidFormatException("The 3rd column in the signature of file " + filename + " must be either ASCII or BINARY", AT);
}

void SMETIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	//Loop through all stations
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		//1. check consitency of station data position -> write location in header or data section
		StationData sd;
		sd.position.setProj(coordout, coordoutparam);
		bool isConsistent = checkConsistency(vecMeteo.at(ii), sd);
		
		if (sd.stationID == ""){
			stringstream ss;
			ss << "Station" << ii+1;

			sd.stationID = ss.str();
		}

		string filename = outpath + "/" + sd.stationID + ".smet";
		if (!IOUtils::validFileName(filename)) //Check whether filename is valid
			throw InvalidFileNameException(filename, AT);

		try {
			fout.open(filename.c_str());
			if (fout.fail()) throw FileAccessException(filename.c_str(), AT);

			fout << "SMET " << smet_version << " ";
			if (outputIsAscii)
				fout << "ASCII" << endl;
			else 
				fout << "BINARY" << endl;

			//2. write header, but first check which meteo parameter fields are actually in use
			vector<bool> vecParamInUse = vector<bool>(MeteoData::nrOfParameters, false); 
			double timezone = IOUtils::nodata;
			checkForUsedParameters(vecMeteo[ii], timezone, vecParamInUse);
			if(out_dflt_TZ!=IOUtils::nodata) {
				writeHeaderSection(isConsistent, sd, out_dflt_TZ, vecParamInUse);
			} else {
				writeHeaderSection(isConsistent, sd, timezone, vecParamInUse);
			}
			
			//3. write data depending on ASCII/BINARY
			if (outputIsAscii) {
				writeDataAscii(isConsistent, vecMeteo[ii], vecParamInUse);
			} else {
				fout << "[DATA]" << endl;
				fout.close();
				fout.open(filename.c_str(), ios::out | ios::app | ios::binary); //reopen as binary file
				if (fout.fail()) throw FileAccessException(filename.c_str(), AT);
				writeDataBinary(isConsistent, vecMeteo[ii], vecParamInUse);
			}
			
			cleanup();

			//4. gzip file or not
		} catch (exception& e){
			cleanup();
			throw;
		}
	}
}

void SMETIO::writeDataBinary(const bool& writeLocationInHeader, const std::vector<MeteoData>& vecMeteo,
                             const std::vector<bool>& vecParamInUse)
{
	char eoln = '\n';

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		float val = 0;

		double julian;
		if(out_dflt_TZ!=IOUtils::nodata) {
			Date tmp_date(vecMeteo[ii].date);
			tmp_date.setTimeZone(out_dflt_TZ);
			julian = tmp_date.getJulianDate();
		} else {
			julian = vecMeteo[ii].date.getJulianDate();
		}		
		fout.write((char*)&julian, sizeof(double));

		if (!writeLocationInHeader){ //Meta data changes
			val = (float)vecMeteo[ii].meta.position.getLat();
			fout.write((char*)&val, sizeof(float));
			val = (float)vecMeteo[ii].meta.position.getLon();
			fout.write((char*)&val, sizeof(float));
			val = (float)vecMeteo[ii].meta.position.getAltitude();
			fout.write((char*)&val, sizeof(float));
		}

		for (unsigned int jj=0; jj<MeteoData::nrOfParameters; jj++){
			if (vecParamInUse[jj]){
				val = (float)vecMeteo[ii].param(jj);
				fout.write((char*)&val, sizeof(float));
			}
		}
		fout.write((char*)&eoln, sizeof(char));
	}
}

void SMETIO::writeDataAscii(const bool& writeLocationInHeader, const std::vector<MeteoData>& vecMeteo,
                            const std::vector<bool>& vecParamInUse)
{
	fout << "[DATA]" << endl;
	fout.fill(' ');
	fout << right;
	fout << fixed;
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if(out_dflt_TZ!=IOUtils::nodata) {
			Date tmp_date(vecMeteo[ii].date);
			tmp_date.setTimeZone(out_dflt_TZ);
			fout << tmp_date.toString(Date::ISO);
		} else {
			fout << vecMeteo[ii].date.toString(Date::ISO);
		}

		if (!writeLocationInHeader){ //Meta data changes
			fout << " " << setw(12) << setprecision(6) << vecMeteo[ii].meta.position.getLat();
			fout << " " << setw(12) << setprecision(6) << vecMeteo[ii].meta.position.getLon();
			fout << " " << setw(8)  << setprecision(2) << vecMeteo[ii].meta.position.getAltitude();
		}

		for (unsigned int jj=0; jj<MeteoData::nrOfParameters; jj++){
			fout << " ";
			if (vecParamInUse[jj]){
				setFormatting(MeteoData::Parameters(jj));
				if (vecMeteo[ii].param(jj) == IOUtils::nodata)
					fout << setprecision(0);
				fout << vecMeteo[ii].param(jj);
			}
		}
		fout << endl;
	}
}

void SMETIO::setFormatting(const MeteoData::Parameters& paramindex)
{
	if ((paramindex == MeteoData::TA) || (paramindex == MeteoData::TSS) || (paramindex == MeteoData::TSG))
		fout << setw(8) << setprecision(2);
	else if ((paramindex == MeteoData::VW) || (paramindex == MeteoData::VW_MAX))
		fout << setw(6) << setprecision(1);
	else if (paramindex == MeteoData::DW)
		fout << setw(5) << setprecision(0);
	else if ((paramindex == MeteoData::ISWR) || (paramindex == MeteoData::RSWR) || (paramindex == MeteoData::ILWR))
		fout << setw(6) << setprecision(0);
	else if (paramindex == MeteoData::HNW)
		fout << setw(6) << setprecision(3);
	else if (paramindex == MeteoData::HS)
		fout << setw(8) << setprecision(3);
	else if (paramindex == MeteoData::RH)
		fout << setw(7) << setprecision(3);
}

void SMETIO::writeHeaderSection(const bool& writeLocationInHeader, const StationData& sd,
                                const double& timezone, const std::vector<bool>& vecParamInUse)
{
	fout << "[HEADER]" << endl;
	fout << "station_id   = " << sd.getStationID() << endl;
	if (sd.getStationName() != "")
		fout << "station_name = " << sd.getStationName() << endl;

	fout << fixed;
	if (writeLocationInHeader){ //TODO: only write if != nodata
		fout << "latitude     = " << setw(14) << setprecision(6) << sd.position.getLat() << "\n";
		fout << "longitude    = " << setw(14) << setprecision(6) << sd.position.getLon() << "\n";
		fout << "altitude     = " << setw(9)  << setprecision(1) << sd.position.getAltitude() << "\n";
		fout << "easting      = " << setw(14) << setprecision(6) << sd.position.getEasting() << "\n";
		fout << "northing     = " << setw(14) << setprecision(6) << sd.position.getNorthing() << "\n";
		fout << "epsg         = " << setw(7)  << setprecision(0) << sd.position.getEPSG() << "\n";
	}

	fout << "nodata       = " << setw(7) << setprecision(0) << IOUtils::nodata << "\n";
	
	if ((timezone != IOUtils::nodata) && (timezone != 0.0))
		fout << "tz           = " << setw(7)  << setprecision(0) << timezone << "\n";

	fout << "fields       = timestamp";

	if (!writeLocationInHeader){
		fout << " latitude longitude altitude";
	}

	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){
		if (vecParamInUse[ii]) {
			std::string column=MeteoData::getParameterName(ii);
			if(column=="RSWR") column="OSWR";
			if(column=="HNW") column="PSUM";
			fout << " " << column;
		}
	}
	fout << endl;
}


void SMETIO::checkForUsedParameters(const std::vector<MeteoData>& vecMeteo, double& timezone,
                                    std::vector<bool>& vecParamInUse)
{
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		for (unsigned int jj=0; jj<MeteoData::nrOfParameters; jj++){
			if (!vecParamInUse[jj])
				if (vecMeteo[ii].param(jj) != IOUtils::nodata)
					vecParamInUse[jj] = true;
		}
	}	

	if (vecMeteo.size() > 0)
		timezone = vecMeteo[0].date.getTimeZone();
}

bool SMETIO::checkConsistency(const std::vector<MeteoData>& vecMeteo, StationData& sd)
{
	if (vecMeteo.size() > 0) //HACK to get the station data even when encoutering bug 87
		sd = vecMeteo[0].meta;

	for (unsigned int ii=1; ii<vecMeteo.size(); ii++){
		if (vecMeteo[ii].meta.position != vecMeteo[ii-1].meta.position)
			return false; //BUG: if one is nodata -> we consider that positions are not consistent
	}

	return true;
}

void SMETIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SMETIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}


#ifndef _METEOIO_JNI
extern "C"
{
#define COMPILE_PLUGIN
#include "exports.h"

	METEOIO_EXPORT void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	METEOIO_EXPORT void* loadObject(const string& classname, const Config& cfg) {
		if(classname == "SMETIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new SMETIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
