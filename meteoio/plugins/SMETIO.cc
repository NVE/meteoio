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
 * - METEOFILE#: input filename and path. As many meteofiles as needed may be specified
 * - METEOPATH: output directory where to write the output meteofiles
 * - METEOPARAM: output file format options (ASCII or BINARY that might be followed by GZIP)
 *
 * Example:
 * @code
 * [Input]
 * METEOFILE1 = ./input/uppper_station.smet
 * METEOFILE2 = ./input/lower_station.smet
 * METEOFILE3 = ./input/outlet_station.smet
 * [Output]
 * METEOPATH = ./output
 * METEOPARAM = ASCII GZIP
 * @endcode
 */

const std::string SMETIO::smet_version = "0.99";
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
}

SMETIO::SMETIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	parseInputOutputSection();
}

SMETIO::SMETIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	parseInputOutputSection();
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
			readHeader(eoln, filename, locationInHeader, timezone, sd, vecDataSequence);
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
	in_dflt_TZ = out_dflt_TZ = 0.;
	cfg.getValue("TZ","Input",in_dflt_TZ,Config::nothrow);
	cfg.getValue("TZ","Output",out_dflt_TZ,Config::nothrow);

	// Parse the [Input] and [Output] sections within Config object cfg
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

	//Parse input section: extract number of files to read and store filenames in vecFiles
	unsigned int counter = 1;
	string filename = "";
	do {
		stringstream ss;
		filename = "";
		
		ss << "METEOFILE" << counter;
		cfg.getValue(ss.str(), "Input", filename, Config::nothrow);

		if (filename != ""){
			if (!IOUtils::validFileName(filename)) //Check whether filename is valid
				throw InvalidFileNameException(filename, AT);

			vecFiles.push_back(filename);
		}
		counter++;
	} while (filename != "");

	nr_stations = counter - 1;

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
                            std::vector< std::vector<StationData> >& vecStation,
                            const unsigned int& stationindex)
{
	//Make sure that vecMeteo/vecStation have the correct dimension and stationindex is valid
	unsigned int startindex=0, endindex=vecFiles.size();	
	if (stationindex != IOUtils::npos){
		if ((stationindex < vecFiles.size()) || (stationindex < vecMeteo.size()) || (stationindex < vecStation.size())){
			startindex = stationindex;
			endindex = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("Invalid stationindex", AT);
		} 

		vecMeteo[stationindex].clear();
		vecStation[stationindex].clear();
	} else {
		vecMeteo.clear();
		vecStation.clear();
		
		vecMeteo = vector< vector<MeteoData> >(vecFiles.size());
		vecStation = vector< vector<StationData> >(vecFiles.size());
		vecMeteo.reserve(nr_stations);
		vecStation.reserve(nr_stations);
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
			readHeader(eoln, filename, locationInHeader, timezone, sd, vecDataSequence);
			SMETIO::checkColumnNames(vecDataSequence, locationInHeader);

			//3. Read DATA
			if (isAscii){
				readDataAscii(eoln, filename, timezone, sd, vecDataSequence, dateStart, dateEnd,
						    vecMeteo[ii], vecStation[ii]);
			} else {
				streampos currpos = fin.tellg();
				fin.close();
				fin.open (filename.c_str(), ios::in|ios::binary);				
				if (fin.fail()) throw FileAccessException(filename, AT);
				fin.seekg(currpos); //jump to binary data section

				readDataBinary(eoln, filename, timezone, sd, vecDataSequence, dateStart, dateEnd,
							vecMeteo[ii], vecStation[ii]);
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
                            const Date& dateStart, const Date& dateEnd,
                            std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation)
{
	const unsigned int nrOfColumns = vecDataSequence.size();

	vecMeteo.reserve(buffer_reserve);
	vecStation.reserve(buffer_reserve);

	while (!fin.eof()){
		MeteoData md;
		if ((timezone != IOUtils::nodata) && (timezone != 0.0))
			md.date.setTimeZone(timezone);
		StationData tmpsd = sd;
		double lat=IOUtils::nodata, lon=IOUtils::nodata, alt=IOUtils::nodata;
		unsigned int poscounter = 0;

		for (unsigned int ii=0; ii<nrOfColumns; ii++){
			if (vecDataSequence[ii] == "timestamp"){
				double tmpval;
				fin.read(reinterpret_cast < char * > (&tmpval), sizeof(double));

				md.date.setDate(tmpval);

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
					SMETIO::getParameter(vecDataSequence[ii], md) = val;
				}
			}
		}

		char c;
		fin.read(&c, sizeof(char));
		if (c != '\n')
			throw InvalidFormatException("Corrupted data in section [DATA]", AT);

		if (poscounter == 3)
			tmpsd.position.setLatLon(lat, lon, alt);

		//cout << "===" << endl;
		//cout << sd << endl;
		//cout << md.date.toString(Date::ISO) << endl;
		if (md.date >= dateStart){
			vecMeteo.push_back(md);
			vecStation.push_back(tmpsd);
		}
	}	
}

void SMETIO::readDataAscii(const char& eoln, const std::string& filename, const double& timezone,
                           const StationData& sd, const std::vector<std::string>& vecDataSequence,
                           const Date& dateStart, const Date& dateEnd,
                           std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation)
{
	string line = "";
	vector<string> tmpvec;
	const unsigned int nrOfColumns = vecDataSequence.size();

	tmpvec.reserve(nrOfColumns);
	vecMeteo.reserve(buffer_reserve);
	vecStation.reserve(buffer_reserve);

	while (!fin.eof()){
		getline(fin, line, eoln);
		IOUtils::stripComments(line);
		IOUtils::trim(line);
		if (line == "") continue; //Pure comment lines and empty lines are ignored

		unsigned int ncols = IOUtils::readLineToVec(line, tmpvec);

		if (ncols != nrOfColumns)
			throw InvalidFormatException("In "+ filename + ": Invalid amount of data in data line "+line, AT);

		MeteoData md;
		if ((timezone != IOUtils::nodata) && (timezone != 0.0))
			md.date.setTimeZone(timezone);
		StationData tmpsd = sd;
		double lat, lon, alt;
		unsigned int poscounter = 0;
		for (unsigned int ii=0; ii<nrOfColumns; ii++){
			if (vecDataSequence[ii] == "timestamp"){
				if (!IOUtils::convertString(md.date, tmpvec[ii]))
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
				if (!IOUtils::convertString(SMETIO::getParameter(vecDataSequence[ii], md), tmpvec[ii]))
					throw InvalidFormatException("In "+filename+": Invalid value for param", AT);
			}
		}

		if (poscounter == 3)
			tmpsd.position.setLatLon(lat, lon, alt);

		if (md.date >= dateStart){
			vecMeteo.push_back(md);
			vecStation.push_back(tmpsd);
		}
	}
}

void SMETIO::readHeader(const char& eoln, const std::string& filename, bool& locationInHeader,
                         double& timezone, StationData& sd, std::vector<std::string>& vecDataSequence)
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
	IOUtils::getValueForKey(mapHeader, "station_id", sd.stationID);
	IOUtils::getValueForKey(mapHeader, "station_name", sd.stationName, IOUtils::nothrow);
	timezone = in_dflt_TZ;
	IOUtils::getValueForKey(mapHeader, "tz", timezone, IOUtils::nothrow);
	sd.position.setProj(coordin, coordinparam); //set the default projection from config file

	//trying to read easting/northing
	double easting=IOUtils::nodata, northing=IOUtils::nodata, alt=IOUtils::nodata;
	short int epsg=IOUtils::snodata;
	IOUtils::getValueForKey(mapHeader, "epsg", epsg, IOUtils::nothrow);
	if(epsg!=IOUtils::snodata) {
		sd.position.setEPSG(epsg);
	}
	
	IOUtils::getValueForKey(mapHeader, "easting", easting, IOUtils::nothrow);
	if (easting != IOUtils::nodata){
		IOUtils::getValueForKey(mapHeader, "northing", northing);
		IOUtils::getValueForKey(mapHeader, "altitude", alt);
		sd.position.setXY(easting, northing, alt);
		locationInHeader = true;
	} else {
		locationInHeader = false;
	}

	//now trying to read lat/long (second, so that it has precedence over east/north coordinates)
	double lat=IOUtils::nodata, lon=IOUtils::nodata;
	IOUtils::getValueForKey(mapHeader, "latitude", lat, IOUtils::nothrow);
	if (lat != IOUtils::nodata){ 
		IOUtils::getValueForKey(mapHeader, "longitude", lon);
		IOUtils::getValueForKey(mapHeader, "altitude", alt);
		sd.position.setLatLon(lat, lon, alt);
		locationInHeader = true;
	}

	IOUtils::getValueForKey(mapHeader, "fields", vecDataSequence);
	//IOUtils::getValueForKey(mapHeader, "units_offset", vecUnitsOffset, IOUtils::nothrow); //TODO: implement these!!
	//IOUtils::getValueForKey(mapHeader, "units_multiplier", vecUnitsMultiplier, IOUtils::nothrow);

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
	if ((version != "0.9") && (version != "0.95") && (version != smet_version))
		throw InvalidFormatException("Unsupported file format version for file " + filename, AT);

	const std::string type = vecSignature[2];
	if (type == "ASCII")
		isAscii = true;
	else if (type == "BINARY")
		isAscii = false;
	else 
		throw InvalidFormatException("The 3rd column in the signature of file " + filename + " must be either ASCII or BINARY", AT);
}

void SMETIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
					    const std::vector< std::vector<StationData> >& vecStation,
					    const std::string&)
{

	//Loop through all stations
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		//1. check consitency of station data position -> write location in header or data section
		StationData sd;
		sd.position.setProj(coordout, coordoutparam);
		bool isConsistent = checkConsistency(vecStation.at(ii), sd);
		
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
			writeHeaderSection(isConsistent, sd, timezone, vecParamInUse);
			
			//3. write data depending on ASCII/BINARY
			if (outputIsAscii) {
				writeDataAscii(isConsistent, vecMeteo[ii], vecStation[ii], vecParamInUse);
			} else {
				fout << "[DATA]" << endl;
				fout.close();
				fout.open(filename.c_str(), ios::out | ios::app | ios::binary); //reopen as binary file
				if (fout.fail()) throw FileAccessException(filename.c_str(), AT);
				writeDataBinary(isConsistent, vecMeteo[ii], vecStation[ii], vecParamInUse);
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
                              const std::vector<StationData>& vecStation, const std::vector<bool>& vecParamInUse)
{
	char eoln = '\n';

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		float val = 0;
		double julian = vecMeteo[ii].date.getJulianDate();
		
		fout.write((char*)&julian, sizeof(double));

		if (!writeLocationInHeader){ //Meta data changes
			val = (float)vecStation[ii].position.getLat();
			fout.write((char*)&val, sizeof(float));
			val = (float)vecStation[ii].position.getLon();
			fout.write((char*)&val, sizeof(float));
			val = (float)vecStation[ii].position.getAltitude();
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
                             const std::vector<StationData>& vecStation, const std::vector<bool>& vecParamInUse)
{
	fout << "[DATA]" << endl;
	fout.fill(' ');
	fout << right;
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		fout << vecMeteo[ii].date.toString(Date::ISO);

		if (!writeLocationInHeader){ //Meta data changes
			fout << " " << setw(12) << setprecision(6) << vecStation[ii].position.getLat();
			fout << " " << setw(12) << setprecision(6) << vecStation[ii].position.getLon();
			fout << " " << setw(8)  << setprecision(2) << vecStation[ii].position.getAltitude();
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
	else if (paramindex == MeteoData::VW)
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

	if (writeLocationInHeader){ //TODO: only write if != nodata
		fout << fixed;
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

bool SMETIO::checkConsistency(const std::vector<StationData>& vecStation, StationData& sd)
{
	for (unsigned int ii=1; ii<vecStation.size(); ii++){
		if (vecStation[ii].position != vecStation[ii-1].position)
			return false;
	}

	if (vecStation.size() > 0)
		sd = vecStation[0];

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
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const string& classname, const Config& cfg) {
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
