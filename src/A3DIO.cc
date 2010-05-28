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
#include "A3DIO.h"

using namespace std;

namespace mio {
/**
 * @page a3d A3D
 * @section a3d_format Format
 * This plugin reads legacy Alpine3D meteorological input files. It reads the meteo1d.txt file that contains, as measured at one unique location, the following fields:
 * - air temperature in field ta
 * - incoming short wave radiation in field iswr
 * - wind velocity in field vw
 * - relative humidity in field rh
 * - incoming long wave radiation in field ea
 * - precipitations in field nswc
 *
 * and optionally a list of stations with their measurements (meteo2d files) for the following parameters (with YYYY being the 4-digits year):
 * - precipitations (file named precYYYY.txt)
 * - relative humidity (file named rhumYYYY.txt)
 * - air temperature (file named tairYYYY.txt)
 * - wind speed (file named wspdYYYY.txt)
 * - and optionnally wind direction (file named wdirYYYY.txt)
 *
 * @section a3d_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 *
 * @section a3d_keywords Keywords
 * This plugin uses the following keywords:
 * - METEOPATH: string containing the path to the meteorological files (ie: where to find meteo1d.txt and meteo2d files)
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - SPECIALPTSFILE: a path+file name to the a file containing grid coordinates of special points of interest (for special outputs)
 */

const double A3DIO::plugin_nodata = -9999.0; //plugin specific nodata value

//A3DIO::A3DIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename)
// {
// 	getProjectionParameters();
// }

//Main constructor
A3DIO::A3DIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

//Copy constructor
A3DIO::A3DIO(const A3DIO& aio) : IOInterface(NULL), cfg(aio.cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

A3DIO::A3DIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

A3DIO::~A3DIO() throw()
{
	cleanup();
}

//Clone function
//A3DIO* A3DIO::clone() const { return new A3DIO(*this); }

void A3DIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
}

void A3DIO::read2DGrid(Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void A3DIO::readDEM(DEMObject& dem_out)
{
	//Nothing so far
	(void)dem_out;
	throw IOException("Nothing implemented here", AT);
}

void A3DIO::readLanduse(Grid2DObject& landuse_out)
{
	//Nothing so far
	(void)landuse_out;
	throw IOException("Nothing implemented here", AT);
}

void A3DIO::readAssimilationData(const Date& date_in, Grid2DObject& da_out)
{
	//Nothing so far
	(void)date_in;
	(void)da_out;
	throw IOException("Nothing implemented here", AT);
}

void A3DIO::write2DGrid(const Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void A3DIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& data,
					  const std::vector< std::vector<StationData> >&vec_stations,
					  const std::string&)
{
	if(vec_stations.size()>0) {
		//A3D format does not support stations changing position over time
		create1DFile(data, vec_stations);
		write2DMeteo(data, vec_stations);
	}
}

void A3DIO::readStationData(const Date&, std::vector<StationData>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void A3DIO::readMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	std::vector< std::vector<StationData> > vecStation;
	readMeteoData(dateStart, dateEnd, vecMeteo, vecStation);
}

void A3DIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
					 std::vector< std::vector<MeteoData> >& vecMeteo,
					 std::vector< std::vector<StationData> >& vecStation,
					 const unsigned int& stationindex)
{
	//if dateStart and dateEnd are the same: return exact match for date
	//if dateStart > dateEnd: return first data set with date > dateStart
	//read in all data starting with dateStart until dateEnd
	//if there is no data at all then the vector will be empty, no exception will be thrown
	(void)stationindex;
	vecMeteo.clear();
	vecStation.clear();

	//first from meteo1d.txt
	read1DMeteo(dateStart, dateEnd, vecMeteo, vecStation);

	//then all corresponding data sets from the 2d meteo files
	//note: they have to be entirely corresponding (for every date)
	try {
		read2DMeteo(vecMeteo, vecStation);
	} catch(std::exception& e){
		std::cerr << "[E] No meteo2d data found or error while reading it, using only Meteo1D data: "
			<< std::endl << "\t" << e.what() << std::endl;
	}
}

void A3DIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts RH to [0,1]
	if(meteo.ta!=IOUtils::nodata) {
		meteo.ta=C_TO_K(meteo.ta);
	}
	if(meteo.tsg!=IOUtils::nodata) {
		meteo.tsg=C_TO_K(meteo.tsg);
	}
	if(meteo.rh!=IOUtils::nodata) {
		meteo.rh /= 100.;
	}
}

void A3DIO::read1DMeteo(const Date& dateStart, const Date& dateEnd,
				std::vector< std::vector<MeteoData> >& vecMeteo,
				std::vector< std::vector<StationData> >& vecStation)
{
	double latitude=IOUtils::nodata, longitude=IOUtils::nodata,
		xcoord=IOUtils::nodata, ycoord=IOUtils::nodata, altitude=IOUtils::nodata;
	std::string tmp="", line="";
	Date tmp_date;
	std::vector<std::string> tmpvec;
	std::map<std::string, std::string> header; // A map to save key value pairs of the file header
	MeteoData tmpdata;
	StationData sd;
	bool eofreached = false;

	cfg.getValue("METEOPATH", "Input", tmp);
	tmp += "/meteo1d.txt";

	if (!IOUtils::fileExists(tmp)) {
		throw FileNotFoundException(tmp, AT);
	}

	fin.clear();
	fin.open (tmp.c_str(), std::ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(tmp,AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	//Go through file, save key value pairs
	try {
		//Read in station meta data
		IOUtils::readKeyValueHeader(header, fin, 5, "="); //Read in 5 lines as header
		IOUtils::getValueForKey(header, "Latitude", latitude);
		IOUtils::getValueForKey(header, "Longitude", longitude);
		IOUtils::getValueForKey(header, "X_Coord", xcoord);
		IOUtils::getValueForKey(header, "Y_Coord", ycoord);
		IOUtils::getValueForKey(header, "Altitude", altitude);

		//HACK!! would it be possible for getValueForKey() to do this transparently? (with a user flag)
		latitude = IOUtils::standardizeNodata(latitude, plugin_nodata);
		longitude = IOUtils::standardizeNodata(longitude, plugin_nodata);
		altitude = IOUtils::standardizeNodata(altitude, plugin_nodata);
		xcoord = IOUtils::standardizeNodata(xcoord, plugin_nodata);
		ycoord = IOUtils::standardizeNodata(ycoord, plugin_nodata);

		//compute/check WGS coordinates (considered as the true reference) according to the projection as defined in cfg
		Coords location(coordin, coordinparam);
		location.setXY(xcoord, ycoord, altitude, false);
		location.setLatLon(latitude, longitude, altitude, false);
		try {
			location.check();
		} catch(...) {
			std::cerr << "[E] Error in geographic coordinates in file " << tmp << " trapped at " << AT << std::endl;
			throw;
		}

		sd.setStationData(location, "");

		//Read one line, construct Date object and see whether date is greater or equal than the date_in object
		IOUtils::skipLines(fin, 1, eoln); //skip rest of line

		//Loop going through the data sequentially until dateStart is found
		do {
			getline(fin, line, eoln); //read complete line
			eofreached = readMeteoDataLine(line, tmpdata, tmp);
			//tmpdata.cleanData();
			convertUnits(tmpdata);

		} while((tmpdata.date < dateStart) && (!eofreached));

		if ((dateEnd < dateStart) && (!eofreached)){ //Special case
			vecMeteo.push_back( std::vector<MeteoData>() );
			vecStation.push_back( std::vector<StationData>() );
			vecMeteo[0].push_back(tmpdata);
			vecStation[0].push_back(sd);
		} else if ((tmpdata.date <= dateEnd)  && (!eofreached)) {
			vecMeteo.push_back( std::vector<MeteoData>() );
			vecStation.push_back( std::vector<StationData>() );
		}

		while ((tmpdata.date <= dateEnd)  && (!eofreached)) {
			//At this point tmpdata.date is >= dateStart
			vecMeteo[0].push_back(tmpdata);
			vecStation[0].push_back(sd);

			getline(fin, line, eoln); //read complete line
			eofreached = readMeteoDataLine(line, tmpdata, tmp);
			//tmpdata.cleanData();
			convertUnits(tmpdata);
		}
		//cout << "Size of buffer: " << vecMeteo[0].size() << "   " << tmp_date << endl;
	} catch(...) {
		std::cout << "[E] " << AT << ": "<< std::endl;
		cleanup();
		throw;
	}

	cleanup();
}

bool A3DIO::readMeteoDataLine(std::string& line, MeteoData& tmpdata, std::string filename)
{
	Date tmp_date;
	int tmp_ymdh[4];
	std::vector<std::string> tmpvec;
	double tmp_values[6];

	if (IOUtils::readLineToVec(line, tmpvec) != 10) {
		return true;
		//throw InvalidFormatException("Premature End of Line or no data for date " + date_in + " found in File " + filename, AT);
	}

	for (int ii=0; ii<4; ii++) {
		if (!IOUtils::convertString(tmp_ymdh[ii], tmpvec.at(ii), std::dec))
			throw InvalidFormatException(filename + ": " + line, AT);
	}

	tmp_date.setDate(tmp_ymdh[0],tmp_ymdh[1],tmp_ymdh[2],tmp_ymdh[3], 0);

	//Read rest of line with values ta, iswr, vw, rh, ea, hnw

	for (int ii=0; ii<6; ii++) { //go through the columns
		if (!IOUtils::convertString(tmp_values[ii], tmpvec.at(ii+4), std::dec)) {
			throw InvalidFormatException(filename + ": " + line, AT);
		}
	}

	tmpdata.setDate(tmp_date);
	tmpdata.setData(MeteoData::TA, tmp_values[0]);
	tmpdata.setData(MeteoData::ISWR, tmp_values[1]);
	tmpdata.setData(MeteoData::VW, tmp_values[2]);
	tmpdata.setData(MeteoData::RH, tmp_values[3]);
	tmpdata.setData(MeteoData::ILWR, tmp_values[4]);
	tmpdata.setData(MeteoData::HNW, tmp_values[5]);

	return false;
}


/*
  Preamble: Files are in METEOFILE directory. 4 types of files:
  prec????.txt == hnw
  rh????.txt == rh
  ta????.txt == ta
  wspd????.txt == vw

  Remarks: The headers of the files may defer - for each unique
  StationData one MeteoData and one StationData object will be created
*/
 void A3DIO::read2DMeteo(std::vector< std::vector<MeteoData> >& vecMeteo, std::vector< std::vector<StationData> >& vecStation)
{

	unsigned int stations=0, bufferindex=0;
	std::map<std::string, unsigned int> hashStations = std::map<std::string, unsigned int>();
	std::vector<std::string> filenames = std::vector<std::string>();

	//Requirement: meteo1D data must exist:
	if ((vecMeteo.size() == 0) || (vecMeteo[0].size() == 0))
		return;

	//1D and 2D data must correspond, that means that if there is 1D data
	//for a certain date (e.g. 1.1.2006) then 2D data must exist (prec2006.txt etc),
	//otherwise throw FileNotFoundException
	Date startDate(vecMeteo[0][0].date);
	Date endDate(vecMeteo[0][vecMeteo[0].size()-1].date);

	constructMeteo2DFilenames(startDate, endDate, filenames);//get all files for all years
	stations = getNrOfStations(filenames, hashStations);

	constructMeteo2DFilenames(startDate, startDate, filenames);//get filenames for current year
	std::cerr << "[I] Number of 2D meteo stations: " << stations << std::endl;

	if (stations < 1) {
		throw InvalidFormatException("[E] No StationData found in 2D Meteo Files", AT);
	}

	std::vector<StationData> tmpvecS = std::vector<StationData>(stations); //stores unique stations

	try {
		for (unsigned int ii=0; ii<filenames.size(); ii++){
			read2DMeteoHeader(filenames[ii], hashStations, tmpvecS);
		}

		//init vecStation with proper StationData, vecMeteo with nodata
		for (unsigned int jj=0; jj<tmpvecS.size(); jj++){
			vecMeteo.push_back( std::vector<MeteoData>() );
			vecStation.push_back( std::vector<StationData>() );
			for (unsigned int ii=0; ii<vecMeteo[0].size(); ii++){
				//NOTE: there needs to be the same amount of 1D and 2D data
				vecStation[jj+1].push_back(tmpvecS[jj]);
				vecMeteo[jj+1].push_back(MeteoData());
			}
		}

		do {
			unsigned int currentindex = bufferindex;
			read2DMeteoData(filenames[0], "nswc", hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[1], "rh", hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[2], "ta", hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[3], "vw", hashStations, vecMeteo, bufferindex);

			if(filenames.size() == 5) { //for keeping dw optional
				bufferindex = currentindex;
				read2DMeteoData(filenames.at(4), "dw", hashStations, vecMeteo, bufferindex);
			}
			//cerr << "bufferindex: " << bufferindex << "  Expected size()" << vecMeteo[0].size() << endl;

			if (bufferindex < (vecMeteo[0].size())) { //number of 1D meteo data
				//construct new filenames for the continued buffering
				constructMeteo2DFilenames(vecMeteo[0][bufferindex].date, vecMeteo[0][bufferindex].date, filenames);
			}
		} while(bufferindex < (vecMeteo[0].size()));
	} catch(...) {
		//clear all 2D meteo data if error occurs
		if (vecMeteo.size() > 1)
			vecMeteo.erase(vecMeteo.begin()+1, vecMeteo.end());

		if (vecStation.size() > 1)
			vecStation.erase(vecStation.begin()+1, vecStation.end());

		cleanup();
		throw;
	}

	//clean data and convert the units
	for (unsigned int ii=1; ii<vecMeteo.size(); ii++) { //loop over all stations except 1D Meteo
		for (unsigned int jj=0; jj<vecMeteo[ii].size(); jj++){ //Meteo1D data already cleaned
			//vecMeteo[ii][jj].cleanData();
			convertUnits(vecMeteo[ii][jj]);
		}
	}
}

void A3DIO::constructMeteo2DFilenames(const Date& startDate, const Date& endDate, std::vector<std::string>& filenames)
{
	int startyear=0, endyear=0, dummy=0;
	std::string tmp;

	filenames.clear();

	startDate.getDate(startyear, dummy, dummy, dummy, dummy);
	endDate.getDate(endyear, dummy, dummy, dummy, dummy);
	cfg.getValue("METEOPATH", "Input", tmp);

	for (int yyyy = startyear; yyyy<=endyear; yyyy++){
		std::stringstream ss;
		ss << yyyy;

		std::string precFilename = tmp + "/prec" + ss.str() + ".txt";
		std::string rhFilename = tmp + "/rhum" + ss.str() + ".txt";
		std::string taFilename = tmp + "/tair" + ss.str() + ".txt";
		std::string wspdFilename = tmp + "/wspd" + ss.str() + ".txt";
		std::string wdirFilename = tmp + "/wdir" + ss.str() + ".txt";

		filenames.push_back(precFilename);
		filenames.push_back(rhFilename);
		filenames.push_back(taFilename);
		filenames.push_back(wspdFilename);

		if (IOUtils::fileExists(wdirFilename)) //keeping wdir optional
			filenames.push_back(wdirFilename);
	}

	for (unsigned int ii=0; ii<filenames.size(); ii++) {
		if (!IOUtils::fileExists(filenames[ii])) {
			throw FileNotFoundException(filenames[ii], AT);
		}
	}
}


unsigned int A3DIO::getNrOfStations(std::vector<std::string>& filenames, std::map<std::string, unsigned int>& hashStations)
{
	std::vector<std::string> tmpvec;
	std::string line_in="";

	for (unsigned int ii=0; ii<filenames.size(); ii++) {
		//cout << *it << endl;
		std::string filename = filenames[ii];

		fin.clear();
		fin.open (filename.c_str(), std::ifstream::in);
		if (fin.fail()) throw FileAccessException(filename, AT);

		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		IOUtils::skipLines(fin, 4, eoln);
		getline(fin, line_in, eoln); //5th line holds the names of the stations
		unsigned int cols = IOUtils::readLineToVec(line_in, tmpvec);
		if ( cols > 4) { // if there are any stations
			//check each station name and whether it's already hashed, otherwise: hash!
			for (unsigned int ii=4; ii<cols; ii++) {
				unsigned int tmp_int = hashStations.count(tmpvec.at(ii));
				if (tmp_int == 0) {
					//cout << "Found station: " << tmpvec.at(ii) << endl;
					hashStations[tmpvec.at(ii)] = hashStations.size();
				}
			}
		}
		cleanup();
	}

	return (hashStations.size());
}

void A3DIO::read2DMeteoData(const std::string& filename, const std::string& parameter,
					std::map<std::string,unsigned int>& hashStations,
					std::vector< std::vector<MeteoData> >& vecM, unsigned int& bufferindex)
{

	std::string line_in = "";
	unsigned int columns;
	std::vector<std::string> tmpvec, vec_names;
	Date tmp_date;
	int tmp_ymdh[4];

	fin.clear();
	fin.open (filename.c_str(), std::ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	IOUtils::skipLines(fin, 4, eoln); //skip first 4 lines
	getline(fin, line_in, eoln); //line containing UNIQUE station names
	columns = IOUtils::readLineToVec(line_in, vec_names);
	if (columns < 4) {
		throw InvalidFormatException("[E] Premature end of line in file " + filename, AT);
	}

	MeteoData& lastMeteoData = vecM[0][vecM[0].size()-1]; //last time stamp in buffer of 1D meteo

	do {
		getline(fin, line_in, eoln);
		std::string tmpline = line_in;
		IOUtils::trim(tmpline);

		if (tmpline=="") {
			break;
		}

		if (IOUtils::readLineToVec(line_in, tmpvec)!=columns) { //Every station has to have its own column
			throw InvalidFormatException("[E] Premature End of Line or no data for date "
			                             + vecM[0][bufferindex].date.toString(Date::FULL) + " found in File "
			                             + filename, AT);
		}

		for (int ii=0; ii<4; ii++) {
			if (!IOUtils::convertString(tmp_ymdh[ii], tmpvec[ii], std::dec)) {
				throw InvalidFormatException("[E] Check date columns in " + filename, AT);
			}
		}
		tmp_date.setDate(tmp_ymdh[0],tmp_ymdh[1],tmp_ymdh[2],tmp_ymdh[3],0);

		MeteoData& currentMeteoData = vecM[0][bufferindex]; //1D Element to synchronize date
		if (tmp_date == currentMeteoData.date) {
			//Read in data
			for (unsigned int ii=4; ii<columns; ii++) {
				unsigned int stationnr = hashStations[vec_names.at(ii)];
				MeteoData& tmpmd = vecM[stationnr][bufferindex];
				tmpmd.date = tmp_date;

				if (parameter == "nswc") {
					if (!IOUtils::convertString(tmpmd.hnw, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For hnw value in " + filename + "  for date " + tmpmd.date.toString(Date::FULL), AT);
					}

				} else if (parameter == "rh") {
					if (!IOUtils::convertString(tmpmd.rh, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For rh value in " + filename + "  for date " + tmpmd.date.toString(Date::FULL), AT);
					}

				} else if (parameter == "ta") {
					if (!IOUtils::convertString(tmpmd.ta, tmpvec[ii], std::dec))
						throw ConversionFailedException("For ta value in " + filename + "  for date " + tmpmd.date.toString(Date::FULL), AT);

				} else if (parameter == "vw") {
					if (!IOUtils::convertString(tmpmd.vw, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For vw value in " + filename + "  for date " + tmpmd.date.toString(Date::FULL), AT);
					}
				} else if (parameter == "dw") {
					if (!IOUtils::convertString(tmpmd.dw, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For dw value in " + filename + "  for date " + tmpmd.date.toString(Date::FULL), AT);
					}
				}
			}

			bufferindex++;
		}
	} while((tmp_date<lastMeteoData.date) && (!fin.eof()));

	cleanup();
}

void A3DIO::read2DMeteoHeader(const std::string& filename, std::map<std::string,unsigned int>& hashStations,
				std::vector<StationData>& vecS)
{
	std::string line_in = "";
	unsigned int columns = 0;
	std::vector<std::string> vec_altitude, vec_xcoord, vec_ycoord, vec_names;

	fin.clear();
	fin.open (filename.c_str(), std::ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	IOUtils::skipLines(fin, 1, eoln);

	//Read all relevant lines in
	getline(fin, line_in, eoln); //Altitude
	columns = IOUtils::readLineToVec(line_in, vec_altitude);

	getline(fin, line_in, eoln); //xcoord
	if (IOUtils::readLineToVec(line_in, vec_xcoord) != columns) {
		throw InvalidFormatException("Column count doesn't match from line to line in " + filename, AT);
	}

	getline(fin, line_in, eoln); //ycoord
	if (IOUtils::readLineToVec(line_in, vec_ycoord) != columns) {
		throw InvalidFormatException("Column count doesn't match from line to line in " + filename, AT);
	}

	getline(fin, line_in, eoln); //names
	if (IOUtils::readLineToVec(line_in, vec_names) != columns) {
		throw InvalidFormatException("Column count doesn't match from line to line in " + filename, AT);
	}

	cleanup();

	//Check for duplicate station names within one file ... station names need to be unique!
	vector<string> vec_dup = vec_names;
	for (unsigned int ii=0; ii<vec_names.size(); ii++){
		const string& tmp = vec_names[ii];
		for (unsigned int jj=0; jj<vec_dup.size(); jj++){
			if (jj != ii){
				if (vec_dup[jj] == tmp)
					throw IOException("Duplicate station names detected in " + filename, AT);
			}
		}
	}

	//Build Coords object to convert easting/northing values to lat/long in WGS84
	Coords coordinate(coordin, coordinparam);

	for (unsigned int ii=4; ii<columns; ii++) {
		unsigned int stationnr = hashStations[vec_names.at(ii)];
		double altitude, easting, northing;
		std::string stationName;
		if ((!IOUtils::convertString(altitude, vec_altitude.at(ii), std::dec))
		    || (!IOUtils::convertString(easting, vec_xcoord.at(ii), std::dec))
		    || (!IOUtils::convertString(northing, vec_ycoord.at(ii), std::dec))
		    || (!IOUtils::convertString(stationName, vec_names.at(ii), std::dec))) {
			throw ConversionFailedException("Conversion of station description failed in " + filename, AT);
		}
		coordinate.setXY(easting, northing, altitude);
		vecS[stationnr-1].stationName = stationName;
		vecS[stationnr-1].position = coordinate;
	}
}

void A3DIO::readSpecialPoints(std::vector<Coords>& pts)
{
	std::string filename="", line_in="";
	std::vector<std::string> tmpvec;
	std::vector< std::pair<int,int> > mypts;

	cfg.getValue("SPECIALPTSFILE", "Input", filename); // cout << tmp << endl;
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}

	fin.clear();
	fin.open (filename.c_str(), std::ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename,AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	while (!fin.eof()) {
		getline(fin, line_in, eoln);

		if (IOUtils::readLineToVec(line_in, tmpvec)==2) { //Try to convert
			int x, y;
			if (!IOUtils::convertString(x, tmpvec.at(0), std::dec)) {
				throw ConversionFailedException("Conversion of a value failed in " + filename + " line: " + line_in, AT);
			}

			if (!IOUtils::convertString(y, tmpvec.at(1), std::dec)) {
				throw ConversionFailedException("Conversion of a value failed in " + filename + " line: " + line_in, AT);
			}

			std::pair<int,int> tmppair(x,y);
			mypts.push_back(tmppair);
		}
	}
	cleanup();

	//Now put everything into the output vector TODO: don't do any intermediate steps... copy directly into vector!
	Coords tmp_pts;
	for (unsigned int jj=0; jj<mypts.size(); jj++) {
		tmp_pts.setGridIndex(mypts.at(jj).first, mypts.at(jj).second, IOUtils::inodata, false);
		pts.push_back(tmp_pts);
	}
}

int A3DIO::create1DFile(const std::vector< std::vector<MeteoData> >& data, const std::vector< std::vector<StationData> >& stations)
{//TODO: add check for stations' positions
	std::string tmp_path;
	cfg.getValue("METEOPATH", "Output", tmp_path);
	const unsigned int sta_nr = stations.size();
	if(sta_nr==0) return EXIT_FAILURE;

	for(unsigned int ii=0; ii<sta_nr; ii++) {
		const unsigned int size = data[ii].size();
		if(size>0) {
			const std::string filename = tmp_path+"/meteo1D_"+stations[ii][0].stationName+".txt";
			std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
			if(!file) {
				throw FileAccessException("[E] Can not open file "+filename, AT);
			}

			file << "Name = " << stations[ii][0].stationName << "\n";
			file << "Latitude = " << stations[ii][0].position.getLat() << "\n";
			file << "Longitude = " << stations[ii][0].position.getLon() << "\n";
			file << "X_Coord = " << stations[ii][0].position.getEasting() << "\n";
			file << "Y_Coord = " << stations[ii][0].position.getNorthing() << "\n";
			file << "Altitude = " << stations[ii][0].position.getAltitude() << "\n";
			file << "YYYY MM DD HH ta iswr vw rh ea nswc\n";

			file.flags ( std::ios::fixed );
			for(unsigned int j=0; j<size; j++) {
				int yyyy, mm, dd, hh;
				data[ii][j].date.getDate(yyyy, mm, dd, hh);
				file.fill('0');
				file << setw(4) << yyyy << " " << setw(2) << mm << " " << setw(2) << dd << " " << setw(2) << hh << " ";
				file.fill(' ');
				if(data[ii][j].ta == IOUtils::nodata)
					file << setw(6) << setprecision(0) <<  IOUtils::nodata << " ";
				else
					file << setw(6) << setprecision(2) <<  K_TO_C(data[ii][j].ta) << " ";
				if(data[ii][j].iswr == IOUtils::nodata)
					file << setw(6) << setprecision(0) << IOUtils::nodata << " ";
				else
					file << setw(6) << setprecision(2) << data[ii][j].iswr << " ";
				if(data[ii][j].vw == IOUtils::nodata)
					file << setw(6) << setprecision(0) << IOUtils::nodata << " ";
				else
					file << setw(6) << setprecision(2) << data[ii][j].vw << " ";
				if(data[ii][j].rh == IOUtils::nodata)
					file << setw(6) << setprecision(0) << IOUtils::nodata << " ";
				else
					file << setw(6) << setprecision(2) << data[ii][j].rh * 100. << " ";
				if(data[ii][j].ilwr == IOUtils::nodata)
					file << setw(6) << setprecision(0) << IOUtils::nodata << " ";
				else
					file << setw(6) << setprecision(2) << data[ii][j].ilwr << " ";
				if(data[ii][j].hnw == IOUtils::nodata)
					file << setw(6) << setprecision(0) << IOUtils::nodata << "\n";
				else
					file << setw(6) << setprecision(2) << data[ii][j].hnw << "\n";
			}
			file.close();
		}
	}
	return EXIT_SUCCESS;
}

int A3DIO::writeHeader(std::ofstream &file, const std::vector< std::vector<StationData> >& stations, const std::string parameter_name)
{
	std::ostringstream str_altitudes;
	std::ostringstream str_eastings;
	std::ostringstream str_northings;
	const unsigned int sta_nr = stations.size();
	if(sta_nr==0) return EXIT_FAILURE;
	
	file << "X:\\filepath " << parameter_name <<endl;
	for(unsigned int ii=0;ii<sta_nr;ii++) {
		if(stations[ii].size() > 0) {
			str_altitudes << stations[ii][0].position.getAltitude() << " ";
			str_eastings << stations[ii][0].position.getEasting() << " ";
			str_northings << stations[ii][0].position.getNorthing() << " ";
		}
	}
	file << "YY MM DD HH " << str_altitudes.str() << "\n"; //altitudes
	file << "YY MM DD HH " << str_eastings.str() << "\n"; //easting
	file << "YY MM DD HH " << str_northings.str() << "\n"; //northing
	file << "YYYY MM DD HH";
	for(unsigned int ii=0; ii<sta_nr; ii++) {
		if(stations[ii].size() > 0) {
			file << " " << stations[ii][0].stationName;
		}
	}
	file << std::endl;
	return EXIT_SUCCESS;
}

void A3DIO::open2DFile(const std::vector< std::vector<StationData> >& stations,
                       const std::string& fileprefix, const std::string& label, const double& year, 
                       std::ofstream& file)
{//creates a meteo2D file according to the specifications
//the header is also written
	std::ostringstream out;
	out << year;

	const std::string filename = fileprefix+out.str()+".txt";

	file.open(filename.c_str(), ios::out | ios::trunc);
	if(!file) {
		throw FileAccessException("Can not create file "+filename, AT);
	}
	writeHeader(file, stations, label);
}

int A3DIO::write2DmeteoFile(const std::vector< std::vector<MeteoData> >& data,
                            const std::vector< std::vector<StationData> >& stations,
                            const unsigned int& parindex, const std::string& fileprefix,
                            const std::string& label)
{//HACK: we assume that all stations have data that is time synchronized...
	const unsigned int sta_nr = stations.size();
	if(sta_nr==0) return EXIT_FAILURE;
	const unsigned int nb_timesteps = data[0].size();
	if(nb_timesteps==0) return EXIT_FAILURE;

	std::ofstream file;
	int startyear, year, month, day, hour;
	data[0][0].date.getDate(startyear, month, day, hour);

	open2DFile(stations, fileprefix, label, startyear, file);
	file.flags ( ios::fixed );

	for(unsigned int ii=0; ii<nb_timesteps; ii++) {
		data[0][ii].date.getDate(year, month, day, hour);
		if(year!=startyear) {
			//if the year has changed, we need to write to a new file
			file.close();
			startyear = year;
			open2DFile(stations, fileprefix, label, year, file);
			file.flags ( ios::fixed );
		}

		file.fill('0');
		file << setw(4) << year << " " << setw(2) << month << " " << setw(2) << day << " " << setw(2) << hour;
		file.fill(' ');
		for(unsigned int j=0; j<sta_nr; j++) {
			double value = data[j][ii].param(parindex);
			if(value==IOUtils::nodata) {
				file << " " << setw(7) << setprecision(0) << IOUtils::nodata;
			} else {
				if(parindex==mio::MeteoData::TA) value = K_TO_C(value);
				if(parindex==mio::MeteoData::RH) value = value*100.;
				file << " " << setw(7) << setprecision(2) << value;
			}
		}
		file << "\n";
	}
	file.close();
	return EXIT_SUCCESS;
}

void A3DIO::write2DMeteo(const std::vector< std::vector<MeteoData> >& data, const std::vector< std::vector<StationData> >& stations)
{
	std::string tmp_path;
	cfg.getValue("METEOPATH", "Output", tmp_path);

	write2DmeteoFile(data, stations, mio::MeteoData::TA, tmp_path+"/tair", "air temperature");
	write2DmeteoFile(data, stations, mio::MeteoData::RH, tmp_path+"/rhum", "relative humidity");
	write2DmeteoFile(data, stations, mio::MeteoData::VW, tmp_path+"/wspd", "wind velocity");
	write2DmeteoFile(data, stations, mio::MeteoData::DW, tmp_path+"/wdir", "wind direction");
	write2DmeteoFile(data, stations, mio::MeteoData::HNW, tmp_path+"/prec", "precipitations");
}

/*extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const std::string& classname, const std::string& filename) {
		if(classname == "A3DIO") {
			cerr << "Creating handle to " << classname << endl;
			//return new A3DIO(deleteObject);
			return new A3DIO(deleteObject, filename);
		}
		cerr << "Could not load " << classname << endl;
		return NULL;
	}
}

*/

} //namespace
