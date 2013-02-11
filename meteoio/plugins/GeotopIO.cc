/***********************************************************************************/
/*  Copyright 2009 EPFL                                                            */
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
#include "GeotopIO.h"

using namespace std;

namespace mio {
/**
 * @page geotop GEOTOP
 * @section geotop_format Format
 * This plugin reads Legacy Geotop meteorological input data.
 *
 * @section geotop_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitation in mm/h
 * - radiation in W/m²
 *
 * @section geotop_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - METAFILE:  string containing the absolute filename of the geotop.inpts file in the [Input] section
 * - METEOPATH: string containing the path to the meteorological files
 * - METEOPREFIX: file name prefix for meteorological files
 * - METEOSEQ:specifiy in which order the columns should be printed out! ONLY relevant for writing out
 */

const double GeotopIO::plugin_nodata = -9999.0; //plugin specific nodata value

GeotopIO::GeotopIO(const std::string& configfile)
         : cfg(configfile), in_tz(0.), out_tz(0.), nr_of_stations(IOUtils::npos),
           fin(), fout(), vec_streampos(), vecStation(), mapColumnNames(),
           coordin(), coordinparam(), coordout(), coordoutparam()
{
           IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
           cfg.getValue("TIME_ZONE", "Input", in_tz, Config::nothrow);
           cfg.getValue("TIME_ZONE", "Output", out_tz, Config::nothrow);
}

GeotopIO::GeotopIO(const Config& cfgreader)
         : cfg(cfgreader), in_tz(0.), out_tz(0.), nr_of_stations(IOUtils::npos),
           fin(), fout(), vec_streampos(), vecStation(), mapColumnNames(),
           coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	cfg.getValue("TIME_ZONE", "Input", in_tz, Config::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_tz, Config::nothrow);
}

GeotopIO::~GeotopIO() throw () {
	cleanup();
}

void GeotopIO::cleanup() throw () {
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void GeotopIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readDEM(DEMObject&) {
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readLanduse(Grid2DObject&) {
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::initParamNames(std::map<std::string, size_t>& mapParam) {
	mapParam["Iprec"] = MeteoData::HNW;
	mapParam["WindSp"] = MeteoData::VW;
	mapParam["WindDir"] = MeteoData::DW;
	mapParam["RH"] = MeteoData::RH;
	mapParam["AirT"] = MeteoData::TA;
	mapParam["AirP"] = MeteoData::P;
	mapParam["Swglob"] = MeteoData::ISWR;
}

void GeotopIO::writeMeteoData(
		const std::vector<std::vector<MeteoData> >& vecMeteo,
		const std::string&) {
	string path;
	vector<string> vecSequence;
	vector<int> ymdhm = vector<int> (5);
	map<string, size_t> mapParam;
	initParamNames(mapParam);
	cfg.getValue("METEOPATH", "Output", path);
	cfg.getValue("METEOSEQ", "Output", vecSequence);

	//Check whether vecSequence is valid, that is the keys are part of mapParam
	for (size_t ii = 0; ii < vecSequence.size(); ii++) {
		map<string, size_t>::iterator it = mapParam.find(vecSequence[ii]);
		if (it == mapParam.end())
			throw InvalidFormatException("Key " + vecSequence[ii]
					+ " invalid in io.ini:METEOSEQ", AT);
	}

	//write the meta data file _meteo.txt
	fout.open(string(path + "/_meteo.txt").c_str());
	fout << "/* Automatically generated by MeteoIO */\n\n";
	fout << "1: double matrix meteo_station{" << vecMeteo.size() << ",13}\n";
	for (size_t ii = 0; ii < vecMeteo.size(); ii++) {
		if (vecMeteo.at(ii).size() > 0) {
			Coords coord = vecMeteo.at(ii).at(0).meta.position;
			coord.setProj(coordout, coordoutparam); //Setting the output projection
			fout.precision(12);
			fout << coord.getEasting() << "\t" << coord.getNorthing() << "\t"
					<< coord.getLat() << "\t" << coord.getLon() << "\t"
					<< coord.getAltitude() << "\t" << plugin_nodata << "\t"
					<< plugin_nodata << "\t" << plugin_nodata << "\t"
					<< plugin_nodata << "\t" << plugin_nodata << "\t"
					<< plugin_nodata << "\t" << plugin_nodata << "\t"
					<< plugin_nodata << endl;
		}
	}
	fout << "\n2: stringbin metocolnames\n"
			<< "{Iprec, WindSp, WindDir, RH, AirT, AirP, Swglob, SWdirect, SWdiffuse, TauCloud, Cloud, LWin, SWnet, Tsup}\n";
	fout.close(); //finished writing meta data

	//Writing actual meteo files
	for (size_t ii = 0; ii < vecMeteo.size(); ii++) {
		stringstream ss;
		ss.fill('0');
		ss << path << "/" << "_meteo" << setw(4) << (ii + 1) << ".txt";

		fout.open(ss.str().c_str());
		if (fout.fail())
			throw FileAccessException(ss.str().c_str(), AT);

		fout << fixed << showpoint << setprecision(5);

		for (size_t jj = 0; jj < vecMeteo.at(ii).size(); jj++) {
			ss.str(""); //clear the stringstream
			Date tmp_date(vecMeteo[ii][jj].date);
			tmp_date.setTimeZone(out_tz);
			tmp_date.getDate(ymdhm[0], ymdhm[1], ymdhm[2], ymdhm[3], ymdhm[4]);

			//the date will be written in the form "DD/MM/YYYY hh:mm"
			ss << setw(2) << ymdhm[2] << "/" << setw(2) << ymdhm[1] << "/"
					<< ymdhm[0] << " " // DD/MM/YYYY
					<< setw(2) << ymdhm[3] << ":" << setw(2) << ymdhm[4]; // hh:mm

			MeteoData tmpmd = vecMeteo[ii][jj];
			convertUnitsBack(tmpmd);

			for (size_t kk = 0; kk < vecSequence.size(); kk++) {
				if (jj == 0) { //This is for writing the header
					if (kk == 0)
						fout << "Date";
					fout << "," << vecSequence[kk];
					if (kk == (vecSequence.size() - 1))
						fout << "\n";
				}

				//Write all the data, make sure to transform the nodata values correctly
				if (tmpmd(mapParam[vecSequence[kk]]) == IOUtils::nodata)
					ss << ", " << plugin_nodata;
				else
					ss << ", " << setprecision(7) << tmpmd(
							mapParam[vecSequence[kk]]);
			}
			fout << ss.str() << "\n";
		}

		fout.close();
	}
}

void GeotopIO::readStationData(const Date&, std::vector<StationData>& vecMeta) {
	string metafile;
	vecMeta.clear();

	if (vecStation.empty()) {
		cfg.getValue("METAFILE", "Input", metafile);
		readMetaData(metafile);
	}

	vecMeta = vecStation;
}

void GeotopIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
		std::vector<std::vector<MeteoData> >& vecMeteo, const size_t& /*stationindex*/) {
	vector<std::string> tmpvec;
	string line, filename, path, prefix;
	vecMeteo.clear();

	cfg.getValue("METEOPATH", "Input", path);
	cfg.getValue("METEOPREFIX", "Input", prefix);

	//read geotop.inpts to find out how many stations exist
	//at what locations they are and what column headers to use
	vector<StationData> myStations;
	readStationData(dateStart, myStations);

	if (vec_streampos.empty()) //the vec_streampos save file pointers for certain dates
		vec_streampos = vector<map<Date, std::streampos> > (vecStation.size());

	if (nr_of_stations == IOUtils::npos)
		nr_of_stations = vecStation.size();

	cerr << "[i] GEOtopIO: Found " << nr_of_stations << " station(s)" << std::endl;

	for (size_t ii = 0; ii < nr_of_stations; ii++) {
		vecMeteo.push_back(vector<MeteoData> ());

		stringstream ss;
		ss.fill('0');
		ss << path << "/" << prefix << setw(4) << (ii + 1) << ".txt";

		filename = ss.str();

		if (!IOUtils::validFileName(filename))
			throw InvalidFileNameException(filename, AT);

		if (!IOUtils::fileExists(filename))
			throw FileNotFoundException(filename, AT);

		fin.clear();
		fin.open(filename.c_str(), std::ifstream::in);
		if (fin.fail())
			throw FileAccessException(filename, AT);

		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		//Go through file, save key value pairs
		try {
			getline(fin, line, eoln); //read complete line meta information
			size_t ncols = IOUtils::readLineToVec(line, tmpvec, ',');

			std::vector<size_t> indices;
			MeteoData md;
			identify_fields(tmpvec, filename, indices, md);
			md.meta = vecStation[ii];

			if (ncols == 0)
				throw InvalidFormatException("No meta data found in " + filename, AT);

			std::vector<double> tmpdata = std::vector<double>(ncols + 1); //one extra for nodata value

			//The following 4 lines are an optimization to jump to the correct position in the file
			streampos current_fpointer = -1; //the filepointer for the current date
			map<Date, streampos>::const_iterator it = vec_streampos.at(ii).find(dateStart);
			if (it != vec_streampos.at(ii).end())
				fin.seekg(it->second); //jump to position in the file

			while (!fin.eof()) {
				streampos tmp_fpointer = fin.tellg();
				getline(fin, line, eoln); //read complete line of data

				if (IOUtils::readLineToVec(line, tmpvec, ',') != ncols) {
					break;
					//throw InvalidFormatException("Premature End " + filename, AT);
				}

				md.reset(); //clear all previous data
				//tmpvec[0] holds the date in many possible formats -> needs to be parsed
				parseDate(tmpvec.at(0), filename + ": " + line, md.date);

				if ((md.date >= dateStart) && (md.date <= dateEnd)) {
					current_fpointer = tmp_fpointer;

					for (size_t jj = 1; jj < ncols; jj++) {
						if (!IOUtils::convertString(tmpdata[jj], tmpvec.at(jj), std::dec))
							throw InvalidFormatException(filename + ": " + line, AT);

						if (indices[jj-1] != IOUtils::npos){
							md(indices[jj-1]) = tmpdata[jj];
						}
					}

					convertUnits(md);
					vecMeteo[ii].push_back(md);
				} else if (md.date > dateEnd) {
					break;
				}
			}

			//save stream position and the corresponding end date
			if (current_fpointer != static_cast<streampos> (-1))
				vec_streampos.at(ii)[dateEnd] = current_fpointer;
		} catch (const std::exception&) {
			cleanup();
			throw;
		}
		fin.close();
	}
}

void GeotopIO::parseDate(const std::string& datestring,
		const std::string& fileandline, Date& date) {
	/*
	 * In order to be more flexible with the date parsing in GEOtop meteo files,
	 * this function will allow any date format common to GEOtop to be accepted
	 * examples for valid dates: 14/4/09 8:1, 14/04/2009 08:01, 14-4-2009 8:01
	 */

	std::vector<int> ymdhm = std::vector<int>(5);

	//parsing the day
	size_t found1 = datestring.find_first_of("/-", 0);
	if (!IOUtils::convertString(ymdhm.at(2), datestring.substr(0, found1),
			std::dec)) //day
		throw InvalidFormatException(fileandline, AT);

	//parsing the month
	size_t found2 = datestring.find_first_of("/-", found1 + 1);
	if (!IOUtils::convertString(ymdhm.at(1), datestring.substr(found1 + 1,
			found2 - found1 - 1), std::dec)) //month
		throw InvalidFormatException(fileandline, AT);

	//parsing the year: possibly prefix of '20' necessary (if the year just states '09' for example)
	size_t found3 = datestring.find_first_of(" ", found2 + 1);
	string year = datestring.substr(found2 + 1, found3 - found2 - 1);
	if (year.length() == 2)
		year = "20" + year; //add year prefix
	if (!IOUtils::convertString(ymdhm.at(0), year, std::dec)) //year
		throw InvalidFormatException(fileandline, AT);

	//parsing hour and minute
	size_t found4 = datestring.find_first_of(":", found3 + 1);
	if (!IOUtils::convertString(ymdhm.at(3), datestring.substr(found3 + 1,
			found4 - found3 - 1), std::dec)) //month
		throw InvalidFormatException(fileandline, AT);
	if (!IOUtils::convertString(ymdhm.at(4), datestring.substr(found4 + 1,
			datestring.length() - found4 - 1), std::dec)) //month
		throw InvalidFormatException(fileandline, AT);

	date.setDate(ymdhm[0], ymdhm[1], ymdhm[2], ymdhm[3], ymdhm[4], in_tz);
}

void GeotopIO::identify_fields(const std::vector<std::string>& tmpvec, const std::string& filename,
                               std::vector<size_t>& indices, MeteoData& md) {
	//Go through the columns and seek out which parameter corresponds with which column
	for (size_t jj = 1; jj < tmpvec.size(); jj++) { //skip field 1, that one is reserved for the date
		std::map<std::string, size_t>::iterator it = mapColumnNames.find(tmpvec[jj]);
		if (it != mapColumnNames.end()) {
			size_t index = it->second;
			if (index == IOUtils::npos){
				index = md.addParameter(it->first);
			}
			indices.push_back(index);
			//cerr << tmpvec[jj] << "=> index " << index << endl;
		} else {
			indices.push_back(IOUtils::npos);
			cerr << "[w] Column '" << tmpvec[jj] << "' in file " << filename << " will be ignored!" << endl;
			//throw InvalidFormatException("The column name " + tmpvec[jj] + " is unknown", AT);
		}
	}
}

void GeotopIO::parseMetaData(const std::string& head, const std::string& datastr, std::vector<std::string>& tmpvec) {
	tmpvec.clear();
	string mdata = datastr.substr(head.length() + 1, datastr.length() + 1);
	IOUtils::readLineToVec(mdata, tmpvec, ',');
}

void GeotopIO::readMetaData(const std::string& metafile) {
	std::string line;

	//vector indexes correspond to meteo data
	int x = 0, y = 1, lat = 2, lon = 3, elv = 4;

	std::vector<std::string> tmpvec;
	unsigned int stationNumber = 1; //Since the stations don't have a name, they will be numbered

	if (!IOUtils::validFileName(metafile)) {
		throw InvalidFileNameException(metafile, AT);
	}
	if (!IOUtils::fileExists(metafile)) {
		throw FileNotFoundException(metafile, AT);
	}

	fin.clear();
	fin.open(metafile.c_str(), std::ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(metafile, AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	std::vector<std::string> vecX, vecY, vecLat, vecLon, vecAlt;

	size_t meta_counter = 0;

	try {
		Coords coordinate(coordin, coordinparam);
		while (!fin.eof()) {
			getline(fin, line, eoln); //read complete line of data

			if (line.find("MeteoStationCoordinateX") == 0) {//in section 1
				parseMetaData("MeteoStationCoordinateX", line, tmpvec);
				vecX = tmpvec;
				meta_counter |= 1;
			} else if (line.find("MeteoStationCoordinateY") == 0) {
				parseMetaData("MeteoStationCoordinateY", line, tmpvec);
				vecY = tmpvec;
				meta_counter |= 2;
			} else if (line.find("MeteoStationLatitude") == 0) {
				parseMetaData("MeteoStationLatitude", line, tmpvec);
				vecLat = tmpvec;
				meta_counter |= 4;
			} else if (line.find("MeteoStationLongitude") == 0) {
				parseMetaData("MeteoStationLongitude", line, tmpvec);
				vecLon = tmpvec;
				meta_counter |= 8;
			} else if (line.find("MeteoStationElevation") == 0) {
				parseMetaData("MeteoStationElevation", line, tmpvec);
				vecAlt = tmpvec;
				meta_counter |= 16;
			} else if (line.find("HeaderIPrec") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::HNW;
			} else if (line.find("HeaderWindVelocity") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::VW;
			} else if (line.find("HeaderWindDirection") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::DW;
			} else if (line.find("HeaderRH") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::RH;
			} else if (line.find("HeaderAirTemp") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::TA;
			} else if (line.find("HeaderSWglobal") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::ISWR;
			} else if (line.find("HeaderCloudSWTransmissivity") == 0) {
				mapColumnNames[getValueForKey(line)] = MeteoData::RSWR;
			} else if (line.find("NumberOfMeteoStations") == 0) {
				size_t pos = line.find("=");
				string val = line.substr(pos+1);
				IOUtils::trim(val);

				if (!IOUtils::convertString(nr_of_stations, val,	std::dec))
					throw InvalidFormatException(metafile + ": " + line, AT);
			}
		}

		//Check that all meta information is available
		if (meta_counter != 31)
			throw InvalidFormatException("Your GEOtop METAFILE " + metafile + " does not contain all required meta data", AT);

		//Check for consistency between the meta data data vectors
		if ((vecX.size() != vecY.size()) || (vecY.size() != vecLat.size()) || (vecLat.size() != vecLon.size()) || (vecLon.size() != vecAlt.size()))
			throw InvalidFormatException("Your GEOtop METAFILE " + metafile + " does not contain a consistent number of meta data fields", AT);

		for (unsigned int i = 0; i < vecX.size(); i++) {
			std::vector<double> tmpdata = std::vector<double>(5);

			if (!IOUtils::convertString(tmpdata.at(x), vecX.at(i),	std::dec))
				throw InvalidFormatException(metafile + ": " + line, AT);
			if (!IOUtils::convertString(tmpdata.at(y), vecY.at(i), std::dec))
				throw InvalidFormatException(metafile + ": " + line, AT);
			if (!IOUtils::convertString(tmpdata.at(lat),	vecLat.at(i), std::dec))
				throw InvalidFormatException(metafile + ": " + line, AT);
			if (!IOUtils::convertString(tmpdata.at(lon),	vecLon.at(i), std::dec))
				throw InvalidFormatException(metafile + ": " + line, AT);
			if (!IOUtils::convertString(tmpdata.at(elv),	vecAlt.at(i), std::dec))
				throw InvalidFormatException(metafile + ": " + line, AT);

			//HACK!! would it be possible for getValueForKey() to do this transparently? (with a user flag)
			tmpdata[0] = IOUtils::standardizeNodata(tmpdata[0], plugin_nodata);
			tmpdata[1] = IOUtils::standardizeNodata(tmpdata[1], plugin_nodata);
			tmpdata[2] = IOUtils::standardizeNodata(tmpdata[2], plugin_nodata);
			tmpdata[3] = IOUtils::standardizeNodata(tmpdata[3], plugin_nodata);
			tmpdata[4] = IOUtils::standardizeNodata(tmpdata[4], plugin_nodata);

			//coordinate.setLatLon(tmpdata[2], tmpdata[3], tmpdata[4], false);
			coordinate.setXY(tmpdata[0], tmpdata[1], tmpdata[4], true);

			try {
				coordinate.check();
			} catch (...) {
				std::cerr << "[E] Error in geographic coordinates in file " << metafile << " trapped at " << AT << std::endl;
				throw;
			}

			stringstream ss;
			ss << "Station_" << stationNumber;
			vecStation.push_back(StationData(coordinate, ss.str()));

			stationNumber++; //Stationnames are simply a sequence of ascending numbers

			IOUtils::trim(line);
		}
	} catch (const std::exception&) {
		cleanup();
		throw;
	}
	cleanup();
}

std::string GeotopIO::getValueForKey(const std::string& line)
{
	size_t pos_start = line.find("\"");
	size_t pos_end = line.find("\"", pos_start+1);

	if ((pos_start != string::npos) && (pos_end != string::npos)) {
		string param_name = line.substr(pos_start+1, pos_end - pos_start - 1);
		IOUtils::trim(param_name);
		return param_name;
	}

	return "";
}

void GeotopIO::readAssimilationData(const Date&, Grid2DObject&) {
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readSpecialPoints(std::vector<Coords>&) {
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::convertUnitsBack(MeteoData& meteo) {
	//converts Kelvin to C, converts RH to [0,100]
	double& ta = meteo(MeteoData::TA);
	if (ta != IOUtils::nodata)
		ta = K_TO_C(ta);

	double& tsg = meteo(MeteoData::TSG);
	if (tsg != IOUtils::nodata)
		tsg = K_TO_C(tsg);

	double& tss = meteo(MeteoData::TSS);
	if (tss != IOUtils::nodata)
		tss = K_TO_C(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
		rh *= 100.;
}

void GeotopIO::convertUnits(MeteoData& meteo) {
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts RH to [0,1]
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
}

} //namespace
