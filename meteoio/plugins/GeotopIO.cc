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
 * - precipitations in mm/h
 * - radiation in W/m²
 *
 * @section geotop_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - METEOPATH: string containing the path to the meteorological files
 * - METEOPREFIX: file name prefix for meteorological files
 */

const double GeotopIO::plugin_nodata = -9999.0; //plugin specific nodata value

GeotopIO::GeotopIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

GeotopIO::GeotopIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

GeotopIO::GeotopIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

GeotopIO::~GeotopIO() throw()
{
	cleanup();
}

void GeotopIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void GeotopIO::read2DGrid(Grid2DObject&, const std::string& filename)
{
	//Nothing so far
	(void)filename;
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::initParamNames(std::map<std::string, unsigned int>& mapParam)
{
	mapParam["Iprec"] = MeteoData::HNW;
	mapParam["WindS"] = MeteoData::VW;
	mapParam["WindDir"] = MeteoData::DW;
	mapParam["RelHum"] = MeteoData::RH;
	mapParam["AirT"] = MeteoData::TA;
	mapParam["AirP"] = MeteoData::P;
	mapParam["SWglobal"] = MeteoData::ISWR;
}


void GeotopIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
						const std::vector< std::vector<StationData> >& vecStation,
						const std::string&)
{
	string path="";
	vector<string> vecSequence;
	vector<int> ymdhm = vector<int>(5);
	map<string, unsigned int> mapParam;
	initParamNames(mapParam);
	cfg.getValue("METEOPATH", "Output", path);
	cfg.getValue("METEOSEQ", "Output", vecSequence);

	//Check whether vecSequence is valid, that is the keys are part of mapParam
	for (unsigned int ii=0; ii<vecSequence.size(); ii++){
		map<string,unsigned int>::iterator it = mapParam.find(vecSequence[ii]);
		if (it == mapParam.end())
			throw InvalidFormatException("Key " + vecSequence[ii] +" invalid in io.ini:METEODESTSEQ", AT);
	}

	//write the meta data file _meteo.txt
	fout.open(string(path + "/_meteo.txt").c_str());
	fout << "/* Automatically generated by MeteoIO */" << endl << endl;
	fout << "index{2}" << endl << endl;
	fout << "1: double matrix meteo_station{" << vecStation.size() << ",13}" << endl;
	for (unsigned int ii=0; ii<vecStation.size(); ii++){
		if (vecStation.at(ii).size()>0){
			Coords coord = vecStation.at(ii).at(0).position;
			coord.setProj(coordout, coordoutparam); //Setting the output projection
			fout.precision(12);
			fout << coord.getEasting() << "\t" << coord.getNorthing() << "\t"
				<< coord.getLat()<< "\t" << coord.getLon() << "\t" << coord.getAltitude() << "\t"
				<< plugin_nodata << "\t" << plugin_nodata << "\t" << plugin_nodata << "\t"
				<< plugin_nodata << "\t" << plugin_nodata << "\t" << plugin_nodata << "\t"
				<< plugin_nodata << "\t" << plugin_nodata << endl;
		}
	}
	fout << endl << "2: stringbin metocolnames" << endl
		<< "{Iprec, WindS, WindDir, RelHum, AirT, AirP, SWglobal, SWdirect, SWdiffuse, TauCloud, Cloud, LWin, SWnet, Tsup}" << endl;
	fout.close(); //finished writing meta data

	//Writing actual meteo files
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		stringstream ss;
		ss.fill('0');
		ss << path << "/" << "_meteo" << setw(4) << (ii+1) << ".txt";

		fout.open(ss.str().c_str());
		if (fout.fail()) throw FileAccessException(ss.str().c_str(), AT);

		fout << fixed << showpoint << setprecision(5);

		for (unsigned int jj=0; jj<vecMeteo.at(ii).size(); jj++){
			ss.str(""); //clear the stringstream
			vecMeteo[ii][jj].date.getDate(ymdhm[0], ymdhm[1], ymdhm[2], ymdhm[3], ymdhm[4]);

			//the date will be written in the form "DD/MM/YYYY hh:mm"
			ss << setw(2) << ymdhm[2] << "/" << setw(2) << ymdhm[1] << "/" << ymdhm[0] << " " // DD/MM/YYYY
			   << setw(2) << ymdhm[3] << ":" << setw(2) << ymdhm[4];                          // hh:mm

			MeteoData tmpmd = vecMeteo[ii][jj];
			convertUnitsBack(tmpmd);

			for (unsigned int kk=0; kk<vecSequence.size(); kk++){
				if (jj==0){ //This is for writing the header
					if (kk==0) fout << "Date";
					fout << "," << vecSequence[kk];
					if (kk==(vecSequence.size()-1)) fout << endl;
				}

				//Write all the data, make sure to transform the nodata values correctly
				if (tmpmd.param(mapParam[vecSequence[kk]]) == IOUtils::nodata)
					ss << ", " << plugin_nodata;
				else
					ss << ", " << setprecision(7) << tmpmd.param(mapParam[vecSequence[kk]]);
			}
			fout << ss.str() << endl;
		}

		fout.close();
	}
}

void GeotopIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	string path="", prefix="";
	vector<string> tmpvec, vecColumnNames;

	vecStation.clear();

	cfg.getValue("METEOPATH", "Input", path);
	cfg.getValue("METEOPREFIX", "Input", prefix);

	readMetaData(vecStation, vecColumnNames, path + "/" + prefix + ".txt");
}

void GeotopIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
							  std::vector< std::vector<MeteoData> >& vecMeteo,
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex)
{
	vector<std::string> tmpvec, vecColumnNames;
	string line="", filename="", path="", prefix="";

	(void)stationindex;
	vecMeteo.clear();
	vecStation.clear();

	cfg.getValue("METEOPATH", "Input", path);
	cfg.getValue("METEOPREFIX", "Input", prefix);

	vector<StationData> myStations;
	/*
	 * read _meteo.txt to find out how many stations exist
	 * at what locations they are and what column headers to use
	 */
	readMetaData(myStations, vecColumnNames, path + "/" + prefix + ".txt");

	std::cout << "[I] GEOtopIO: Found " << myStations.size() << " station(s)" << std::endl;

	for (unsigned int ii=0; ii<myStations.size(); ii++) {
		vecMeteo.push_back( vector<MeteoData>() );
		vecStation.push_back( vector<StationData>() );

		stringstream ss;
		ss.fill('0');
		ss << path << "/" << prefix << setw(4) << (ii+1) << ".txt";

		filename = ss.str();
		//cout << ss.str() << endl;

		if (!IOUtils::validFileName(filename)) {
			throw InvalidFileNameException(filename, AT);
		}

		if (!IOUtils::fileExists(filename)) {
			throw FileNotFoundException(filename, AT);
		}

		fin.clear();
		fin.open (filename.c_str(), std::ifstream::in);
		if (fin.fail()) {
			throw FileAccessException(filename, AT);
		}

		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		//Go through file, save key value pairs
		try {
			getline(fin, line, eoln); //read complete line meta information
			unsigned int ncols = IOUtils::readLineToVec(line, tmpvec, ',');

			std::map<std::string, unsigned int> mapHeader;
			makeColumnMap(tmpvec, vecColumnNames, mapHeader);

			if (ncols == 0)
				throw InvalidFormatException("No meta data found in " + filename, AT);

			std::vector<double> tmpdata = std::vector<double>(ncols+1); //one extra for nodata value
			std::vector<int> ymdh = std::vector<int>(4);
			while (!fin.eof()){
				getline(fin, line, eoln); //read complete line of data

				MeteoData md;
				if (IOUtils::readLineToVec(line, tmpvec, ',') != ncols) {
					break;
					//throw InvalidFormatException("Premature End " + filename, AT);
				}

				//tmpvec[0] holds the date in many possible formats -> needs to be parsed
				parseDate(tmpvec.at(0), filename+": "+line, md.date);

				for (unsigned int jj=1; jj<ncols; jj++) {
					if (!IOUtils::convertString(tmpdata[jj], tmpvec.at(jj), std::dec))
						throw InvalidFormatException(filename + ": " + line, AT);
				}
				tmpdata[ncols] = IOUtils::nodata;

				md.setData(MeteoData::TA, tmpdata[mapHeader["ta"]]);
				md.setData(MeteoData::ISWR, tmpdata[mapHeader["iswr"]]);
				md.setData(MeteoData::VW, tmpdata[mapHeader["vw"]]);
				md.setData(MeteoData::DW, tmpdata[mapHeader["dw"]]);
				md.setData(MeteoData::RH, tmpdata[mapHeader["rh"]]);
				md.setData(MeteoData::ILWR, tmpdata[mapHeader["ilwr"]]);
				md.setData(MeteoData::HNW, tmpdata[mapHeader["hnw"]]);
				md.setData(MeteoData::P, tmpdata[mapHeader["p"]]);

				if ((md.date >= dateStart) && (md.date <= dateEnd)){
					convertUnits(md);
					vecMeteo[ii].push_back(md);
					vecStation[ii].push_back(myStations[ii]);
				}
			}
		} catch(std::exception& e) {
			cleanup();
			throw;
		}
		fin.close();
	}
}

void GeotopIO::parseDate(const std::string& datestring, const std::string& fileandline, Date& date)
{
	/*
	 * In order to be more flexible with the date parsing in GEOtop meteo files,
	 * this function will allow any date format common to GEOtop to be accepted
	 * examples for valid dates: 14/4/09 8:1, 14/04/2009 08:01, 14-4-2009 8:01
	 */

	std::vector<int> ymdhm = std::vector<int>(5);

	//parsing the day
	size_t found1 = datestring.find_first_of("/-",0);
	if (!IOUtils::convertString(ymdhm.at(2), datestring.substr(0, found1), std::dec)) //day
		throw InvalidFormatException(fileandline, AT);

	//parsing the month
	size_t found2 = datestring.find_first_of("/-", found1+1);
	if (!IOUtils::convertString(ymdhm.at(1), datestring.substr(found1+1, found2-found1-1), std::dec)) //month
		throw InvalidFormatException(fileandline, AT);

	//parsing the year: possibly prefix of '20' necessary (if the year just states '09' for example)
	size_t found3 = datestring.find_first_of(" ", found2+1);
	string year = datestring.substr(found2+1, found3-found2-1);
	if (year.length() == 2) year = "20" + year; //add year prefix
	if (!IOUtils::convertString(ymdhm.at(0), year, std::dec)) //year
		throw InvalidFormatException(fileandline, AT);

	//parsing hour and minute
	size_t found4 = datestring.find_first_of(":", found3+1);
	if (!IOUtils::convertString(ymdhm.at(3), datestring.substr(found3+1,found4-found3-1), std::dec)) //month
		throw InvalidFormatException(fileandline, AT);
	if (!IOUtils::convertString(ymdhm.at(4), datestring.substr(found4+1,datestring.length()-found4-1), std::dec)) //month
		throw InvalidFormatException(fileandline, AT);

	date.setDate(ymdhm[0], ymdhm[1], ymdhm[2], ymdhm[3], ymdhm[4]);
}

void GeotopIO::makeColumnMap(const std::vector<std::string>& tmpvec,
					    const std::vector<std::string>& vecColumnNames,
					    std::map<std::string, unsigned int>& mapHeader)
{
	/*
	  #1 Precipitation intensity (mm/h)
	  #2 Wind speed (m/s)
	  #3 Direction from which wind comes from (degree from North, clockwise)
	  #4 Relative humidity (%)
	  #5 Air temperature (C)
	  #6 Air pressure (mbar)
	  #7 Global shortwave radiation (W/m2)
	  #8 Direct shortwave radiation (W/m2)
	  #9 Diffuse shortwave radiation (W/m2)
	  #10 Cloudiness transmissivity
	  #11 Cloudiness (fraction from 0 to 1)
	  #12 Incoming longwave radiation (W/m2)
	  #13 Net shortwave radiation (W/m2)
	  #14 Temperature of the soil surface (for Dirichlet conditions)
	*/

	for (unsigned int ii=0; ii<vecColumnNames.size(); ii++){
		std::string current="";
		switch(ii){
		case 0: mapHeader["hnw"] = tmpvec.size(); current="hnw"; break;
		case 1: mapHeader["vw"] = tmpvec.size(); current="vw"; break;
		case 2: mapHeader["dw"] = tmpvec.size(); current="dw"; break;
		case 3: mapHeader["rh"] = tmpvec.size(); current="rh"; break;
		case 4: mapHeader["ta"] = tmpvec.size(); current="ta"; break;
		case 5: mapHeader["p"] = tmpvec.size(); current="p"; break;
		case 6: mapHeader["iswr"] = tmpvec.size(); current="iswr"; break;
		case 7: mapHeader["dirswr"] = tmpvec.size(); current="dirswr"; break;
		case 8: mapHeader["diffswr"] = tmpvec.size(); current="diffswr"; break;
		case 9: mapHeader["cloudt"] = tmpvec.size(); current="cloudt"; break;
		case 10: mapHeader["cloudi"] = tmpvec.size(); current="cloudi"; break;
		case 11: mapHeader["ilwr"] = tmpvec.size(); current="ilwr"; break;
		case 12: mapHeader["nswr"] = tmpvec.size(); current="nswr"; break;
		case 13: mapHeader["tsup"] = tmpvec.size(); current="tsup"; break;
		default: throw IOException("GEOtopIO can only deal with 14 meteo parameters", AT); break;
		}

		//Go through the columns and seek out which parameter corresponds with which column
		for (unsigned int jj=1; jj<tmpvec.size(); jj++){
			if (tmpvec[jj].length() >= vecColumnNames[ii].length())
				if (vecColumnNames[ii] == tmpvec[jj].substr(0, vecColumnNames[ii].length()))
					mapHeader[current] = jj;
		}
	}
}

void GeotopIO::readMetaData(std::vector<StationData>& vecStation, std::vector<std::string>& vecColumnNames,
					   const std::string& metafile)
{
	std::string line="";
	std::vector<std::string> tmpvec;
	unsigned int stationNumber = 1; //Since the stations don't have a name, they will be numbered

	if (!IOUtils::validFileName(metafile)) {
		throw InvalidFileNameException(metafile, AT);
	}
	if (!IOUtils::fileExists(metafile)) {
		throw FileNotFoundException(metafile, AT);
	}

	fin.clear();
	fin.open (metafile.c_str(), std::ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(metafile, AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	try {
		Coords coordinate(coordin, coordinparam);
		while (!fin.eof()){
			getline(fin, line, eoln); //read complete line of data

			if (line.substr(0,2) == "1:"){//in section 1
				do {
					getline(fin, line, eoln);
					unsigned int ncols = IOUtils::readLineToVec(line, tmpvec);

					if (ncols != 13){
						break;
					} else {
						std::vector<double> tmpdata=std::vector<double>(ncols);
						for (unsigned int jj=0; jj<ncols; jj++) {
							if (!IOUtils::convertString(tmpdata.at(jj), tmpvec.at(jj), std::dec))
								throw InvalidFormatException(metafile + ": " + line, AT);
						}
						//HACK!! would it be possible for getValueForKey() to do this transparently? (with a user flag)
						tmpdata[0] = IOUtils::standardizeNodata(tmpdata[0], plugin_nodata);
						tmpdata[1] = IOUtils::standardizeNodata(tmpdata[1], plugin_nodata);
						tmpdata[2] = IOUtils::standardizeNodata(tmpdata[2], plugin_nodata);
						tmpdata[3] = IOUtils::standardizeNodata(tmpdata[3], plugin_nodata);
						tmpdata[4] = IOUtils::standardizeNodata(tmpdata[4], plugin_nodata);

						coordinate.setLatLon(tmpdata[2], tmpdata[3], tmpdata[4], false);
						coordinate.setXY(tmpdata[0], tmpdata[1], tmpdata[4], false);
						try {
							coordinate.check();
						} catch(...) {
							std::cerr << "[E] Error in geographic coordinates in file " << metafile << " trapped at " << AT << std::endl;
							throw;
						}

						stringstream ss;
						ss << "Station_" << stationNumber;
						vecStation.push_back( StationData(coordinate, ss.str()) );

						stationNumber++; //Stationnames are simply a sequence of ascending numbers
					}
					IOUtils::trim(line);
				} while ((line.substr(0,2) != "/*") && (line!="") && (!fin.eof()));
			} else if (line.substr(0,2) == "2:"){//in section 2
				getline(fin, line, eoln);
				IOUtils::trim(line);
				if (line.length()>2) line = line.substr(1,line.length()-2);

				IOUtils::readLineToVec(line,vecColumnNames, ',');
				for (unsigned int ii=0; ii<vecColumnNames.size(); ii++){
					IOUtils::trim(vecColumnNames[ii]); //trim the column headers
				}
			}
		}
	} catch(std::exception& e) {
		cleanup();
		throw;
	}
	cleanup();
}

void GeotopIO::readAssimilationData(const Date&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::write2DGrid(const Grid2DObject&, const std::string& name)
{
	//Nothing so far
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::convertUnitsBack(MeteoData& meteo)
{
	//converts Kelvin to C, converts RH to [0,100]
	if(meteo.ta!=IOUtils::nodata) {
		meteo.ta=K_TO_C(meteo.ta);
	}

	if(meteo.tsg!=IOUtils::nodata) {
		meteo.tsg=K_TO_C(meteo.tsg);
	}

	if(meteo.tss!=IOUtils::nodata) {
		meteo.tss=K_TO_C(meteo.tss);
	}

	if(meteo.rh!=IOUtils::nodata) {
		meteo.rh *= 100.;
	}
}

void GeotopIO::convertUnits(MeteoData& meteo)
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

#ifndef _METEOIO_JNI
extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const std::string& classname, const Config& cfg) {
		if(classname == "GeotopIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GeotopIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
