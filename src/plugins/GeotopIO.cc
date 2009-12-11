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
 * - METEOPATH: string containing the path to the meteorological files
 * - METEOPREFIX: file name prefix for meteorological files
 */

using namespace std;

GeotopIO::GeotopIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

GeotopIO::GeotopIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	//Nothing else so far
}

GeotopIO::GeotopIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	//Nothing else so far
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
}

void GeotopIO::read2DGrid(Grid2DObject&, const string& filename)
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

void GeotopIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex)
{

	vector<string> tmpvec, vecColumnNames;
	string line="", filename="", path="", prefix="";

	(void)stationindex;
	vecMeteo.clear();
	vecStation.clear();

	cfg.getValue("METEOPATH", path); 
	cfg.getValue("METEOPREFIX", prefix);

	vector<StationData> myStations;
	/*
	 * read _meteo.txt to find out how many stations exist
	 * at what locations they are and what column headers to use
	 */
	readMetaData(myStations, vecColumnNames, path + "/" + prefix + ".txt");

	cout << "[i] GEOtopIO: Found " << myStations.size() << " station(s)" << endl;

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
		fin.open (filename.c_str(), ifstream::in);
		if (fin.fail()) {
			throw FileAccessException(filename, AT);
		}
  
		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
		
		//Go through file, save key value pairs
		try {
			getline(fin, line, eoln); //read complete line meta information
			unsigned int ncols = IOUtils::readLineToVec(line, tmpvec, ',');

			map<std::string, unsigned int> mapHeader;
			makeColumnMap(tmpvec, vecColumnNames, mapHeader);
			
			if (ncols == 0)
				throw InvalidFormatException("No meta data found in " + filename, AT);

			vector<double> tmpdata = vector<double>(ncols+1); //one extra for nodata value
			vector<int> ymdh = vector<int>(4);
			while (!fin.eof()){
				getline(fin, line, eoln); //read complete line of data
				//cout << line << endl;

				MeteoData md;
				if (IOUtils::readLineToVec(line, tmpvec, ',') != ncols) {
					break;
					//throw InvalidFormatException("Premature End " + filename, AT);
				}

				if (!IOUtils::convertString(ymdh[0], "20" + tmpvec[0].substr(6,2), std::dec)) //day
					throw InvalidFormatException(filename + ": " + line, AT);
				if (!IOUtils::convertString(ymdh[1], tmpvec[0].substr(3,2), std::dec)) //month
					throw InvalidFormatException(filename + ": " + line, AT);
				if (!IOUtils::convertString(ymdh[2], tmpvec[0].substr(0,2), std::dec)) //year
					throw InvalidFormatException(filename + ": " + line, AT);
				if (!IOUtils::convertString(ymdh[3], tmpvec[0].substr(9,2), std::dec)) //hour
					throw InvalidFormatException(filename + ": " + line, AT);

				for (unsigned int jj=1; jj<ncols; jj++) {
					if (!IOUtils::convertString(tmpdata[jj], tmpvec.at(jj), std::dec))
						throw InvalidFormatException(filename + ": " + line, AT);
				}
				tmpdata[ncols] = IOUtils::nodata;

				md.setMeteoData(Date_IO(ymdh[0],ymdh[1],ymdh[2],ymdh[3]), 
							 tmpdata[mapHeader["ta"]], 
							 tmpdata[mapHeader["iswr"]], 
							 tmpdata[mapHeader["vw"]], 
							 tmpdata[mapHeader["dw"]], 
							 tmpdata[mapHeader["rh"]], 
							 tmpdata[mapHeader["lwr"]], 
							 tmpdata[mapHeader["hnw"]], 
							 IOUtils::nodata,
							 IOUtils::nodata,
							 IOUtils::nodata,
							 IOUtils::nodata,
							 tmpdata[mapHeader["p"]]); 

				if ((md.date >= dateStart) && (md.date <= dateEnd)){
					convertUnits(md);
					vecMeteo[ii].push_back(md);
					vecStation[ii].push_back(myStations[ii]);
				} else {
					//cout << "Ignoring " << md.date.toString() << "  date.Start = " << dateStart.toString()
					//	<< "  date.End = " << dateEnd.toString() << endl;
				}
			}
		} catch(std::exception& e) {
			cleanup();
			throw;
		}
		fin.close();
	}
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
		string current="";
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
		case 11: mapHeader["lwr"] = tmpvec.size(); current="lwr"; break; 
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
	string line="";
	vector<string> tmpvec;

	if (!IOUtils::validFileName(metafile)) {
		throw InvalidFileNameException(metafile, AT);
	}
	if (!IOUtils::fileExists(metafile)) {
		throw FileNotFoundException(metafile, AT);
	}
  
	fin.clear();
	fin.open (metafile.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(metafile, AT);
	}
  
	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	try {
		while (!fin.eof()){
			getline(fin, line, eoln); //read complete line of data

			if (line.substr(0,2) == "1:"){//in section 1
				do {
					getline(fin, line, eoln); 
					unsigned int ncols = IOUtils::readLineToVec(line, tmpvec);

					if (ncols != 13){ 
						break;
					} else {
						vector<double> tmpdata=vector<double>(ncols);
						for (unsigned int jj=0; jj<ncols; jj++) {
							if (!IOUtils::convertString(tmpdata.at(jj), tmpvec.at(jj), std::dec))
								throw InvalidFormatException(metafile + ": " + line, AT);
						}
						vecStation.push_back(StationData(tmpdata[0], tmpdata[1], tmpdata[4], "", tmpdata[2], tmpdata[3]));
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

void GeotopIO::readAssimilationData(const Date_IO&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readSpecialPoints(CSpecialPTSArray&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::write2DGrid(const Grid2DObject&, const string& name)
{
	//Nothing so far
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::convertUnits(MeteoData& meteo)
{
	//converts C to Kelvin, converts lwr to ea, converts RH to [0,1]
	if(meteo.ta==nodata) {
		meteo.ta=nodata;
	} else {
		meteo.ta=C_TO_K(meteo.ta);
	}
	
	if(meteo.tsg==nodata) {
		meteo.tsg=nodata;
	} else {
		meteo.tsg=C_TO_K(meteo.tss);
	}
	
	if(meteo.tss==nodata) {
		meteo.tss=nodata;
	} else {
		meteo.tss=C_TO_K(meteo.tss);
	}

	if(meteo.rh==nodata) {
		meteo.rh=nodata;
	} else {
		meteo.rh /= 100.;
	}
}

extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
  
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "GeotopIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GeotopIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
