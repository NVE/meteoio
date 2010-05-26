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
#include "SNIO.h"

using namespace std;
using namespace mio;

namespace mio {
/**
 * @page snowpack SNIO
 * @section snowpack_format Format
 * This is for reading grid data in the SNOWPACK meteo format
 *
 * @section snowpack_units Units
* - temperatures in celsius (input and output) or in kelvin (input only)
 * - relative humidity in % (input and output) or in [0;1] (input only)
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 *
 * @section snowpack_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords); [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords); [Input] section
 * - COORDSYS: output coordinate system (see Coords); [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords); [Output] section
 * - METEOPATH: path to the output directory; [Output] section
 * - METEOFILE#: input meteo data file, e.g. METEOFILE1, METEOFILE2; [Input] section
 * - STATION#: station name as listed in the METAFILE, e.g. STATION1, STATION2; [Input] section
 * - METAFILE: filename of the meta data file; [Input] section
 * - NROFSTATIONS: integer, the number of stations for which meteo files are provided; [Input] section
 */
}

const int SNIO::sn_julian_offset = 2415021;
const double SNIO::plugin_nodata = 0.0; //plugin specific nodata value

SNIO::SNIO(void (*delObj)(void*), const std::string& filename) : IOInterface(delObj), cfg(filename)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

SNIO::SNIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

SNIO::SNIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

SNIO::~SNIO() throw()
{
	cleanup();
}

void SNIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void SNIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*filename*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	//the meta data cannot change for the stations in dependence of time
	string strNrOfStations="";
	unsigned int nrOfStations = 0;

	cfg.getValue("NROFSTATIONS", "Input", strNrOfStations);
	if (!IOUtils::convertString(nrOfStations, strNrOfStations, std::dec))
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);

	if (vecAllStations.size() == 0)
		readMetaData(nrOfStations);
	
	vecStation = vecAllStations;
}

void SNIO::readMetaData(unsigned int& nrOfStations)
{	
	string stationname, metafile;
	cfg.getValue("METAFILE", "Input", metafile);
	if (!IOUtils::validFileName(metafile))
		throw InvalidFileNameException(metafile, AT);
	if (!IOUtils::fileExists(metafile))
		throw FileNotFoundException(metafile, AT);
	fin.clear();
	
	//Loop over all stations
	for (unsigned int ii=0; ii<nrOfStations; ii++){
		string line="";
		stringstream ss;
		ss << ii+1;

		cfg.getValue("STATION" + ss.str(), "Input", stationname);

		fin.open (metafile.c_str(), std::ifstream::in);	
		if (fin.fail())
			throw FileAccessException(metafile, AT);
		
		try{ 
			char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
			
			unsigned int linenr = 0;
			vector<string> tmpvec;

			while (!fin.eof()){
				getline(fin, line, eoln); //read complete line of data
				
				linenr++;
				ss.str("");
				ss << linenr;
			
				unsigned int ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)

				if (ncols==0){
					//Ignore empty lines
				} else if ((ncols<6) || (ncols>6)){
					throw InvalidFormatException(metafile+": each line must have 6 columns", AT);
				} else {

					//6 columns exist
					if (tmpvec.at(0) == stationname){
						StationData sd;
						parseMetaDataLine(tmpvec, sd);
						vecAllStations.push_back(sd);
					}
				}
			}
		} catch(std::exception& e){
			cleanup();
			throw;
		}
		cleanup();
	}
}

void SNIO::parseMetaDataLine(const vector<string>& vecLine, StationData& sd)
{
	if (vecLine.size() != 6)
		throw InvalidFormatException("While reading metadata: each line must have 6 columns", AT);	

	//Extract all data as double values
	vector<double> tmpdata = vector<double>(vecLine.size());
	for (unsigned int ii=2; ii<5; ii++) {
		if (!IOUtils::convertString(tmpdata[ii], vecLine[ii], std::dec))
			throw ConversionFailedException("While reading meta data for station " + vecLine[0], AT);
	}

	Coords stationcoord(coordin, coordinparam);
	stationcoord.setLatLon(tmpdata[3], tmpdata[4], tmpdata[2]);
	sd.setStationData(stationcoord, vecLine[0]);
}


void SNIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
					 std::vector< std::vector<MeteoData> >& vecMeteo, 
					 std::vector< std::vector<StationData> >& vecStation,
					 const unsigned int&)
{
	/*
	 * Read the meteorological snowpack input file, formatted as follows:
	 * M Date(ISO) Date(Julian) TA RH VW VDir ISWR RSWR ILWR TSS TSG HNW HS 
	 * The first line may be a comment, it won't start with "M", but with "MTO"
	 * The meteo data is terminated by a singular "END" on a line of its own
	 */

	vector<string> tmpvec;
	string strNrOfStations;
	unsigned int nrOfStations = 0;

	cfg.getValue("NROFSTATIONS", "Input", strNrOfStations);

	if (!IOUtils::convertString(nrOfStations, strNrOfStations, std::dec))
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);

	if (vecAllStations.size() == 0)
		readMetaData(nrOfStations);

	vecMeteo.clear();
	vecStation.clear();
	vecMeteo.insert(vecMeteo.begin(), vecAllStations.size(), vector<MeteoData>());
	vecStation.insert(vecStation.begin(), vecAllStations.size(), vector<StationData>());

	for (unsigned int ii=0; ii<vecAllStations.size(); ii++){
		string filename="", line="";
		stringstream ss;

		ss << ii+1;
		cfg.getValue("METEOFILE"+ss.str(), "Input", filename);

		if (!IOUtils::validFileName(filename))
			throw InvalidFileNameException(filename, AT);
		if (!IOUtils::fileExists(filename))
			throw FileNotFoundException(filename, AT);
  
		fin.clear();
		fin.open (filename.c_str(), std::ifstream::in);
	
		if (fin.fail())
			throw FileAccessException(filename, AT);
		if (fin.eof())
			throw InvalidFileNameException(filename + ": Empty file", AT);
	
		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
		
		try {
			getline(fin, line, eoln);      //read complete line meta information, ignore it
			if (line.length()>=3){
				if (line.substr(0,3) != "MTO") //if its not meta information rewind to the beginning
					fin.seekg (0, ios::beg);
			}else {
				throw InvalidFormatException(filename + ": first line in invalid format", AT);
			}
		

			unsigned int linenr = 0;

			while (!fin.eof()){
				getline(fin, line, eoln); //read complete line of data

				stringstream ss;
				linenr++;
				ss << linenr;

				unsigned int ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)
			
				if (ncols >= 15){//valid length for MeteoData
					MeteoData md;
					parseMeteoLine(tmpvec, filename + ":" + ss.str(), md);
					
					if ((md.date >= dateStart) && (md.date <= dateEnd)){//check date and add to vectors
						convertUnits(md);
						vecMeteo[ii].push_back(md);
						vecStation[ii].push_back(vecAllStations.at(ii));
					}
				} else if (ncols == 1){
					if (tmpvec.at(0) == "END") {
						break; //reached end of MeteoData
					} else {
						throw InvalidFormatException(filename + ":line " + ss.str() + " premature end of line", AT);
					}
				} else if (ncols == 0){
					//Ignore empty lines
				} else {
					throw InvalidFormatException(filename + ":line " + ss.str() + " premature end of line", AT);
				}
			}				
		} catch (std::exception& e){
			cleanup();
			throw;
		}	
		cleanup();
	}
}

void SNIO::parseMeteoLine(const vector<string>& vecLine, const string& filepos, MeteoData& md)
{
	if (vecLine.size() < 15)
		throw InvalidFormatException("At " + filepos + " line is too short", AT);

	if (vecLine[0] != "M")
		throw InvalidFormatException("At " + filepos + " meteo input lines must start with 'M'", AT);

	//deal with the date
	if (vecLine[1].length() != 10)
		throw InvalidFormatException("At " + filepos + " date format must be DD.MM.YYYY", AT);
	string year  = vecLine[1].substr(6,4);
	string month = vecLine[1].substr(3,2);
	string day   = vecLine[1].substr(0,2);
	
	if (!IOUtils::convertString(md.date, year+"-"+month+"-"+day+"T"+vecLine[2]))
		throw InvalidFormatException("At " + filepos + " date format invalid", AT);
	
	//Extract all data as double values
	vector<double> tmpdata = vector<double>(vecLine.size());
	for (unsigned int ii=4; ii<vecLine.size(); ii++) {
		if (!IOUtils::convertString(tmpdata[ii], vecLine[ii], std::dec))
			throw ConversionFailedException("At " + filepos, AT);
	}

	//Copy data into MeteoData object
	md.setData(MeteoData::TA, tmpdata[4]);
	md.setData(MeteoData::RH, tmpdata[5]);
	md.setData(MeteoData::VW, tmpdata[6]);
	md.setData(MeteoData::DW, tmpdata[7]);
	md.setData(MeteoData::ISWR, tmpdata[8]);
	md.setData(MeteoData::RSWR, tmpdata[9]);
	md.setData(MeteoData::ILWR, tmpdata[10]);
	md.setData(MeteoData::TSS, tmpdata[11]);
	md.setData(MeteoData::TSG, tmpdata[12]);
	md.setData(MeteoData::HNW, tmpdata[13]);
	md.setData(MeteoData::HS, tmpdata[14]);
}

void SNIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
                          const std::vector< std::vector<StationData> >& vecStation,
                          const std::string&)
{
	string path="";
	cfg.getValue("METEOPATH", "Output", path);

	for(unsigned int ii=0; ii<vecMeteo.size(); ii++) {
		if(vecStation[ii].size()>0) {
			std::string station_name = vecStation[ii][0].stationName;
			if (station_name == "") station_name = "UNKNOWN";
			const std::string output_name = path + "/" + station_name + ".inp";
			if( !IOUtils::fileExists(output_name) ) {
				fout.open(output_name.c_str());
				writeStationHeader(vecMeteo[ii], station_name);
			} else {
				fout.open(output_name.c_str());
			}
			writeStationMeteo(vecMeteo[ii], output_name);
			fout.close();
		}
	}
}

void SNIO::writeStationHeader(const std::vector<MeteoData>& Meteo, const std::string station_name)
{
	//writing the (very basic) metadata
	fout << "MTO <" << station_name << "> " << Meteo.size() << "\n";
}

void SNIO::writeStationMeteo(const std::vector<MeteoData>& Meteo, const std::string& file_name)
{ //write out the data for 1 station
	unsigned int failure_count = 0;
	unsigned int Dirichlet_failure_count = 0;

	for(unsigned int ii=0; ii<Meteo.size(); ii++) {
		int YYYY, MM, DD, HH, MI;
		Meteo[ii].date.getDate(YYYY, MM, DD, HH, MI);
		const double sn_julian = Meteo[ii].date.getJulianDate() - sn_julian_offset + 0.5;
		const double ta = Meteo[ii].ta;
		const double rh = Meteo[ii].rh;
		const double hnw = Meteo[ii].hnw;
		const double vw = Meteo[ii].vw;
		const double dw = Meteo[ii].dw;
		const double iswr = Meteo[ii].iswr;
		const double rswr = Meteo[ii].rswr;
		const double ilwr = Meteo[ii].ilwr;
		const double tss = Meteo[ii].tss;
		const double tsg = Meteo[ii].tsg;
		const double hs = Meteo[ii].hs;

		fout.fill('0');
		fout << "M " << setw(2) << DD << "." << setw(2) << MM << "." << setw(4) << YYYY << " " << setw(2) << HH << ":" << setw(2) << MI << " ";
		fout.flags ( ios::fixed );
		fout << setprecision(6) << setw(12) << sn_julian << " ";

		//default formatting parameters for the measurements
		fout.flags ( ios::fixed );
		fout.fill(' ');
		fout.width(6);

		//TA, RH, VW, DW
		if(ta==IOUtils::nodata) {
			failure_count++;
			fout << setw(6) << setprecision(0) << ta << " ";
		} else
			fout << setw(6) << setprecision(2) << K_TO_C(ta) << " ";
		if(rh==IOUtils::nodata) {
			failure_count++;
			fout << setw(5) << setprecision(0) << rh << " ";
		} else
			fout << setw(5) << setprecision(1) << rh * 100. << " ";
		if(vw==IOUtils::nodata) {
			failure_count++;
			fout << setw(4) << setprecision(0) << vw << " ";
		} else {
			fout << setw(4) << setprecision(1) << vw << " ";
		}
		if(dw==IOUtils::nodata)
			failure_count++;
		fout << setw(4) << setprecision(0) << dw << " ";

		//ISWR, RSWR
		if(iswr==IOUtils::nodata && rswr==IOUtils::nodata) {
			failure_count++;
			fout << setw(4) << setprecision(0) << iswr << " " << setprecision(0) << rswr << " ";
		} else {
			if(iswr==IOUtils::nodata)
				fout << setw(4) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(4) << setprecision(0) << iswr << " ";
			if(rswr==IOUtils::nodata)
				fout << setw(4) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(4) << setprecision(0) << rswr << " ";
		}

		//LWR
		if(ilwr==IOUtils::nodata) {
			failure_count++;
			fout << setw(4) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(4) << setprecision(0) << ilwr << " ";
		}

		//TSS, TSG (only required for Dirichlet)
		if(tss==IOUtils::nodata) {
			Dirichlet_failure_count++;
			fout << setw(6) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(6) << setprecision(2) << K_TO_C(tss) << " ";
		}
		if(tsg==IOUtils::nodata) {
			Dirichlet_failure_count++;
			fout << setw(6) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(6) << setprecision(2) << K_TO_C(tsg) << " ";
		}

		//HNW, HS
		if(hnw==IOUtils::nodata && hs==IOUtils::nodata) {
			failure_count++;
			fout << setw(5) << setprecision(2) << hnw << " " << setprecision(3) << hs << " ";
		} else {
			if(hnw==IOUtils::nodata)
				fout << setw(5) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(5) << setprecision(2) << hnw << " ";
			if(hs==IOUtils::nodata)
				fout << setw(5) << setprecision(1) << "0.0";
			else
				fout << setw(5) << setprecision(3) << hs;
		}

		//we don't write any snow depth temperatures.
		//we can not write wind velocity at the wind station, but since it is optional...

		fout << endl;
	}

	fout << "END" << endl;

	if(failure_count>0 || Dirichlet_failure_count>0) {
		std::cout << "[W] " << failure_count << " (and potentially " << Dirichlet_failure_count <<
		" more) errors found when writing " << file_name << std::endl;
	}
}


void SNIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

/*void SNIO::convertUnitsBack(MeteoData& meteo)
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
}
*/
void SNIO::convertUnits(MeteoData& meteo)
{
	//converts C to Kelvin, converts ilwr to ea, converts RH to [0,1]
	if(meteo.ta!=IOUtils::nodata){
		if (meteo.ta < 100)
			meteo.ta=C_TO_K(meteo.ta);
	}
	
	if(meteo.tsg!=IOUtils::nodata){
		if (meteo.tsg < 100)
			meteo.tsg=C_TO_K(meteo.tsg);
	}
	
	if(meteo.tss!=IOUtils::nodata){
		if (meteo.tss < 100)
			meteo.tss=C_TO_K(meteo.tss);
	}

	if (meteo.rh!=IOUtils::nodata){
		if (meteo.rh>1.2)
			meteo.rh /= 100;
	}
}


#ifndef _METEOIO_JNI
extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const string& classname, const string& filename) {
		if(classname == "SNIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new SNIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif
