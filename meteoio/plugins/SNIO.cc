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
#include <meteoio/meteolaws/Atmosphere.h>

using namespace std;

namespace mio {
/**
 * @page snowpack SNIO
 * @section snowpack_format Format
 * This is for reading meteo data in the SNOWPACK meteo format. The metadata has to be provided
 * in a separate file that might contain multiple stations, one per line. Each line has the following structure:\n
 * ALI2 Allieres:Chenau 1767 6.993 46.489 1.22 \n
 * where the first field is the short name, followed by the fullname and the location, then the elevation,
 * the longitude, the latitude and a wind coefficient (unused by MeteoIO). The short name is used for
 * identifying the station and matching it with the data file (name given in io.ini). If no such metadata file is
 * provided, the metadata will be left nodata. This only makes sense if the metadata would be later filled by another way
 * (like a merge).
 *
 * @section snowpack_units Units
 * - temperatures in degrees Celsius (input and output) or in kelvins (input only)
 * - relative humidity in % (input and output) or in [0;1] (input only)
 * - wind speed in m/s
 * - precipitations in mm/h (kg/m²)
 * - radiation in W/m²
 *
 * @section snowpack_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - METEOPATH: path to the meteo files directory; [Input] and [Output] sections
 * - METEOFILE#: input meteo data file, e.g. METEOFILE1, METEOFILE2; [Input] section
 * - STATION#: station name as listed in the METAFILE, e.g. STATION1, STATION2; [Input] section
 * - METAFILE: filename of the meta data file (in METEOPATH); [Input] section (optional but recommended)
 * - NROFSTATIONS: integer, the number of stations for which meteo files are provided; [Input] section
 * - optional, data must follow order given below but may be missing:
 * 	- NUMBER_MEAS_TEMPERATURES: integer, the number of measured snow temperatures provided; [Input] section \n
 * 		The depths of the sensors can be given under FIXED_SENSOR_DEPTHS (default: 0.25, 0.5, 1.0, 1.5, -0.1 m)
 * 	- NUMBER_OF_SOLUTES: integer, the number of solutes for which input data are provided; [Input] section
 * 	- VW_DRIFT: bool, a wind velocity to use for blowing and drifting snow is provided; [Input] section
 * 	- RHO_HN: bool, measured new snow density is provided; [Input] section
 */

const int SNIO::sn_julian_offset = 2415021;
const double SNIO::plugin_nodata = -999.0; //plugin specific nodata value
const unsigned int SNIO::min_nr_meteoData = 15;

SNIO::SNIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = out_tz = 0.;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_tz,Config::nothrow);
}

SNIO::SNIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = out_tz = 0.;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_tz,Config::nothrow);
}

SNIO::SNIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = out_tz = 0.;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_tz,Config::nothrow);
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
	vecStation.clear();

	cfg.getValue("NROFSTATIONS", "Input", strNrOfStations);
	if (!IOUtils::convertString(nrOfStations, strNrOfStations, std::dec))
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);

	if (vecAllStations.size() == 0)
		readMetaData(nrOfStations);
	
	vecStation = vecAllStations; //vecAllStations is a global vector that holds all meta data
}

bool SNIO::readStationMetaData(const std::string& metafile, const std::string& stationname, StationData& sd)
{
	fin.open (metafile.c_str(), std::ifstream::in);
	if (fin.fail())
		throw FileAccessException(metafile, AT);

	try{
		string line="";
		const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		unsigned int linenr = 0;
		vector<string> tmpvec;

		while (!fin.eof()){
			getline(fin, line, eoln); //read complete line of data

			linenr++;
			stringstream ss;
			ss << linenr;

			unsigned int ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)

			if (ncols==0){
				//Ignore empty lines
			} else if ((ncols<6) || (ncols>6)){
				throw InvalidFormatException(metafile+":"+ss.str() + " each line must have 6 columns", AT);
			} else {
				//6 columns exist
				if (tmpvec.at(0) == stationname){
					parseMetaDataLine(tmpvec, sd);
					return(true);
				}
			}
		}
		return(false);
	} catch(std::exception& e){
		cleanup();
		throw;
	}
}

void SNIO::readMetaData(unsigned int& nrOfStations)
{
	string stationID, metafile="", inpath;
	cfg.getValue("METAFILE", "Input", metafile, Config::nothrow);
	cfg.getValue("METEOPATH", "Input", inpath);

	fin.clear();
	
	//Loop over all stations
	for (unsigned int ii=0; ii<nrOfStations; ii++){
		stringstream snum;
		snum << ii+1;

		cfg.getValue("STATION" + snum.str(), "Input", stationID);

		StationData sd;
		if(metafile!="") { //a metafile has been provided, so get metadata
			stringstream meta_with_path;
			meta_with_path << inpath << "/" << metafile;
			if ( !IOUtils::validFileName(meta_with_path.str()) )
				throw InvalidFileNameException(meta_with_path.str(), AT);
			if ( !IOUtils::fileExists(meta_with_path.str()) )
				throw FileNotFoundException(meta_with_path.str(), AT);
			if( readStationMetaData(meta_with_path.str(), stationID, sd)==false) {
				stringstream ss;
				ss << "No metadata found for station " << stationID << " in " << metafile;
				throw NoAvailableDataException(ss.str(), AT);
			}
		}
		vecAllStations.push_back(sd);
		cleanup();
	}
}

void SNIO::parseMetaDataLine(const std::vector<std::string>& vecLine, StationData& sd)
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
	sd.setStationData(stationcoord, vecLine[0], vecLine[1]);
}


void SNIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                         std::vector< std::vector<MeteoData> >& vecMeteo, const unsigned int&)
{
	/*
	 * Read the meteorological snowpack input file, formatted as follows:
	 * M Date Time Date(Julian) TA RH VW DW ISWR RSWR ILWR TSS TSG HNW HS (TS1 TS2 ... TSN CONC0 CONC1 ... CONCM rho_hn)
	 * The first line may be a comment, it won't start with "M", but with "MTO"
	 * The meteo data is terminated by a singular "END" on a line of its own
	 */

	vector<string> tmpvec;
	string strNrOfStations, inpath;
	unsigned int nrOfStations = 0;

	cfg.getValue("NROFSTATIONS", "Input", strNrOfStations);
	cfg.getValue("METEOPATH", "Input", inpath);

	if (!IOUtils::convertString(nrOfStations, strNrOfStations, std::dec))
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);

	if (vecAllStations.size() == 0)
		readMetaData(nrOfStations);

	vecMeteo.clear();
	vecMeteo.insert(vecMeteo.begin(), vecAllStations.size(), vector<MeteoData>());

	for (unsigned int ii=0; ii<vecAllStations.size(); ii++){
		string filename="", line="";
		stringstream ss, file_with_path;

		ss << ii+1;
		cfg.getValue("METEOFILE"+ss.str(), "Input", filename);
		file_with_path << inpath << "/" << filename;

		if ( !IOUtils::validFileName(file_with_path.str()) )
			throw InvalidFileNameException(file_with_path.str(), AT);
		if ( !IOUtils::fileExists(file_with_path.str()) )
			throw FileNotFoundException(file_with_path.str(), AT);
  
		fin.clear();
		fin.open (file_with_path.str().c_str(), std::ifstream::in);

		if (fin.fail())
			throw FileAccessException(file_with_path.str(), AT);
		if (fin.eof())
			throw InvalidFileNameException(file_with_path.str() + ": Empty file", AT);
	
		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		try {
			getline(fin, line, eoln);      //read complete line meta information, ignore it
			if (line.length()>=3){
				if (line.substr(0,3) != "MTO") //if its not meta information rewind to the beginning
					fin.seekg (0, ios::beg);
			}else {
				throw InvalidFormatException(file_with_path.str() + ": first line in invalid format", AT);
			}

			unsigned int linenr = 0;

			while (!fin.eof()){
				getline(fin, line, eoln); //read complete line of data

				stringstream ss;
				linenr++;
				ss << linenr;

				unsigned int ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)
			
				if (ncols >= min_nr_meteoData){
					MeteoData md;
					md.meta = vecAllStations[ii];
					parseMeteoLine(tmpvec, file_with_path.str() + ":" + ss.str(), md);
					
					if ((md.date >= dateStart) && (md.date <= dateEnd)){//check date and add to vectors
						convertUnits(md);
						vecMeteo[ii].push_back(md);
					}
				} else if (ncols == 1){
					if (tmpvec.at(0) == "END") {
						break; //reached end of MeteoData
					} else {
						throw InvalidFormatException(file_with_path.str() + ":line " + ss.str() + " premature end of line", AT);
					}
				} else if (ncols == 0){
					//Ignore empty lines
				} else {
					throw InvalidFormatException(file_with_path.str() + ":line " + ss.str() + " premature end of line", AT);
				}
			}				
		} catch (std::exception& e){
			cleanup();
			throw;
		}
		cleanup();
	}
}

void SNIO::parseMeteoLine(const std::vector<std::string>& vecLine, const std::string& filepos, MeteoData& md)
{
	/*
	 * This function takes a meteo line, extracts the date (ignores Julian) and then converts
	 * all meteo parameters to doubles and finally copies them into the MeteoData object md
	 */
	std::string err_msg("");
	try {
		if (vecLine.size() < min_nr_meteoData)
			throw err_msg = " line is too short\n";

		if (vecLine[0] != "M")
			throw err_msg = " meteo input lines must start with 'M'\n";

		//deal with the date
		if (vecLine[1].length() != 10)
			throw err_msg = " date format must be DD.MM.YYYY\n";
		const string year  = vecLine[1].substr(6,4);
		const string month = vecLine[1].substr(3,2);
		const string day   = vecLine[1].substr(0,2);
	
		if (!IOUtils::convertString(md.date, year+"-"+month+"-"+day+"T"+vecLine[2], in_tz, std::dec))
			throw err_msg = " date format invalid\n";
	
		//Extract all data as double values
		vector<double> tmpdata = vector<double>(vecLine.size());
		for (unsigned int ii=4; ii<vecLine.size(); ii++) {
			if (!IOUtils::convertString(tmpdata[ii], vecLine[ii], std::dec))
				throw err_msg = " Conversion failed!\n";
		}

		//Copy data into MeteoData object
		md.setData(MeteoData::TA, tmpdata[4]);
		md.setData(MeteoData::RH, tmpdata[5]);
		md.setData(MeteoData::VW, tmpdata[6]);
		md.setData(MeteoData::DW, tmpdata[7]);
		md.setData(MeteoData::ISWR, tmpdata[8]);
		md.setData(MeteoData::RSWR, tmpdata[9]);
	
	double& ea = tmpdata[10];
	if ((ea <= 1) && (ea != plugin_nodata)){
		if ((md.ta != plugin_nodata) && (md.rh != plugin_nodata)) {
			if(ea==0.)
				ea = Atmosphere::Brutsaert_ilwr(md.rh/100., C_TO_K(md.ta));
			else
				ea = Atmosphere::Omstedt_ilwr(md.rh/100., C_TO_K(md.ta), ea); //calculate ILWR from cloudiness
		} else {
			ea = plugin_nodata;
		}
	}

	md.setData(MeteoData::ILWR, ea);
		md.setData(MeteoData::TSS, tmpdata[11]);
		md.setData(MeteoData::TSG, tmpdata[12]);
		md.setData(MeteoData::HNW, tmpdata[13]);
		md.setData(MeteoData::HS, tmpdata[14]);

	// Read optional values
	 	// TS[]: snow temperatures
		unsigned int ii = min_nr_meteoData, jj;
		unsigned int number_meas_temperatures = 0;
		cfg.getValue("NUMBER_MEAS_TEMPERATURES", "Input", number_meas_temperatures, Config::nothrow);
		if (vecLine.size() < min_nr_meteoData + number_meas_temperatures)
			throw err_msg = ": Not enough measured temperatures data\n";
		stringstream ss("");
		for (jj=1; jj<=number_meas_temperatures; jj++) {
			ss.str("");
			ss << "TS" << (jj);
			md.addParameter(ss.str());
			md.param(ss.str()) = tmpdata[ii++];
		}
		// CONC[]: solute concentrations
		unsigned int number_of_solutes = 0;
		cfg.getValue("NUMBER_OF_SOLUTES", "Input", number_of_solutes, Config::nothrow);
		if (vecLine.size() < min_nr_meteoData + number_meas_temperatures + number_of_solutes)
			throw err_msg = ": Not enough solute data\n";
		jj = 0;
		for (jj = 0 ; jj<number_of_solutes; jj++) {
			ss.str("");
			ss << "CONC" << jj;
			md.addParameter(ss.str());
			md.param(ss.str()) = tmpdata[ii++];
		}
		// VW_DRIFT: optional wind velocity for blowing and drifting snow
		bool vw_drift = false;
		cfg.getValue("VW_DRIFT", "Input", vw_drift, Config::nothrow);
		if (vw_drift) {
			if (vecLine.size() < ii)
				throw err_msg = ": No data for vw_drift\n";
			md.addParameter("VW_DRIFT");
			md.param("VW_DRIFT") = tmpdata[ii++];
		}
		// RHO_HN: measured new snow density
		bool rho_hn = false;
		cfg.getValue("RHO_HN", "Input", rho_hn, Config::nothrow);
		if (rho_hn) {
			if (vecLine.size() < ii)
				throw err_msg = ": No data for rho_hn\n";
			md.addParameter("RHO_HN");
			md.param("RHO_HN") = tmpdata[ii++];
		}
		if (vecLine.size() > ii)
			throw err_msg = ": Too many data\n";
	} catch(string& str) {
		throw InvalidFormatException("At " + filepos + str, AT);
	}
}

void SNIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	string outpath="";
	cfg.getValue("METEOPATH", "Output", outpath);

	for(unsigned int ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].size() > 0) {
			std::string station_id = vecMeteo[ii][0].meta.getStationID();
			if (station_id == "") station_id = "UNKNOWN";
			const std::string output_name = outpath + "/" + station_id + ".inp";
			if( !IOUtils::fileExists(output_name) ) {
				fout.open(output_name.c_str());
				writeStationHeader(vecMeteo[ii], station_id);
			} else {
				fout.open(output_name.c_str());
			}
			writeStationMeteo(vecMeteo[ii], output_name);
			fout.close();
		}
	}
}

void SNIO::writeStationHeader(const std::vector<MeteoData>& vecmd, const std::string station_id)
{
	//writing the (very basic) metadata
	fout << "MTO <" << station_id << "> " << vecmd.size() << "\n";
}

void SNIO::writeStationMeteo(const std::vector<MeteoData>& vecmd, const std::string& file_name)
{ //write out the data for 1 station
	unsigned int failure_count = 0;
	unsigned int Dirichlet_failure_count = 0;
	unsigned int optional_failure_count = 0;

	for(unsigned int jj=0; jj<vecmd.size(); jj++) {
		int YYYY, MM, DD, HH, MI;
		Date tmp_date(vecmd[jj].date);
		tmp_date.setTimeZone(out_tz);
		tmp_date.getDate(YYYY, MM, DD, HH, MI);
		const double sn_julian = tmp_date.getJulianDate() - sn_julian_offset + 0.5;
		const double ta = vecmd[jj].ta;
		const double rh = vecmd[jj].rh;
		const double hnw = vecmd[jj].hnw;
		const double vw = vecmd[jj].vw;
		const double dw = vecmd[jj].dw;
		const double iswr = vecmd[jj].iswr;
		const double rswr = vecmd[jj].rswr;
		const double ilwr = vecmd[jj].ilwr;
		const double tss = vecmd[jj].tss;
		const double tsg = vecmd[jj].tsg;
		const double hs = vecmd[jj].hs;

		fout.fill('0');
		fout << "M " << setw(2) << DD << "." << setw(2) << MM << "." << setw(4) << YYYY << " " << setw(2) << HH << ":" << setw(2) << MI << " ";
		fout.flags ( ios::fixed );
		fout << setprecision(6) << setw(12) << sn_julian << " ";

		//default formatting parameters for the measurements
		fout.flags(ios::fixed);
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
			fout << setw(6) << setprecision(0) << iswr << " " << setprecision(0) << rswr << " ";
		} else {
			if(iswr==IOUtils::nodata)
				fout << setw(6) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(6) << setprecision(1) << iswr << " ";
			if(rswr==IOUtils::nodata)
				fout << setw(6) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(6) << setprecision(1) << rswr << " ";
		}

		//LWR
		if(ilwr==IOUtils::nodata) {
			failure_count++;
			fout << setw(5) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(5) << setprecision(1) << ilwr << " ";
		}

		//TSS, TSG (only required for Dirichlet)
		if(tss==IOUtils::nodata) {
			Dirichlet_failure_count++;
			fout << setw(7) << setprecision(1) << "0.0" << " ";
		} else {
			fout << setw(7) << setprecision(2) << K_TO_C(tss) << " ";
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
			fout << setw(7) << setprecision(2) << hnw << " " << setw(6) << setprecision(3) << hs << " ";
		} else {
			if(hnw==IOUtils::nodata)
				fout << setw(7) << setprecision(1) << "0.0" << " ";
			else
				fout << setw(7) << setprecision(4) << hnw << " ";
			if(hs==IOUtils::nodata)
				fout << setw(6) << setprecision(1) << "0.0";
			else
				fout << setw(6) << setprecision(3) << hs;
		}

		// Write optional values
		//TS[]: snow temperatures
		stringstream ss;
		for (unsigned int kk=1; kk<100; kk++) {
			ss.str("");
			ss << "TS" << kk;
			if (vecmd[jj].param_exists(ss.str())){
				const double ts = vecmd[jj].param(ss.str());
				if (ts == IOUtils::nodata) {
					optional_failure_count++;
					fout << setw(7) << setprecision(0) << ts << " ";
				} else {
					fout << setw(7) << setprecision(2) << K_TO_C(ts) << " ";
				}
			} else {
				break;
			}
		}
		//CONC[]: solute concentrations
		for (unsigned int kk=0; kk<100; kk++) {
			ss.str("");
			ss << "CONC" << kk;
			if (vecmd[jj].param_exists(ss.str())) {
				const double conc = vecmd[jj].param(ss.str());
				if (conc == IOUtils::nodata) {
					optional_failure_count++;
					fout << setw(6) << setprecision(0) << conc << " ";
				} else {
					fout << setw(6) << setprecision(4) << conc << " ";
				}
			} else {
				break;
			}
		}
		// VW_DRIFT: optional wind velocity for blowing and drifting snow
		if (vecmd[jj].param_exists("VW_DRIFT")) {
			const double vw_drift = vecmd[jj].param("VW_DRIFT");
			if (vw_drift == IOUtils::nodata) {
				optional_failure_count++;
				fout << setw(4) << setprecision(0) << vw_drift << " ";
			} else {
				fout << setw(4) << setprecision(1) << vw_drift << " ";
			}
		}
		// RHO_HN: measured new snow density
		if (vecmd[jj].param_exists("RHO_HN")) {
			const double rho_hn = vecmd[jj].param("RHO_HN");
			if (rho_hn == IOUtils::nodata) {
				optional_failure_count++;
				fout << setw(6) << setprecision(0) << rho_hn << " ";
			} else {
				fout << setw(6) << setprecision(1) << rho_hn << " ";
			}
		}

		fout << endl;
	}

	fout << "END" << endl;

	if ((failure_count > 0) || (Dirichlet_failure_count > 0) || (optional_failure_count > 0)) {
		std::cout << "[W] " << failure_count << " basic input data, " << Dirichlet_failure_count <<
				" Dirichlet boundary condition data, and " << optional_failure_count <<
				" optional data found missing when writing " << file_name << std::endl;
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

	stringstream ss;
	for (unsigned int ii=1; ii<50; ii++){
		ss.str("");
		ss << "TS" << ii;
		if (meteo.param_exists(ss.str())){
			double& value = meteo.param(ss.str());
			value = C_TO_K(value);
		} else {
			break;
		}
	}
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
		if(classname == "SNIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new SNIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
