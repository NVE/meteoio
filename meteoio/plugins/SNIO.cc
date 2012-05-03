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
 * This is for reading meteo data in the SNOWPACK meteo format. The metadata might be provided
 * in a separate file that can contain multiple stations, one per line. Each line has the following structure:
 * - ALI2 Allieres:Chenau 1767 6.993 46.489 1.22 \n
 *
 * where the first field is the short name, followed by the fullname and the location, then the elevation,
 * the longitude, the latitude and a wind coefficient (unused by MeteoIO). The short name is used for
 * identifying the station (stationID) and matching it with the data file (name given in io.ini).
 * If no such metadata file is provided, the metadata will be left nodata. This only makes sense
 * if the metadata would be later filled by another way (like a merge).
 *
 * In any case, the header line of the data file MUST contain a short name suitable to be used as station_id.
 *
 * Finally, when writing to a file, the header line will only be created if the file does not already exist. If the
 * file already exists on the disk, new data will be appended without attempting to write a header line.
 *
 * @section snowpack_units Units
 * - temperatures in degrees Celsius (input and output) or in kelvins (input only)
 * - relative humidity in % (input and output) or in [0;1] (input only)
 * - wind speed in m/s
 * - precipitations in mm w.e. (kg/m²) per meteo time step
 * - radiation in W/m²
 *
 * @section snowpack_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - METEOPATH: path to the meteo files directory; [Input] and [Output] sections
 * - STATION#: input meteo data file, e.g. STATION1, STATION2; [Input] section
 * - METAFILE: filename of the meta data file (in METEOPATH); [Input] section (optional but recommended)
 * - optional:
 * 	- ISWR_INP or RSWR_INP: if one of these data is missing, set corresponding switch to false. The setting
 *                        will be valid for all stations.
 * 	- additional data must follow order given below but may be missing:
 * 		- NUMBER_MEAS_TEMPERATURES: integer, the number of measured snow temperatures provided; [Input] section \n
 * 		- NUMBER_OF_SOLUTES: integer, the number of solutes for which input data are provided; [Input] section
 * 		- VW_DRIFT: bool, a wind velocity to use for blowing and drifting snow is provided; [Input] section
 * 		- RHO_HN: bool, measured new snow density is provided; [Input] section
 *
 * @section snowpack_errors Errors
 * When writing in Snowpack format, potential errors in the data set are written out. The error count is split between the different
 * types of errors:
 *	- basic input data: each mandatory parameter that is missing at a timestep increments this counter;
 *	- Dirichlet boundary condition data: each TSS or TSG parameter that is missing at a timestep increments this counter;
 *	- optional data: each snow temperature, solutes, snow density, wind drift optional parameter that is missing (if it
 *	 was previously written out) increments this counter;
 * Overall, this means that the count of basic errors should be zero while Dirichlet errors might be tolerable as well as optional data errors.
 */

const int SNIO::sn_julian_offset = 2415021;
const double SNIO::plugin_nodata = -999.0; //plugin specific nodata value
const size_t SNIO::min_nr_meteoData = 15;
const std::string SNIO::dflt_extension = ".inp";

SNIO::SNIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = out_tz = 0.;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_tz,Config::nothrow);
	iswr_inp = rswr_inp = true;
	cfg.getValue("ISWR_INP","Input",iswr_inp,Config::nothrow);
	cfg.getValue("RSWR_INP","Input",rswr_inp,Config::nothrow);
	nr_meteoData = min_nr_meteoData;
	if (!iswr_inp || !rswr_inp)
		nr_meteoData = min_nr_meteoData - 1;
}

SNIO::SNIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = out_tz = 0.;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_tz,Config::nothrow);
	iswr_inp = rswr_inp = true;
	cfg.getValue("ISWR_INP","Input",iswr_inp,Config::nothrow);
	cfg.getValue("RSWR_INP","Input",rswr_inp,Config::nothrow);
	nr_meteoData = min_nr_meteoData;
	if (!iswr_inp || !rswr_inp)
		nr_meteoData = min_nr_meteoData - 1;
}

SNIO::SNIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	in_tz = out_tz = 0.;
	cfg.getValue("TIME_ZONE","Input",in_tz,Config::nothrow);
	cfg.getValue("TIME_ZONE","Output",out_tz,Config::nothrow);
	iswr_inp = rswr_inp = true;
	cfg.getValue("ISWR_INP","Input",iswr_inp,Config::nothrow);
	cfg.getValue("RSWR_INP","Input",rswr_inp,Config::nothrow);
	nr_meteoData = min_nr_meteoData;
	if (!iswr_inp || !rswr_inp)
		nr_meteoData = min_nr_meteoData - 1;
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

void SNIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
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
	if (vecAllStations.size() == 0)
		readMetaData();

	vecStation = vecAllStations; //vecAllStations is a global vector that holds all meta data
}

bool SNIO::readStationMetaData(const std::string& metafile, const std::string& stationID, StationData& sd)
{
	if (!IOUtils::validFileName(metafile)) {
		std::stringstream ss;
		ss << "\"" << metafile << "\" is not a valid file name. Please check your METAFILE key!";
		throw InvalidFileNameException(metafile, AT);
	}

	if (!IOUtils::fileExists(metafile)) {
		std::stringstream ss;
		ss << "File \"" << metafile << "\" does not exist. Please check your METAFILE key!";
		throw FileNotFoundException(ss.str(), AT);
	}


	fin.clear();
	fin.open (metafile.c_str(), std::ifstream::in);
	if (fin.fail())
		throw FileAccessException(metafile, AT);

	try {
		string line="";
		const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		size_t linenr = 0;
		vector<string> tmpvec;

		while (!fin.eof()) {
			linenr++;
			getline(fin, line, eoln); //read complete line of data

			const size_t ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)

			if (ncols==0) {
				//Ignore empty lines
			} else if (ncols == 6) { //valid line, e.g. "MST96 Weissfluhjoch:StudyPlot_MST 2540 9.81 46.831 1.00"
				if (tmpvec.at(0) == stationID) {
					parseMetaDataLine(tmpvec, sd);
					cleanup();
					return true;
				}
			} else {
				stringstream ss;
				ss << linenr;
				throw InvalidFormatException(metafile+":"+ss.str() + " each line must have 6 columns", AT);
			}
		}
		cleanup();
		return false;
	} catch(const std::exception&){
		cleanup();
		throw;
	}
}

std::string SNIO::getStationID(const std::string& filename)
{
	/**
	 * This function will return the station name as retrieved from
	 * the first line of a SNIO formatted meteo file
	 */
	if ( !IOUtils::validFileName(filename) )
		throw InvalidFileNameException(filename, AT);
	if ( !IOUtils::fileExists(filename) )
		throw FileNotFoundException(filename, AT);

	fin.clear();
	fin.open (filename.c_str(), std::ifstream::in);

	if (fin.fail())
		throw FileAccessException(filename, AT);
	if (fin.eof())
		throw InvalidFileNameException(filename + ": Empty file", AT);

	const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	string station_id = "";
	try {
		string line = "";
		vector<string> tmpvec;

		getline(fin, line, eoln);      //read complete line meta information, parse it
		const size_t ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)
		if ((ncols != 3) || (tmpvec.at(0) != "MTO") || (tmpvec.at(1).size() < 3) || (tmpvec.at(1)[0] != '<') || (tmpvec.at(1)[tmpvec.at(1).length()-1] != '>'))
			throw InvalidFormatException(filename + ": first line in invalid format", AT);

		//Now get the 2nd column looking something like <{STATIONNAME}Data>
		station_id = tmpvec[1].substr(1, tmpvec[1].length()-1); //leaving away the
		size_t pos = station_id.find("Data");
		if (pos != string::npos) {
			station_id = station_id.substr(0, pos);
		} else {
			throw InvalidFormatException(filename + ": first line in invalid format", AT);
		}
	} catch (...) {
		cleanup();
		throw;
	}

	cleanup();
	return station_id;
}

void SNIO::readMetaData()
{
	/**
	 * Parse through the io.ini file and read the desired file names STATION#
	 */
	vecAllStations.clear();
	string metafile="", inpath="";
	cfg.getValue("METAFILE", "Input", metafile, Config::nothrow);
	cfg.getValue("METEOPATH", "Input", inpath);

	size_t current_stationnr = 1;
	string current_station;
	do { //Look for STATION1, STATION2, etc
		current_station = "";
		stringstream ss;
		ss << "STATION" << current_stationnr;
		cfg.getValue(ss.str(), "Input", current_station, Config::nothrow);

		if ((current_stationnr == 1) && (current_station == ""))
			throw InvalidFormatException("Missing key \"STATION1\" in config: Please specify a SNOWPACK formatted meteo data file", AT);

		if (current_station != ""){
			if(IOUtils::getExtension(current_station)=="") current_station += dflt_extension; //default extension
			string station_id = getStationID(inpath+ "/" +current_station);

			StationData sd(Coords(), station_id);
			if (metafile!="") { //a metafile has been provided, so get metadata
				if (readStationMetaData(inpath+ "/" +metafile, station_id, sd) == false) {
					stringstream msg;
					msg << "No metadata found for station " << station_id << " in " << metafile;
					throw NoAvailableDataException(msg.str(), AT);
				}
			}
			vecAllStations.push_back(sd);
			cerr << "\t[i] Read meta data for station ID \"" << station_id << "\"\n";
		}

		current_stationnr++;
	} while (current_station != "");
}

void SNIO::parseMetaDataLine(const std::vector<std::string>& vecLine, StationData& sd)
{
	if (vecLine.size() != 6)
		throw InvalidFormatException("While reading metadata: each line must have 6 columns", AT);

	//Extract all data as double values
	vector<double> tmpdata = vector<double>(vecLine.size());
	for (size_t ii=2; ii<5; ii++) {
		if (!IOUtils::convertString(tmpdata[ii], vecLine[ii], std::dec))
			throw ConversionFailedException("While reading meta data for station " + vecLine[0], AT);
	}

	Coords stationcoord(coordin, coordinparam);
	stationcoord.setLatLon(tmpdata[3], tmpdata[4], tmpdata[2]);
	sd.setStationData(stationcoord, vecLine[0], vecLine[1]);
}


void SNIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                         std::vector< std::vector<MeteoData> >& vecMeteo, const size_t&)
{
	/*
	 * Read the meteorological snowpack input file, formatted as follows:
	 * M Date Time Date(Julian) TA RH VW DW ISWR RSWR ILWR TSS TSG HNW HS (TS1 TS2 ... TSN CONC0 CONC1 ... CONCM rho_hn)
	 * The first line may be a comment, it won't start with "M", but with "MTO"
	 * The meteo data is terminated by a singular "END" on a line of its own
	 */
	vector<string> tmpvec;
	string inpath = "";

	cfg.getValue("METEOPATH", "Input", inpath);

	if (vecAllStations.size() == 0)
		readMetaData();

	vecMeteo.clear();
	vecMeteo.insert(vecMeteo.begin(), vecAllStations.size(), vector<MeteoData>());
	if (vec_streampos.size() == 0) //the vec_streampos save file pointers for certain dates
		vec_streampos = vector< map<Date, std::streampos> >(vecAllStations.size());

	for (size_t ii=0; ii<vecAllStations.size(); ii++){
		string filename="", line="";
		stringstream ss, file_with_path;

		ss << ii+1;
		cfg.getValue("STATION"+ss.str(), "Input", filename);
		if(IOUtils::getExtension(filename)=="") filename += dflt_extension; //default extension
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

		const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		try {
			getline(fin, line, eoln);      //read complete line meta information, ignore it
			if (line.length()>=3){
				if (line.substr(0,3) != "MTO") //if its not meta information rewind to the beginning
					fin.seekg (0, ios::beg);
			} else {
				throw InvalidFormatException(file_with_path.str() + ": first line in invalid format", AT);
			}

			size_t linenr = 0;

			//The following 4 lines are an optimization to jump to the correct position in the file
			streampos current_fpointer = -1;  //the filepointer for the current valid date
			map<Date,streampos>::const_iterator it = vec_streampos.at(ii).find(dateStart);
			if (it != vec_streampos.at(ii).end())
				fin.seekg(it->second); //jump to position in the file

			while (!fin.eof()){
				streampos tmp_fpointer = fin.tellg();
				getline(fin, line, eoln); //read complete line of data

				stringstream ss;
				linenr++;
				ss << linenr;

				const size_t ncols = IOUtils::readLineToVec(line, tmpvec); //split up line (whitespaces are delimiters)

				if (ncols >= nr_meteoData){
					MeteoData md;
					md.meta = vecAllStations[ii];
					parseMeteoLine(tmpvec, file_with_path.str() + ":" + ss.str(), dateStart, dateEnd, md);

					if ((md.date >= dateStart) && (md.date <= dateEnd)){//check date and add to vectors
						convertUnits(md);
						vecMeteo[ii].push_back(md);
						current_fpointer = tmp_fpointer; //save this file pointer, it's a valid one for sure
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

			//save stream position and the corresponding end date
			if (current_fpointer != ((ifstream::pos_type)-1)) vec_streampos.at(ii)[dateEnd] = current_fpointer;
		} catch (const std::exception&){
			cleanup();
			throw;
		}
		cleanup();
	}
}

void SNIO::parseMeteoLine(const std::vector<std::string>& vecLine, const std::string& filepos,
                          const Date& dateStart, const Date& dateEnd, MeteoData& md)
{
	/*
	 * This function takes a meteo line, extracts the date (ignores Julian) and then converts
	 * all meteo parameters to doubles and finally copies them into the MeteoData object md
	 */
	if (vecLine.size() < nr_meteoData)
		throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": line is too short", AT);

	if (vecLine[0] != "M")
		throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": meteo input lines must start with 'M'", AT);

	//deal with the date
	if (vecLine[1].length() != 10)
		throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": date format must be DD.MM.YYYY", AT);
	const string year  = vecLine[1].substr(6,4);
	const string month = vecLine[1].substr(3,2);
	const string day   = vecLine[1].substr(0,2);

	if (!IOUtils::convertString(md.date, year+"-"+month+"-"+day+"T"+vecLine[2], in_tz, std::dec))
		throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": invalid date format", AT);

	if ((md.date < dateStart) || (md.date > dateEnd)) //stop parsing data for dates out of the scope
		return;

	//Extract all data as double values
	vector<double> tmpdata = vector<double>(vecLine.size());
	for (size_t ii=4; ii<vecLine.size(); ii++) {
		if (!IOUtils::convertString(tmpdata[ii], vecLine[ii], std::dec))
			throw ConversionFailedException("Reading station "+md.meta.stationID+", at "+filepos+": can not convert  '"+vecLine[ii]+"' to double", AT);
	}

	//Copy data into MeteoData object
	size_t ii = 4;
	md(MeteoData::TA) = tmpdata[ii++];
	md(MeteoData::RH) = tmpdata[ii++];
	md(MeteoData::VW) = tmpdata[ii++];
	md(MeteoData::DW) = tmpdata[ii++];
	if (iswr_inp)
		md(MeteoData::ISWR) = tmpdata[ii++];
	else
		md(MeteoData::ISWR) = IOUtils::nodata;
	if (rswr_inp)
		md(MeteoData::RSWR) = tmpdata[ii++];
	else
		md(MeteoData::RSWR) = IOUtils::nodata;

	double& ea = tmpdata[ii++];
	if ((ea <= 1) && (ea != plugin_nodata)){
		if ((md(MeteoData::TA) != plugin_nodata) && (md(MeteoData::RH) != plugin_nodata)) {
			if(ea==0.)
				ea = Atmosphere::Brutsaert_ilwr(md(MeteoData::RH)/100., C_TO_K(md(MeteoData::TA)));
			else
				ea = Atmosphere::Omstedt_ilwr(md(MeteoData::RH)/100., C_TO_K(md(MeteoData::TA)), ea); //calculate ILWR from cloudiness
		} else {
			ea = plugin_nodata;
		}
	}

	md(MeteoData::ILWR) = ea;
	md(MeteoData::TSS)  = tmpdata[ii++];
	md(MeteoData::TSG)  = tmpdata[ii++];
	md(MeteoData::HNW)  = tmpdata[ii++];
	md(MeteoData::HS)   = tmpdata[ii++]; // nr_meteoData

	// Read optional values
	size_t jj;
	// TS[]: snow temperatures
	size_t number_meas_temperatures = 0;
	cfg.getValue("NUMBER_MEAS_TEMPERATURES", "Input", number_meas_temperatures, Config::nothrow);
	if (vecLine.size() < nr_meteoData + number_meas_temperatures)
		throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": not enough measured temperatures data", AT);
	stringstream ss("");
	for (jj = 1; jj <= number_meas_temperatures; jj++) {
		ss.str("");
		ss << "TS" << (jj);
		md.addParameter(ss.str());
		md(ss.str()) = tmpdata[ii++];
	}
	// CONC[]: solute concentrations
	size_t number_of_solutes = 0;
	cfg.getValue("NUMBER_OF_SOLUTES", "Input", number_of_solutes, Config::nothrow);
	if (vecLine.size() < nr_meteoData + number_meas_temperatures + number_of_solutes)
		throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": not enough solute data", AT);
	jj = 0;
	for (jj = 0 ; jj < number_of_solutes; jj++) {
		ss.str("");
		ss << "CONC" << jj;
		md.addParameter(ss.str());
		md(ss.str()) = tmpdata[ii++];
	}
	// VW_DRIFT: optional wind velocity for blowing and drifting snow
	bool vw_drift = false;
	cfg.getValue("VW_DRIFT", "Input", vw_drift, Config::nothrow);
	if (vw_drift) {
		if (vecLine.size() < ii+1)
			throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": no data for vw_drift", AT);
		md.addParameter("VW_DRIFT");
		md("VW_DRIFT") = tmpdata[ii++];
	}
	// RHO_HN: measured new snow density
	bool rho_hn = false;
	cfg.getValue("RHO_HN", "Input", rho_hn, Config::nothrow);
	if (rho_hn) {
		if (vecLine.size() < ii+1)
			throw InvalidFormatException("Reading station "+md.meta.stationID+", at "+filepos+": no data for rho_hn", AT);
		md.addParameter("RHO_HN");
		md("RHO_HN") = tmpdata[ii++];
	}
	if (vecLine.size() > ii) {
		std::stringstream ss;
		ss << "Reading station " << md.meta.stationID << ", at " << filepos << ": too many fields.\n";
		ss << "Looking for " << nr_meteoData << " standard fields + " << number_meas_temperatures << " snow temperatures + ";
		ss << number_of_solutes << " solutes";

		size_t nb_fields = nr_meteoData + number_meas_temperatures + number_of_solutes;
		if(vw_drift) {
			ss << " + 1 VW_DRIFT";
			nb_fields++;
		}
		if(rho_hn) {
			ss << " + 1 RHO_HN";
			nb_fields++;
		}
		ss << " = " << nb_fields << " fields, but found " << vecLine.size() << " fields";

		throw InvalidFormatException(ss.str(), AT);
	}
}

void SNIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	string outpath="";
	cfg.getValue("METEOPATH", "Output", outpath);

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
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

	for(size_t jj=0; jj<vecmd.size(); jj++) {
		int YYYY, MM, DD, HH, MI;
		Date tmp_date(vecmd[jj].date);
		tmp_date.setTimeZone(out_tz);
		tmp_date.getDate(YYYY, MM, DD, HH, MI);
		const double sn_julian = tmp_date.getJulianDate() - sn_julian_offset + 0.5;
		const double ta = vecmd[jj](MeteoData::TA);
		const double rh = vecmd[jj](MeteoData::RH);
		const double hnw = vecmd[jj](MeteoData::HNW);
		const double vw = vecmd[jj](MeteoData::VW);
		const double dw = vecmd[jj](MeteoData::DW);
		const double iswr = vecmd[jj](MeteoData::ISWR);
		const double rswr = vecmd[jj](MeteoData::RSWR);
		const double ilwr = vecmd[jj](MeteoData::ILWR);
		const double tss = vecmd[jj](MeteoData::TSS);
		const double tsg = vecmd[jj](MeteoData::TSG);
		const double hs = vecmd[jj](MeteoData::HS);

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
			if(tss==IOUtils::nodata) failure_count++; //if we have tss, we can compute the local ilwr
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
		for (size_t kk=1; kk<100; kk++) {
			ss.str("");
			ss << "TS" << kk;
			if (vecmd[jj].param_exists(ss.str())){
				const double ts = vecmd[jj](ss.str());
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
		for (size_t kk=0; kk<100; kk++) {
			ss.str("");
			ss << "CONC" << kk;
			if (vecmd[jj].param_exists(ss.str())) {
				const double conc = vecmd[jj](ss.str());
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
			const double vw_drift = vecmd[jj]("VW_DRIFT");
			if (vw_drift == IOUtils::nodata) {
				optional_failure_count++;
				fout << setw(4) << setprecision(0) << vw_drift << " ";
			} else {
				fout << setw(4) << setprecision(1) << vw_drift << " ";
			}
		}
		// RHO_HN: measured new snow density
		if (vecmd[jj].param_exists("RHO_HN")) {
			const double rho_hn = vecmd[jj]("RHO_HN");
			if (rho_hn == IOUtils::nodata) {
				optional_failure_count++;
				fout << setw(6) << setprecision(0) << rho_hn << " ";
			} else {
				fout << setw(6) << setprecision(1) << rho_hn << " ";
			}
		}

		fout << "\n";
	}

	fout << "END" << endl;

	if ((failure_count > 0) || (Dirichlet_failure_count > 0) || (optional_failure_count > 0)) {
		cerr << "[W] " << failure_count << " basic input data, " << Dirichlet_failure_count <<
		        " Dirichlet boundary condition data, and " << optional_failure_count <<
		        " optional data found missing when writing " << file_name << "\n";
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

void SNIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void SNIO::convertUnits(MeteoData& meteo)
{
	//converts C to Kelvin, converts ilwr to ea, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	if (ta != IOUtils::nodata) {
		if (ta < 100)
			ta = C_TO_K(ta);
	}

	double& tsg = meteo(MeteoData::TSG);
	if ( tsg != IOUtils::nodata) {
		if (tsg < 100)
			tsg = C_TO_K(tsg);
	}

	double& tss = meteo(MeteoData::TSS);
	if (tss != IOUtils::nodata) {
		if (tss < 100)
			tss = C_TO_K(tss);
	}

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata) {
		if (rh > 1.2)
			rh /= 100;
	}

	stringstream ss;
	for (size_t ii=1; ii<50; ii++){
		ss.str("");
		ss << "TS" << ii;
		if (meteo.param_exists(ss.str())){
			double& value = meteo(ss.str());
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
