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
 * The Station meteo data files is a station centered, ascii file format that has been designed with flexibility and ease of use in mind. Please refer to its <a href="../SMET_specifications.pdf">official format specification</a> for more information.
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
{//HACK: It should support coordinates in the data, ie: it should use the given date! (and TZ)
	vecStation.clear();
	vecStation.reserve(nr_stations);

	//Now loop through all requested stations, open the respective files and parse them
	for (size_t ii=0; ii<vec_smet_reader.size(); ii++){
		StationData sd;
		smet::SMETReader& myreader = vec_smet_reader[ii];

		read_meta_data(myreader, sd);
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
		size_t counter = 1;
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

		for (size_t ii=0; ii<vecFiles.size(); ii++){
			vec_smet_reader.push_back(smet::SMETReader(vecFiles[ii]));
		}
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

void SMETIO::identify_fields(const std::vector<std::string>& fields, std::vector<size_t>& indexes,
                             bool& julian_present, MeteoData& md)
{
	/*
	 * This function associates a parameter index for MeteoData objects with the
	 * lineup of field types in a SMET header. The following SMET fields are treated
	 * exceptionally:
	 * - julian, associated with IOUtils::npos
	 * - latitude, associated with IOUtils::npos-1
	 * - longitude, associated with IOUtils::npos-2
	 * - easting, associated with IOUtils::npos-3
	 * - norhting, associated with IOUtils::npos-4
	 * - altitude, associated with IOUtils::npos-5
	 * If a paramter is unknown in the fields section, then it is added as separate field to MeteoData
	 */
	for (size_t ii=0; ii<fields.size(); ii++){
		if (fields[ii] == "TA") {
			indexes.push_back(md.getParameterIndex("TA"));
		} else if (fields[ii] == "TSS") {
			indexes.push_back(md.getParameterIndex("TSS"));
		} else if (fields[ii] == "TSG") {
			indexes.push_back(md.getParameterIndex("TSG"));
		} else if (fields[ii] == "RH") {
			indexes.push_back(md.getParameterIndex("RH"));
		} else if (fields[ii] == "VW") {
			indexes.push_back(md.getParameterIndex("VW"));
		} else if (fields[ii] == "VW_MAX") {
			indexes.push_back(md.getParameterIndex("VW_MAX"));
		} else if (fields[ii] == "DW") {
			indexes.push_back(md.getParameterIndex("DW"));
		} else if (fields[ii] == "ISWR") {
			indexes.push_back(md.getParameterIndex("ISWR"));
		} else if (fields[ii] == "OSWR") {
			indexes.push_back(md.getParameterIndex("RSWR"));
		} else if (fields[ii] == "ILWR") {
			indexes.push_back(md.getParameterIndex("ILWR"));
		} else if (fields[ii] == "OLWR") {
			md.addParameter("OLWR");
			indexes.push_back(md.getParameterIndex("OLWR"));
		} else if (fields[ii] == "PINT") {
			md.addParameter("PINT");
			indexes.push_back(md.getParameterIndex("PINT"));
		} else if (fields[ii] == "PSUM") {
			indexes.push_back(md.getParameterIndex("HNW"));
		} else if (fields[ii] == "HS") {
			indexes.push_back(md.getParameterIndex("HS"));
		} else if (fields[ii] == "julian") {
			julian_present = true;
			indexes.push_back(IOUtils::npos);
		} else if (fields[ii] == "latitude") {
			indexes.push_back(IOUtils::npos-1);
		} else if (fields[ii] == "longitude") {
			indexes.push_back(IOUtils::npos-2);
		} else if (fields[ii] == "easting") {
			indexes.push_back(IOUtils::npos-3);
		} else if (fields[ii] == "northing") {
			indexes.push_back(IOUtils::npos-4);
		} else if (fields[ii] == "altitude") {
			indexes.push_back(IOUtils::npos-5);
		} else {
			md.addParameter(fields[ii]);
			indexes.push_back(md.getParameterIndex(fields[ii]));
		}
	}
}

void SMETIO::read_meta_data(const smet::SMETReader& myreader, StationData& meta)
{
	/*
	 * This function reads in the header data provided by a SMETReader object.
	 * SMETReader objects read all the header info upon construction and can subsequently
	 * be queried for that info
	 */
	double nodata_value = myreader.get_header_doublevalue("nodata");

	meta.position.setProj(coordin, coordinparam); //set the default projection from config file
	if (myreader.location_in_header(smet::WGS84)){
		double lat = myreader.get_header_doublevalue("latitude");
		double lon = myreader.get_header_doublevalue("longitude");
		double alt = myreader.get_header_doublevalue("altitude");
		meta.position.setLatLon(lat, lon, alt);
	}

	if (myreader.location_in_header(smet::EPSG)){
		double east  = myreader.get_header_doublevalue("easting");
		double north = myreader.get_header_doublevalue("northing");
		double alt   = myreader.get_header_doublevalue("altitude");
		short int epsg  = (short int)(floor(myreader.get_header_doublevalue("epsg") + 0.1));
		meta.position.setEPSG(epsg); //this needs to be set before calling setXY(...)
		meta.position.setXY(east, north, alt);
	}

	meta.stationID = myreader.get_header_value("station_id");
	meta.stationName = myreader.get_header_value("station_name");

	bool data_epsg = myreader.location_in_data(smet::EPSG);
	if (data_epsg){
		double d_epsg = myreader.get_header_doublevalue("epsg");
		short int epsg = IOUtils::snodata;
		if (d_epsg != nodata_value)
			epsg  = (short int)(floor(d_epsg + 0.1));

		meta.position.setEPSG(epsg);
	}
}

void SMETIO::copy_data(const smet::SMETReader& myreader,
                       const std::vector<std::string>& timestamps,
                       const std::vector<double>& mydata, std::vector<MeteoData>& vecMeteo)
{
	/*
	 * This function parses the data read from a SMETReader object, a vector<double>,
	 * and copies the values into their respective places in the MeteoData structure
	 * Meta data, whether in header or in data is also handled
	 */
	vector<size_t> indexes;
	MeteoData md;
	bool julian_present = false;
	vector<string> fields;
	string myfields = myreader.get_header_value("fields");
	IOUtils::readLineToVec(myfields, fields);
	identify_fields(fields, indexes, julian_present, md);

	if ((timestamps.size() == 0) && (!julian_present)) return; //nothing to do

	bool olwr_present = md.param_exists("OLWR");

	bool data_wgs84 = myreader.location_in_data(smet::WGS84);
	bool data_epsg = myreader.location_in_data(smet::EPSG);

	read_meta_data(myreader, md.meta);

	double nodata_value = myreader.get_header_doublevalue("nodata");
	double current_timezone = myreader.get_header_doublevalue("timezone");
	if (current_timezone == nodata_value)
		current_timezone = in_dflt_TZ;

	Date tmp_date;
	size_t nr_of_lines = mydata.size() / indexes.size();
	bool timestamp_present = myreader.contains_timestamp();

	double lat=IOUtils::nodata, lon=IOUtils::nodata, east=IOUtils::nodata, north=IOUtils::nodata, alt=IOUtils::nodata;
	size_t current_index = 0; //index to vec_data
	for (size_t ii = 0; ii<nr_of_lines; ii++){
		if (timestamp_present)
			IOUtils::convertString(md.date, timestamps[ii], current_timezone);

		//cout << md.date.toString(Date::ISO) << ": ";

		vecMeteo.push_back(md);
		MeteoData& tmp_md = vecMeteo.back();

		//Copy data points
		for (size_t jj=0; jj<indexes.size(); jj++){
			const double& current_data = mydata[current_index];
			//cout << current_data << " ";
			if (indexes[jj] >= IOUtils::npos-5){ //the special fields have high indexes
				if (indexes[jj] == IOUtils::npos){
					if (!timestamp_present){
						if (current_data != nodata_value)
							tmp_md.date.setDate(current_data, current_timezone);
					}
				} else if (indexes[jj] == IOUtils::npos-1){
					lat = current_data;
				} else if (indexes[jj] == IOUtils::npos-2){
					lon = current_data;
				} else if (indexes[jj] == IOUtils::npos-3){
					east = current_data;
				} else if (indexes[jj] == IOUtils::npos-4){
					north = current_data;
				} else if (indexes[jj] == IOUtils::npos-5){
					alt = current_data;
				}
			} else {
				if (current_data == nodata_value)
					tmp_md.param(indexes[jj]) = IOUtils::nodata;
				else
					tmp_md.param(indexes[jj]) = current_data;
			}

			if (data_epsg)
				tmp_md.meta.position.setXY(east, north, alt);

			if (data_wgs84)
				tmp_md.meta.position.setXY(lat, lon, alt);

			current_index++;
		}

		if ((olwr_present) && (md.tss == IOUtils::nodata)) //HACK
			md.tss = olwr_to_tss(md.param("OLWR"));

		//cout << tmp_md << endl;
		//cout << endl;
	}
}

double SMETIO::olwr_to_tss(const double& olwr) {
	const double ea = 1.;
	if(olwr==IOUtils::nodata) return IOUtils::nodata;
	return pow( olwr / ( ea * Cst::stefan_boltzmann ), 0.25);
}

void SMETIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                           std::vector< std::vector<MeteoData> >& vecMeteo,
                           const size_t& stationindex)
{
	//Make sure that vecMeteo have the correct dimension and stationindex is valid
	size_t startindex=0, endindex=vecFiles.size();
	if (stationindex != (size_t)IOUtils::npos){ //HACK do we really still need stationindex??
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
	for (size_t ii=startindex; ii<endindex; ii++){
		string filename = vecFiles.at(ii); //filename of current station

		if (!IOUtils::fileExists(filename))
			throw FileNotFoundException(filename, AT);

		smet::SMETReader& myreader = vec_smet_reader.at(ii);
		myreader.convert_to_MKSA(true); // we want converted values for MeteoIO

		vector<double> mydata; //sequentially store all data in the smet file
		vector<string> mytimestamps;

		if (myreader.contains_timestamp()){
			myreader.read(dateStart.toString(Date::ISO), dateEnd.toString(Date::ISO), mytimestamps, mydata);
		} else {
			myreader.read(dateStart.getJulianDate(), dateEnd.getJulianDate(), mydata);
		}

		copy_data(myreader, mytimestamps, mydata, vecMeteo[ii]);
	}
}

void SMETIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	//Loop through all stations
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		//1. check consitency of station data position -> write location in header or data section
		StationData sd;
		sd.position.setProj(coordout, coordoutparam);
		const bool isConsistent = checkConsistency(vecMeteo.at(ii), sd);

		if (sd.stationID == ""){
			stringstream ss;
			ss << "Station" << ii+1;
			sd.stationID = ss.str();
		}

		const string filename = outpath + "/" + sd.stationID + ".smet";
		if (!IOUtils::validFileName(filename)) //Check whether filename is valid
			throw InvalidFileNameException(filename, AT);

		//2. check which meteo parameter fields are actually in use
		vector<bool> vecParamInUse = vector<bool>(MeteoData::nrOfParameters, false);
		double timezone = IOUtils::nodata;
		checkForUsedParameters(vecMeteo[ii], timezone, vecParamInUse);

		stringstream ss;
		try {
			smet::SMETType type = smet::ASCII;
			if (!outputIsAscii) type = smet::BINARY;

			smet::SMETWriter mywriter(filename, type, outputIsGzipped);
			mywriter.set_header_value("station_id", sd.stationID);
			if (sd.stationName != "")
				mywriter.set_header_value("station_name", sd.stationName);
			mywriter.set_header_value("nodata", IOUtils::nodata);

			vector<int> myprecision, mywidth; //set meaningful precision/width for each column

			if (outputIsAscii) {
				ss << "timestamp";
			} else {
				ss << "julian";
				myprecision.push_back(8);
				mywidth.push_back(16);
			}

			if (isConsistent) {
				mywriter.set_header_value("latitude", sd.position.getLat());
				mywriter.set_header_value("longitude", sd.position.getLon());
				mywriter.set_header_value("easting", sd.position.getEasting());
				mywriter.set_header_value("northing", sd.position.getNorthing());
				mywriter.set_header_value("altitude", sd.position.getAltitude());
				mywriter.set_header_value("epsg", (double)sd.position.getEPSG());

				if ((timezone != IOUtils::nodata) && (timezone != 0.0))
					mywriter.set_header_value("tz", timezone);
			} else {
				ss << " latitude longitude altitude";
				myprecision.push_back(8); //for latitude
				mywidth.push_back(11);    //for latitude
				myprecision.push_back(8); //for longitude
				mywidth.push_back(11);    //for longitude
				myprecision.push_back(1); //for altitude
				mywidth.push_back(7);     //for altitude
			}

			//Add all other used parameters
			int tmpwidth, tmpprecision;
			for (size_t ll=0; ll<MeteoData::nrOfParameters; ll++){
				if (vecParamInUse[ll]) {
					std::string column=MeteoData::getParameterName(ll);
					if(column=="RSWR") column="OSWR";
					if(column=="HNW") column="PSUM";
					ss << " " << column;

					getFormatting(ll, tmpprecision, tmpwidth);
					myprecision.push_back(tmpprecision);
					mywidth.push_back(tmpwidth);
				}
			}
			mywriter.set_header_value("fields", ss.str());
			mywriter.set_width(mywidth);
			mywriter.set_precision(myprecision);

			vector<string> vec_timestamp;
			vector<double> vec_data;
			for (size_t jj=0; jj<vecMeteo[ii].size(); jj++){
				if (outputIsAscii){
					if(out_dflt_TZ!=IOUtils::nodata) {
						Date tmp_date(vecMeteo[ii][jj].date);
						tmp_date.setTimeZone(out_dflt_TZ);
						vec_timestamp.push_back(tmp_date.toString(Date::ISO));
					} else {
						vec_timestamp.push_back(vecMeteo[ii][jj].date.toString(Date::ISO));
					}
				} else {
					double julian;
					if(out_dflt_TZ!=IOUtils::nodata) {
						Date tmp_date(vecMeteo[ii][jj].date);
						tmp_date.setTimeZone(out_dflt_TZ);
						julian = tmp_date.getJulianDate();
					} else {
						julian = vecMeteo[ii][jj].date.getJulianDate();
					}
					vec_data.push_back(julian);
				}

				if (!isConsistent){ //Meta data changes
					vec_data.push_back(vecMeteo[ii][jj].meta.position.getLat());
					vec_data.push_back(vecMeteo[ii][jj].meta.position.getLon());
					vec_data.push_back(vecMeteo[ii][jj].meta.position.getAltitude());
				}

				for (size_t kk=0; kk<MeteoData::nrOfParameters; kk++){
					if (vecParamInUse[kk])
						vec_data.push_back(vecMeteo[ii][jj].param(kk)); //add data value
				}
			}

			if (outputIsAscii) mywriter.write(vec_timestamp, vec_data);
			else mywriter.write(vec_data);

		} catch(const exception& e) {
			throw;
		}
	}
}

void SMETIO::getFormatting(const size_t& param, int& prec, int& width)
{
	if ((param == MeteoData::TA) || (param == MeteoData::TSS) || (param == MeteoData::TSG)){
		prec = 2;
		width = 8;
	} else if ((param == MeteoData::VW) || (param == MeteoData::VW_MAX)){
		prec = 1;
		width = 6;
	} else if (param == MeteoData::DW){
		prec = 0;
		width = 5;
	} else if ((param == MeteoData::ISWR) || (param == MeteoData::RSWR) || (param == MeteoData::ILWR)){
		prec = 0;
		width = 6;
	} else if (param == MeteoData::HNW){
		prec = 3;
		width = 6;
	} else if (param == MeteoData::HS){
		prec = 3;
		width = 8;
	} else if (param == MeteoData::RH){
		prec = 3;
		width = 7;
	} else {
		prec = 3;
		width = 8;
	}
}

void SMETIO::checkForUsedParameters(const std::vector<MeteoData>& vecMeteo, double& timezone,
                                    std::vector<bool>& vecParamInUse)
{
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		for (size_t jj=0; jj<MeteoData::nrOfParameters; jj++){
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

	for (size_t ii=1; ii<vecMeteo.size(); ii++){
		const Coords& p1 = vecMeteo[ii-1].meta.position;
		const Coords& p2 = vecMeteo[ii].meta.position;
		if (p1 != p2) {
			//we don't mind if p1==nodata or p2==nodata
			if(p1.isNodata()==false && p2.isNodata()==false) return false;
		}
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
