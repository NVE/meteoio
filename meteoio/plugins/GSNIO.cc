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
#include "GSNIO.h"

using namespace std;

namespace mio {
/**
 * @page gsn GSN
 * @section gsn_format Format
 * This plugin reads meteorological data from GSN (Global Sensor Network, see <a href="http://sourceforge.net/apps/trac/gsn/"> GSN home page</a>) as a web service. It therefore requires GSoap.
 * @subsection gsn_fields Field mapping
 * The following GSN fields are read from GSN and mapped to MeteoData attributes:
 * RELATIVE_HUMIDITY mapped to MeteoData::RH
 * AIR_TEMPERATURE mapped to MeteoData::TA
 * WIND_DIRECTION mapped to MeteoData::DW
 * WIND_SPEED_MAX mapped to MeteoData::VW_MAX
 * WIND_SPEED_SCALAR_AV to MeteoData::VW
 * INCOMING_SHORTWAVE_RADIATION mapped to MeteoData::ISWR
 * INCOMING_LONGWAVE_RADIATION mapped to MeteoData::ILWR
 * OUTGOING_SHORTWAVE_RADIATION mapped to MeteoData::RSWR
 * OUTGOING_LONGWAVE_RADIATION mapped to parameter named "OLWR"
 * SNOW_HEIGHT mapped to MeteoData::HS
 * HACK: CURRENTLY NOTHING is mapped to MeteoData::HNW
 *
 * @section gsn_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 * - time is provided as a Unix timestamp, which is always in UTC
 *
 * @section gsn_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - ENDPOINT: The URL of the web service e.g. http://192.33.210.10:22201/services/A3DWebService/
 * - PROXY: an IP address or a resolveable hostname
 * - PROXYPORT: the port the proxy is listening on
 * - PROXYUSER: (if necessary) a proxy username
 * - PROXYPASS: (if necessary) a proxy password
 * - STATION#: station code for the given number #
 *
 * @section license Licensing
 * This software is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Part of the software embedded in this product is gSOAP software.
 * Portions created by gSOAP are Copyright (C) 2001-2009 Robert A. van Engelen, Genivia inc. All Rights Reserved.
 * THE SOFTWARE IN THIS PRODUCT WAS IN PART PROVIDED BY GENIVIA INC AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

const double GSNIO::plugin_nodata = -999.0; //plugin specific nodata value

GSNIO::GSNIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::GSNIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::GSNIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::~GSNIO() throw(){}

void GSNIO::initGSNConnection(){
	//default timezone
	string tmp_timezone = "";
	cfg.getValue("TIME_ZONE", "Input", tmp_timezone, Config::nothrow);
	if (tmp_timezone != ""){
		IOUtils::convertString(default_timezone, tmp_timezone);
	} else {
		default_timezone = 1.0;
	}

	endpoint = hostname = port = userid = passwd = "";
	proxyport = -1;

	//soap_init(&gsn);
	//soap_init2(&gsn, SOAP_IO_KEEPALIVE, SOAP_IO_KEEPALIVE);

	cfg.getValue("ENDPOINT", "INPUT", endpoint, Config::nothrow);
	if (endpoint != ""){
		gsn.soap_endpoint = endpoint.c_str();
		cout << "\tUsing GSN Endpoint: " << endpoint << endl;
	}

	/*
	 * Trying to read proxy settings:
	 * - Firstly the hostname and port (both have to be provided).
	 * - If this succeeds then the username and password will be read
	 * - parameters not set will be set to ""
	 */
	try {
		cfg.getValue("PROXY", "INPUT", hostname, Config::nothrow);
		if (hostname == "") return;
		cfg.getValue("PROXYPORT", "INPUT", port, Config::nothrow);
		if (port == "") return;

		if (!IOUtils::convertString(proxyport, port, std::dec))
			throw ConversionFailedException("", AT);
		if (proxyport < 1)
			throw IOException("",AT);

		gsn.proxy_host = hostname.c_str();
		gsn.proxy_port = proxyport;

		cfg.getValue("PROXYUSER", "INPUT", userid);
		gsn.proxy_userid = userid.c_str();
		cfg.getValue("PROXYPASS", "INPUT", passwd);
		gsn.proxy_passwd = passwd.c_str();
	} catch(...){
		//Whatever happens here can be ignored, because the proxy settings are optional
	}
}

void GSNIO::read2DGrid(Grid2DObject&, const std::string& filename)
{
	//Nothing so far
	(void)filename;
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&,
					  const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	if (vecMeta.size() == 0)
		readMetaData();

	vecStation = vecMeta;
}

void GSNIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                          std::vector< std::vector<MeteoData> >& vecMeteo,
                          const size_t& stationindex)
{
	if (vecMeta.size() == 0)
		readMetaData();

	if (vecMeta.size() == 0) //if there are no stations -> return
		return;

	unsigned int indexStart=0, indexEnd=vecMeta.size();

	//The following part decides whether all the stations are rebuffered or just one station
	if (stationindex == IOUtils::npos){
		vecMeteo.clear();
		vecMeteo.insert(vecMeteo.begin(), vecMeta.size(), vector<MeteoData>());
	} else {
		if (stationindex < vecMeteo.size()){
			indexStart = stationindex;
			indexEnd   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("You tried to access a stationindex in readMeteoData that is out of bounds", AT);
		}
	}

	for (unsigned int ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		readData(dateStart, dateEnd, vecMeteo[ii], ii);
		reverse(vecMeteo[ii].begin(), vecMeteo[ii].end());
		//cout << "vecMeteo[" <<ii << "].size() = " <<  vecMeteo[ii].size() << endl;
		for (size_t jj=0; jj<vecMeteo[ii].size(); jj++){
			//cout << vecMeteo[ii][jj] << endl;
		}
	}
}

void GSNIO::readMetaData()
{
	vecMeta.clear();

	if (vecStationName.size() == 0)
		readStationNames(); //reads station names into vector<string> vecStationName
	
	for (size_t ii=0; ii<vecStationName.size(); ii++){
		//retrieve meta info current station
		double lat=IOUtils::nodata, lon=IOUtils::nodata, alt=IOUtils::nodata;
		size_t info_complete = 0; //this variable stores whether lat,lon,alt were successfully read
		_ns1__getVirtualSensorsDetailsResponse metadata;
		_ns1__getVirtualSensorsDetails metadata_req;
		
		//Set up the SOAP request for addressing information about this virtual sensor
		ns2__GSNWebService_USCOREFieldSelector tmp;
		tmp.vsname = vecStationName[ii];
		metadata_req.detailsType.push_back(ns2__GSNWebService_USCOREDetailsType__ADDRESSING);
		metadata_req.fieldSelector.push_back(&tmp);
		
		if (gsn.getVirtualSensorsDetails(&metadata_req, &metadata) == SOAP_OK){
			size_t details = metadata.virtualSensorDetails.size();

			if (details == 0)
				throw IOException("No meta data for sensor " + vecStationName[ii], AT);

			for (size_t jj=0; jj<details; jj++){
				size_t predicates = metadata.virtualSensorDetails[jj]->addressing->predicates.size();
				for (size_t kk=0; kk<predicates; kk++){
					string field_name = metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->name;
					string field_val  = metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->__item;
					IOUtils::toUpper(field_name);
					IOUtils::toUpper(field_val);

					if (field_val != "NULL"){
						if (field_name == "LATITUDE"){
							IOUtils::convertString(lat, field_val);
							info_complete |= 1;
						} else if (field_name == "LONGITUDE"){
							IOUtils::convertString(lon, field_val);
							info_complete |= 2;
						} else if (field_name == "ALTITUDE"){
							IOUtils::convertString(alt, field_val);
							info_complete |= 4;
						}
					}

					//cout << metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->name << " -> " 
					//	<< metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->__item  << endl;
				}

				if (info_complete != 7){
					;//throw IOException("Incomplete meta data (location info) for sensor " + vecStationName[ii], AT);
				}
			}
		} else {
			soap_print_fault(&gsn, stderr);
			throw IOException("Error in communication with GSN while retrieving virtual sensor meta data", AT);
		}

		//Save the meta data in StationData objects
		Coords current_coord(coordin, coordinparam);
		current_coord.setLatLon(lat, lon, alt);
		StationData sd(current_coord, vecStationName[ii], vecStationName[ii]);

		vecMeta.push_back(sd);
	}
}
/*
void GSNIO::parseString(const std::string& in_string, MeteoData& md){
	std::stringstream ss(in_string);
	std::string tmpstring;

	while (std::getline(ss, tmpstring, ';')){
		std::stringstream data(tmpstring);
		while (std::getline(data, tmpstring, '=')){
			const string key = tmpstring;
			if (!(std::getline(data, tmpstring, '=')))
				throw InvalidFormatException("",AT);

			if (key == "LIGHT") convertStringToDouble(md.iswr, tmpstring, "ISWR");
			else if (key == "AIR_TEMP") convertStringToDouble(md.ta, tmpstring, "Air Temperature");
			else if (key == "WIND_SPEED") convertStringToDouble(md.vw, tmpstring, "Wind Velocity");
			else if (key == "WIND_DIRECTION") convertStringToDouble(md.dw, tmpstring, "Wind Velocity");
			else if (key == "INCOMING_RADIATION") {
				convertStringToDouble(md.iswr, tmpstring, "incoming_radiation");
				//convertStringToDouble(md.ilwr, tmpstring, "solar_rad");
			}
			else if (key == "RELATIVE_HUMIDITY") convertStringToDouble(md.rh, tmpstring, "relative_humidity");
			else if (key == "SURFACE_TEMP") convertStringToDouble(md.tss, tmpstring, "soil_temp");
			else if (key == "GROUND_TEMP_TNX") convertStringToDouble(md.tsg, tmpstring, "ground_temp_tnx");
			else if (key == "RAIN_METER") convertStringToDouble(md.hnw, tmpstring, "rain_meter");
			else if (key == "TIMED") {
				tmpstring = tmpstring.substr(0,tmpstring.length()-3); //cut away the seconds
				time_t measurementTime;
				if (!IOUtils::convertString(measurementTime, tmpstring, std::dec)) {
					stringstream ss;
					ss << "Conversion failed for value TIMED=" << tmpstring;
					throw ConversionFailedException(ss.str(), AT);
				}
				md.date.setDate(measurementTime);
				md.date.rnd(60., Date::CLOSEST); //HACK: round to closest minute
			}
		}
	}
}
*/

void GSNIO::map_parameters(const std::vector<ns2__GSNWebService_USCOREDataField*>& field,
                           MeteoData& md, std::vector<size_t>& index)
{
	for (size_t ii=0; ii<field.size(); ii++) {
		string field_name = *field.at(ii)->name;
		IOUtils::toUpper(field_name);

		if (field_name == "RELATIVE_HUMIDITY"){
			index.push_back(MeteoData::RH);
		} else if (field_name == "AIR_TEMPERATURE"){
			index.push_back(MeteoData::TA);
		} else if (field_name == "WIND_DIRECTION"){
			index.push_back(MeteoData::DW);
		} else if (field_name == "WIND_SPEED_MAX"){
			index.push_back(MeteoData::VW_MAX);
		} else if (field_name == "WIND_SPEED_SCALAR_AV"){
			index.push_back(MeteoData::VW);
		} else if (field_name == "INCOMING_SHORTWAVE_RADIATION"){
			index.push_back(MeteoData::ISWR);
		} else if (field_name == "INCOMING_LONGWAVE_RADIATION"){
			index.push_back(MeteoData::ILWR);
		} else if (field_name == "OUTGOING_SHORTWAVE_RADIATION"){
			index.push_back(MeteoData::RSWR);
		} else if (field_name == "OUTGOING_LONGWAVE_RADIATION"){ //is used to calculate TSS
			md.addParameter("OLWR");
			index.push_back(md.getParameterIndex("OLWR"));
		} else if (field_name == "SNOW_HEIGHT"){
			index.push_back(MeteoData::HS);
		} else if (field_name == "RAIN_METER"){
			index.push_back(MeteoData::HNW);
		} else if (field_name == "SURFACE_TEMP"){
			index.push_back(MeteoData::TSS);
		} else if (field_name == "SOLAR_RAD"){
			index.push_back(MeteoData::ISWR);
		} else {
			index.push_back(IOUtils::npos);
		}

		//cout << *field.at(ii)->name << " ";
		//cout << "(" << *field.at(ii)->type << ")  ";
	}
}

void GSNIO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	_ns1__getMultiData data_req;
	_ns1__getMultiDataResponse data;

	bool multipage = false; //there are multiple pages if there are more than 1000 tuples returned
	LONG64 start_date(dateStart.getUnixDate());
	LONG64 end_date(dateEnd.getUnixDate());
			
	start_date *= 1000; //GSN is using ms, not seconds
	end_date   *= 1000; //GSN is using ms, not seconds

	ns2__GSNWebService_USCOREFieldSelector current_station;
	current_station.vsname = vecStationName.at(stationindex);
	data_req.fieldSelector.push_back(&current_station);

	data_req.from = &start_date;
	data_req.to   = &end_date;
	//int nb = 1005; //this is a way to set the maximum amount of rows read
	//data_req.nb   = &nb;

	if (gsn.getMultiData(&data_req, &data) != SOAP_OK){
		soap_print_fault(&gsn, stderr);
		throw IOException("Error in communication with GSN while retrieving sensor data for " + current_station.vsname, AT);
	}

	if (data.queryResult.size() != 1)
		throw IOException("Incompatible data format retrieved for " + current_station.vsname, AT);

	MeteoData tmpmeteo;
	bool olwr_present = false;
	tmpmeteo.meta = vecMeta.at(stationindex);
	vector<size_t> index;

	if (data.queryResult.at(0)->format != NULL){ //data arrived
		multipage = data.queryResult[0]->hasNext;
		map_parameters(data.queryResult[0]->format->field, tmpmeteo, index);
		olwr_present = tmpmeteo.param_exists("OLWR");

		for (size_t ii=0; ii < data.queryResult.at(0)->streamElements.size(); ii++) {
			double tt;
			IOUtils::convertString(tt, *data.queryResult.at(0)->streamElements.at(ii)->timed);
			tmpmeteo.date.setUnixDate((time_t)(floor(tt/1000.0)));
			tmpmeteo.date.setTimeZone(default_timezone);

			for (size_t jj=0; jj < data.queryResult.at(0)->streamElements.at(ii)->field.size(); jj++){
				string value = data.queryResult.at(0)->streamElements.at(ii)->field.at(jj)->__item;
				IOUtils::toUpper(value);
				//cout << value << "  ";
				if (index[jj] != IOUtils::npos){
					if (value != "NULL"){
						IOUtils::convertString(tmpmeteo.param(index[jj]), value);
					} else {
						tmpmeteo.param(index[jj]) = IOUtils::nodata;
					}
				}
			}
			convertUnits(tmpmeteo);
			if ((olwr_present) && (tmpmeteo.tss == IOUtils::nodata))
				tmpmeteo.tss = olwr_to_tss(tmpmeteo.param("OLWR"));

			//cout << endl << tmpmeteo << endl;
			vecMeteo.push_back(tmpmeteo);
			tmpmeteo.tss = IOUtils::nodata; //if tss has been set, then it needs to be reset manually
		}
	}
	
	if (multipage) {
		string sid = data.queryResult.at(0)->sid;
		_ns1__getNextData requestNext;
		_ns1__getNextDataResponse responseNext;

		requestNext.sid = sid;
		bool still_pages = true;
		int result = SOAP_OK;
		
		while (still_pages && (result == SOAP_OK)){
			result = gsn.getNextData(&requestNext, &responseNext);
			still_pages = responseNext.queryResult.at(0)->hasNext;
			
			for (size_t ii=0; ii < responseNext.queryResult.at(0)->streamElements.size(); ii++) {
				double tt;
				IOUtils::convertString(tt, *responseNext.queryResult.at(0)->streamElements.at(ii)->timed);
				tmpmeteo.date.setUnixDate((time_t)(floor(tt/1000.0)));
				tmpmeteo.date.setTimeZone(default_timezone);
				
				for (size_t jj=0; jj < responseNext.queryResult.at(0)->streamElements.at(ii)->field.size(); jj++) {
					string value = responseNext.queryResult.at(0)->streamElements.at(ii)->field.at(jj)->__item;
					IOUtils::toUpper(value);
					//cout << value << "  ";
					if (index[jj] != IOUtils::npos){
						if (value != "NULL"){
							IOUtils::convertString(tmpmeteo.param(index[jj]), value);
						} else {
							tmpmeteo.param(index[jj]) = IOUtils::nodata;
						}
					}
				}
				//cout << endl << tmpmeteo << endl;
				convertUnits(tmpmeteo);
				if ((olwr_present) && (tmpmeteo.tss == IOUtils::nodata))
					tmpmeteo.tss = olwr_to_tss(tmpmeteo.param("OLWR"));

				vecMeteo.push_back(tmpmeteo);
				tmpmeteo.tss = IOUtils::nodata; //if tss has been set, then it needs to be reset manually
			}
		}		
	}			
}

double GSNIO::olwr_to_tss(const double& olwr) {
	const double ea = 1.;
	if (olwr == IOUtils::nodata) 
		return IOUtils::nodata;

	return pow( olwr / ( ea * Cst::stefan_boltzmann ), 0.25);
}

void GSNIO::readStationNames()
{
	vecStationName.clear();

	//cfg.getValue("NROFSTATIONS", "Input", str_stations, Config::nothrow);

	size_t current_stationnr = 1;
	string current_station;
	do {
		current_station = "";
		stringstream ss;
		ss << "STATION" << current_stationnr;
		cfg.getValue(ss.str(), "Input", current_station, Config::nothrow);	

		if (current_station != ""){
			vecStationName.push_back(current_station); //add station name to vector of all station names
			cout << "\tRead io.ini stationname: '" << current_station << "'" << endl;	
		}

		current_stationnr++;
	} while (current_station != "");

	if (vecStationName.size() == 0){
		//just take all sensors available
		listSensors(vecStationName);
	}
}

void GSNIO::listSensors(std::vector<std::string>& vec_names)
{
	_ns1__listVirtualSensorNamesResponse sensor_names;
	_ns1__listVirtualSensorNames sensor_req;

	if (gsn.listVirtualSensorNames(&sensor_req, &sensor_names) == SOAP_OK){
		cout << "Number of sensors accessible thorugh GSN: " << sensor_names.virtualSensorName.size() << endl;

		for (size_t ii=0; ii<sensor_names.virtualSensorName.size(); ii++){
			cout << "\tRead GSN sensor " << ii << "  name: " << sensor_names.virtualSensorName[ii] << endl;
			vec_names.push_back(sensor_names.virtualSensorName[ii]);
		}
	} else {
		soap_print_fault(&gsn, stderr);
		throw IOException("Error in communication with GSN",AT);
	}
}

void GSNIO::readAssimilationData(const Date&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::write2DGrid(const Grid2DObject&, const std::string& name)
{
	//Nothing so far
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::convertUnits(MeteoData& meteo)
{
	meteo.standardizeNodata(plugin_nodata);

	//converts C to Kelvin, converts ilwr to ea, converts RH to [0,1]
	if(meteo.ta != IOUtils::nodata) {
		meteo.ta = C_TO_K(meteo.ta);
	}

	if(meteo.tsg != IOUtils::nodata) {
		meteo.tsg = C_TO_K(meteo.tsg);
	}

	if(meteo.tss != IOUtils::nodata) {
		meteo.tss = C_TO_K(meteo.tss);
	}

	if(meteo.rh != IOUtils::nodata) {
		meteo.rh /= 100.;
	}
}

extern "C"
{
#define COMPILE_PLUGIN
#include "exports.h"

	METEOIO_EXPORT void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	METEOIO_EXPORT void* loadObject(const string& classname, const Config& cfg) {
		if(classname == "GSNIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GSNIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}

} //namespace
