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
 * <center><table border="0">
 * <tr><td>
 * <table border="1">
 * <tr><th>GSN attribute</th><th>MeteoData field</th></tr>
 * <tr><td>RELATIVE_HUMIDITY</td><td>MeteoData::RH</td></tr>
 * <tr><td>AIR_TEMPERATURE</td><td>MeteoData::TA</td></tr>
 * <tr><td>WIND_DIRECTION</td><td>MeteoData::DW</td></tr>
 * <tr><td>WIND_SPEED_MAX</td><td>MeteoData::VW_MAX</td></tr>
 * <tr><td>WIND_SPEED_SCALAR_AV</td><td>MeteoData::VW</td></tr>
 * <tr><td>INCOMING_SHORTWAVE_RADIATION</td><td>MeteoData::ISWR</td></tr>
 * <tr><td>INCOMING_LONGWAVE_RADIATION</td><td>MeteoData::ILWR</td></tr>
 * <tr><td>OUTGOING_SHORTWAVE_RADIATION</td><td>MeteoData::RSWR</td></tr>
 * <tr><td>OUTGOING_LONGWAVE_RADIATION</td><td>equivalent MeteoData::TSS</td></tr>
 * <tr><td>SNOW_HEIGHT</td><td>MeteoData::HS</td></tr>
 * <tr><td>RAIN_METER</td><td>MeteoData::HNW</td></tr>
 * <tr><td>SURFACE_TEMP</td><td>MeteoData::TSS</td></tr>
 * <tr><td>SOLAR_RAD</td><td>MeteoData::ISWR</td></tr>
 * </table></td></tr>
 * </table></center>
 * Please keep in mind that the names in GSN have currently not been standardized. This means that any sensor that does not use the above names will not be properly supported (some fields might be missing)!
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
 * - ENDPOINT: The URL of the web service e.g. http://planetdata.epfl.ch:22001/services/GSNWebService/
 * - PROXY: an IP address or a resolveable hostname
 * - PROXYPORT: the port the proxy is listening on
 * - PROXYUSER: (if necessary) a proxy username
 * - PROXYPASS: (if necessary) a proxy password
 * - STATION#: station code for the given number #
 *
 * If no STATION keys are given, the full list of ALL stations in GSN will be printed out and used.
 *
 * @section license Licensing
 * This software is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Part of the software embedded in this product is gSOAP software.
 * Portions created by gSOAP are Copyright (C) 2001-2009 Robert A. van Engelen, Genivia inc. All Rights Reserved.
 * THE SOFTWARE IN THIS PRODUCT WAS IN PART PROVIDED BY GENIVIA INC AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
		reverse(vecMeteo[ii].begin(), vecMeteo[ii].end()); //this is necessary because GSN data comes sorted descending by date
		/*//The following block can be commented in for testing purposes
		cout << "vecMeteo[" <<ii << "].size() = " <<  vecMeteo[ii].size() << endl;
		for (size_t jj=0; jj<vecMeteo[ii].size(); jj++){
			cout << vecMeteo[ii][jj] << endl;
		}*/
	}
}

void GSNIO::readMetaData()
{
	vecMeta.clear();

	if (vecStationName.size() == 0)
		readStationNames(); //reads station names into vector<string> vecStationName

	for (size_t ii=0; ii<vecStationName.size(); ii++) {
		//retrieve meta info current station
		std::string name;
		double lat=IOUtils::nodata, lon=IOUtils::nodata, alt=IOUtils::nodata;
		double slope_angle = IOUtils::nodata;
		double slope_azi = IOUtils::nodata;
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
				throw IOException("No meta data for sensor " + vecStationName[ii] + ". Is it a valid sensor?", AT);

			for (size_t jj=0; jj<details; jj++){
				size_t predicates = metadata.virtualSensorDetails[jj]->addressing->predicates.size();
				for (size_t kk=0; kk<predicates; kk++){
					string field_name = metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->name;
					string field_val  = metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->__item;
					IOUtils::toUpper(field_name);
					IOUtils::toUpper(field_val);

					if (field_val != "NULL") {
						if (field_name == "NAME") {
							IOUtils::convertString(name, field_val);
							info_complete |= 1;
						} else if (field_name == "LATITUDE") {
							IOUtils::convertString(lat, field_val);
							info_complete |= 2;
						} else if (field_name == "LONGITUDE") {
							IOUtils::convertString(lon, field_val);
							info_complete |= 4;
						} else if (field_name == "ALTITUDE") {
							IOUtils::convertString(alt, field_val);
							info_complete |= 8;
						} else if (field_name == "SLOPE") {
							IOUtils::convertString(slope_angle, field_val);
						} else if (field_name == "EXPOSITION") {
							std::string tmp;
							IOUtils::convertString(tmp, field_val);
							if (IOUtils::isNumeric(tmp)) IOUtils::convertString(slope_azi, field_val);
							else slope_azi=IOUtils::bearing(tmp);
							info_complete |= 16;
						}
					}

					//cout << metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->name << " -> "
					//	<< metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->__item  << endl;
				}

				if (info_complete != 31){
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
		StationData sd(current_coord, vecStationName[ii], name);

		if (slope_angle != IOUtils::nodata){
			if ((slope_angle==0.) && (slope_azi==IOUtils::nodata)) {
				sd.setSlope(slope_angle, 0.); //expostion: north assumed
			} else {
				sd.setSlope(slope_angle, slope_azi);
			}
		}

		vecMeta.push_back(sd);
	}
}

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
			parse_streamElement(index, olwr_present, vecMeteo, tmpmeteo, data.queryResult.at(0)->streamElements.at(ii));
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
				parse_streamElement(index, olwr_present, vecMeteo, tmpmeteo, responseNext.queryResult.at(0)->streamElements.at(ii));
			}
		}
	}
}

void GSNIO::parse_streamElement(const std::vector<size_t>& index, const bool& olwr_present,
                                std::vector<MeteoData>& vecMeteo, MeteoData& tmpmeteo, ns2__GSNWebService_USCOREStreamElement* streamElement)
{
	/**
	 * This procedure takes a streamElement pointer from either the _ns1__getNextDataResponse
	 * or _ns1__getMultiDataResponse and parses the field elements into a MeteoData object
	 * (tmpmeteo). Finally it adjusts the units and calculates TSS from OLWR if necessary and
	 * possible.
	 */
	double tt;
	IOUtils::convertString(tt, *streamElement->timed);
	tmpmeteo.date.setUnixDate((time_t)(floor(tt/1000.0)));
	tmpmeteo.date.setTimeZone(default_timezone);

	for (size_t jj=0; jj < streamElement->field.size(); jj++){
		string value = streamElement->field.at(jj)->__item;
		IOUtils::toUpper(value);
		//cout << value << "  ";
		if (index[jj] != IOUtils::npos){
			if (value != "NULL"){
				IOUtils::convertString(tmpmeteo(index[jj]), value);
			} else {
				tmpmeteo(index[jj]) = IOUtils::nodata;
			}
		}
	}

	convertUnits(tmpmeteo);
	if ((olwr_present) && (tmpmeteo(MeteoData::TSS) == IOUtils::nodata))
		tmpmeteo(MeteoData::TSS) = olwr_to_tss(tmpmeteo("OLWR"));

	//cout << endl << tmpmeteo << endl;
	vecMeteo.push_back(tmpmeteo);
	tmpmeteo(MeteoData::TSS) = IOUtils::nodata; //if tss has been set, then it needs to be reset manually
}

double GSNIO::olwr_to_tss(const double& olwr) {
	const double ea = 1.;
	if (olwr == IOUtils::nodata)
		return IOUtils::nodata;

	return pow( olwr / ( ea * Cst::stefan_boltzmann ), 0.25);
}

void GSNIO::readStationNames()
{
	/**
	 * Parse through the io.ini file and copy the desired station names STATION#
	 * into vecStationName, if no stations are configured explicitly simply take
	 * all stations that are available in the current GSN instance
	 */
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
			cout << "\tRead stationname '" << current_station << "'" << endl;
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
	/**
	 * Retrieve all station names, that are available in the current GSN instance
	 */
	_ns1__listVirtualSensorNamesResponse sensor_names;
	_ns1__listVirtualSensorNames sensor_req;

	if (gsn.listVirtualSensorNames(&sensor_req, &sensor_names) == SOAP_OK){
		cout << "[I] No STATIONS specified... Using all " << sensor_names.virtualSensorName.size() << " sensors available through GSN" << endl;
		for (size_t ii=0; ii<sensor_names.virtualSensorName.size(); ii++){
			cout << "\tSTATION" << ii+1 << " = " << sensor_names.virtualSensorName[ii] << endl;
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
