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
 * - URL: The URL of the RESTful web service e.g. http://planetdata.epfl.ch:22001/rest
 * - USER: The username to access the service
 * - PASS: The password to authenticate the USER
 * - STATION#: station code for the given number #, e. g. la_fouly_1034 (case sensitive!)
 *
 * If no STATION keys are given, the full list of ALL stations available to the user in GSN will be used!
 * This may result in a long download.
 *
 */

const int GSNIO::http_timeout = 20; // seconds until connect time out for libcurl
const std::string GSNIO::sensors_endpoint = "sensors";

GSNIO::GSNIO(const std::string& configfile)
      : cfg(configfile), vecStationName(), vecMeta(), vecAllMeta(), coordin(), coordinparam(), coordout(), coordoutparam(),
        endpoint(), userid(), passwd(), default_timezone(1.)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::GSNIO(const Config& cfgreader)
      : cfg(cfgreader), vecStationName(), vecMeta(), vecAllMeta(), coordin(), coordinparam(), coordout(), coordoutparam(),
        endpoint(), userid(), passwd(), default_timezone(1.)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::~GSNIO() throw(){}

void GSNIO::initGSNConnection() {
	curl_global_init(CURL_GLOBAL_ALL);

	default_timezone = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", default_timezone, IOUtils::nothrow);

	cfg.getValue("URL", "Input", endpoint, IOUtils::nothrow);
	if (!endpoint.empty()){
		if (*endpoint.rbegin() != '/') endpoint += "/";
		cerr << "\tUsing GSN Endpoint: " << endpoint << endl;
	}

	cfg.getValue("USER", "Input", userid);
	cfg.getValue("PASS", "Input", passwd);
}

void GSNIO::read2DGrid(Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
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

	if (vecMeta.empty())
		readMetaData();

	vecStation = vecMeta;
}
void GSNIO::readMetaData()
{	
	vecMeta.clear();

	if (vecStationName.empty())
		readStationNames(); //reads station names into vector<string> vecStationName

	
	//Get Meta Data for all stations first
	getAllStations();

	if (!vecStationName.empty()) { //if the user has specified a subset of stations
		for (size_t ii=0; ii<vecStationName.size(); ii++) {
			for (size_t jj=0; jj<vecAllMeta.size(); jj++) {
				if (vecAllMeta[jj].stationID == vecStationName[ii]) {
					vecMeta.push_back(vecAllMeta[jj]);
				}
			}
		}
	} else { //otherwise use all available stations
		vecMeta = vecAllMeta;
	}
}

void GSNIO::readStationNames()
{
	/**
	 * Parse through the io.ini file and copy the desired station names STATION#
	 * into vecStationName, if no stations are configured all stations with 
	 * meteo data will be taken into account
	 */
	vecStationName.clear();

	size_t current_stationnr = 1;
	string current_station;

	do {
		current_station = string("");
		ostringstream ss;
		ss << "STATION" << current_stationnr;
		cfg.getValue(ss.str(), "Input", current_station, IOUtils::nothrow);

		if (!current_station.empty()){
			vecStationName.push_back(current_station); //add station name to vector of all station names
			cerr << "\tRead stationname '" << current_station << "'\n";
		}

		current_stationnr++;
	} while (!current_station.empty());

	/*
	if (vecStationName.empty()){
		//just take all sensors available
		listSensors(vecStationName);
	}
	*/
}

void GSNIO::save_station(const std::string& id, const std::string& name, const double& lat, const double& lon, 
                         const double& alt, const double& slope_angle, const double& slope_azi)
{
	Coords current_coord(coordin, coordinparam);
	current_coord.setLatLon(lat, lon, alt);
	StationData sd(current_coord, id, name);
					
	if (slope_angle != IOUtils::nodata){
		if ((slope_angle == 0.) && (slope_azi == IOUtils::nodata)) {
			sd.setSlope(slope_angle, 0.); //expostion: north assumed
		} else {
			sd.setSlope(slope_angle, slope_azi);
		}
	}

	vecAllMeta.push_back(sd);
}

void GSNIO::getAllStations()
{
	/**
	 * Retrieve all station names, that are available in the current GSN instance
	 * and which are accessible for the current user (see Input::USER)
	 */
	const string vsname_str = "# vsname:";
	const string altitude_str = "# altitude:";
	const string longitude_str = "# longitude:";
	const string latitude_str = "# latitude:";
	const string slope_str = "# slope:";
	const string exposition_str = "# exposition:";
	const string name_str = "# name:";

	stringstream ss;
	string line = "";

	vecAllMeta.clear();	

	if (curl_read(sensors_endpoint + "?username=" + userid + "&password=" + passwd, ss) == CURLE_OK) {
		string name="", id="", azi="";
		double lat=0., lon=0., alt=0., slope_angle=IOUtils::nodata, slope_azi=IOUtils::nodata;
		unsigned int valid = 0;

		while (getline(ss, line)) {
			if (!line.compare(0, vsname_str.size(), vsname_str)) {

				if (valid == 15) { // Last station was valid: store StationData
					save_station(id, name, lat, lon, alt, slope_angle, slope_azi);
				}

				id = line.substr(vsname_str.size());
				IOUtils::trim(id);
				slope_angle = slope_azi = IOUtils::nodata;
				name  = azi = "";
				valid = 1;
			} else if (!line.compare(0, altitude_str.size(), altitude_str)) {
				IOUtils::convertString(alt, line.substr(altitude_str.size()));
				valid |= 2;
			} else if (!line.compare(0, latitude_str.size(), latitude_str)) {
				IOUtils::convertString(lat, line.substr(latitude_str.size()));
				valid |= 4;
			} else if (!line.compare(0, longitude_str.size(), longitude_str)) {
				IOUtils::convertString(lon, line.substr(longitude_str.size()));
				valid |= 8;
			} else if (!line.compare(0, name_str.size(), name_str)) { // optional
				name = line.substr(name_str.size());
				IOUtils::trim(name);
			} else if (!line.compare(0, slope_str.size(), slope_str)) { //optional
				IOUtils::convertString(slope_angle, line.substr(slope_str.size()));
			} else if (!line.compare(0, exposition_str.size(), exposition_str)) { //optional
				azi = line.substr(exposition_str.size());
				if (IOUtils::isNumeric(azi)) {
					IOUtils::convertString(slope_azi, azi);
				} else {
					slope_azi = IOUtils::bearing(azi);
				}
			}
		}

		if (valid == 15) { // Last station was valid: store StationData
			save_station(id, name, lat, lon, alt, slope_angle, slope_azi);
		}
	} else {
		throw IOException("Could not retrieve list of sensors", AT);
	}
}

void GSNIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                          std::vector< std::vector<MeteoData> >& vecMeteo,
                          const size_t& stationindex)
{
	if (vecMeta.empty())
		readMetaData();

	if (vecMeta.empty()) //if there are no stations -> return
		return;

	size_t indexStart=0, indexEnd=vecMeta.size();

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

	for (size_t ii=indexStart; ii<indexEnd; ii++){ //loop through stations
		readData(dateStart, dateEnd, vecMeteo[ii], ii);
	}
}

void GSNIO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	stringstream ss;
	string line = "";
	bool fields_detected = false;
	cout << "End date: " << dateEnd.toString(Date::ISO) << endl;
	string request = sensors_endpoint + "/" + vecMeta[stationindex].stationID + "?from=" + dateStart.toString(Date::ISO) + ":00"
	                 + "&to=" + dateEnd.toString(Date::ISO) + ":00" + "&username=" + userid + "&password=" + passwd;

	cout << "Requesting: " << request << endl;
	
	if (curl_read(request, ss) == CURLE_OK) {
		while (getline(ss, line)) {
			MeteoData tmpmeteo;
			vector<size_t> index;
			bool olwr_present = false;
			tmpmeteo.meta = vecMeta.at(stationindex);

			if (!fields_detected && (line.size() && (line[0] == '#'))) { //Detect fields, it's the first line without a colon
				size_t found = line.find(":");
				if (found == string::npos) {
					fields_detected = true;
					map_parameters(line, tmpmeteo, index);
					olwr_present = tmpmeteo.param_exists("OLWR");
				}
			}

			if (line.size() && (line[0] != '#')) { // data
				parse_streamElement(line, index, olwr_present, vecMeteo, tmpmeteo);
			}
		}
	} else {
		throw IOException("Could not retrieve data for station " + vecMeta[stationindex].stationID, AT);
	}
}

void GSNIO::map_parameters(const std::string& fields, MeteoData& md, std::vector<size_t>& index)
{
	std::vector<std::string> field;

	cout << fields << endl;
	IOUtils::readLineToVec(fields.substr(2), field, ',');

	for (size_t ii=0; ii<field.size(); ii++) {
		const string field_name = IOUtils::strToUpper(field[ii]);

		if (field_name == "RELATIVE_HUMIDITY" || field_name == "RH" || field_name == "AIR_HUMID"){
			index.push_back(MeteoData::RH);
		} else if (field_name == "AIR_TEMPERATURE" || field_name == "TA" || field_name == "AIR_TEMP"){
			index.push_back(MeteoData::TA);
		} else if (field_name == "WIND_DIRECTION" || field_name == "DW"){
			index.push_back(MeteoData::DW);
		} else if (field_name == "WIND_SPEED_MAX" || field_name == "VW_MAX"){
			index.push_back(MeteoData::VW_MAX);
		} else if (field_name == "WIND_SPEED_SCALAR_AV" || field_name == "VW" || field_name == "WIND_SPEED"){
			index.push_back(MeteoData::VW);
		} else if (field_name == "INCOMING_SHORTWAVE_RADIATION" || field_name == "ISWR" || field_name == "SOLAR_RAD"){
			index.push_back(MeteoData::ISWR);
		} else if (field_name == "INCOMING_LONGWAVE_RADIATION" || field_name == "ILWR"){
			index.push_back(MeteoData::ILWR);
		} else if (field_name == "OUTGOING_SHORTWAVE_RADIATION" || field_name == "RSWR"){
			index.push_back(MeteoData::RSWR);
		} else if (field_name == "OUTGOING_LONGWAVE_RADIATION" || field_name == "RLWR"){ //is used to calculate TSS
			md.addParameter("OLWR");
			index.push_back(md.getParameterIndex("OLWR"));
		} else if (field_name == "SNOW_HEIGHT" || field_name == "HS1"){
			index.push_back(MeteoData::HS);
		} else if (field_name == "RAIN_METER" || field_name == "PINT"){
			index.push_back(MeteoData::HNW);
		} else if (field_name == "SURFACE_TEMP" || field_name == "TSS"){
			index.push_back(MeteoData::TSS);
		} else {
			index.push_back(IOUtils::npos);
		}
	}
}

void GSNIO::parse_streamElement(const std::string& line, const std::vector<size_t>& index, const bool& olwr_present, std::vector<MeteoData>& vecMeteo, MeteoData& tmpmeteo)
{
	vector<string> data;
	IOUtils::readLineToVec(line, data, ',');	
	if (data.size() < 2) return; // Malformed for sure

	cout << "Attempting to parse: " << line << endl;

	//The first or the second element is a UNIX timestamp, let's see which one
	double timestamp;
	if (IOUtils::isNumeric(data[0])) {
		IOUtils::convertString(timestamp, data[0]);
	} else if (IOUtils::isNumeric(data[1])) {
		IOUtils::convertString(timestamp, data[1]);
	} else {
		throw InvalidFormatException("No date string present when reading meteo data for station + " + tmpmeteo.meta.stationID, AT);
	}

	tmpmeteo.date.setUnixDate((time_t)(floor(timestamp/1000.0)));
	tmpmeteo.date.setTimeZone(default_timezone);

	for (size_t jj=2; jj<data.size(); jj++){
		const string value = IOUtils::strToUpper(data[jj]);
		if (index[jj-2] != IOUtils::npos){
			if (value != "NULL"){
				IOUtils::convertString(tmpmeteo(index[jj-2]), value);
			} else {
				tmpmeteo(index[jj-2]) = IOUtils::nodata;
			}
		}
	}

	convertUnits(tmpmeteo);
	if ((olwr_present) && (tmpmeteo(MeteoData::TSS) == IOUtils::nodata))
		tmpmeteo(MeteoData::TSS) = olwr_to_tss(tmpmeteo("OLWR"));

	vecMeteo.push_back(tmpmeteo);
	tmpmeteo(MeteoData::TSS) = IOUtils::nodata; //if tss has been set, then it needs to be reset manually

	/**
	 * This procedure takes a streamElement pointer from either the _ns1__getNextDataResponse
	 * or _ns1__getMultiDataResponse and parses the field elements into a MeteoData object
	 * (tmpmeteo). Finally it adjusts the units and calculates TSS from OLWR if necessary and
	 * possible.
	 */
	/*
	double tt;
	IOUtils::convertString(tt, *streamElement->timed);
	tmpmeteo.date.setUnixDate((time_t)(floor(tt/1000.0)));
	tmpmeteo.date.setTimeZone(default_timezone);

	for (size_t jj=0; jj < streamElement->field.size(); jj++){
		const string value = IOUtils::strToUpper( streamElement->field.at(jj)->__item );
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

	vecMeteo.push_back(tmpmeteo);
	tmpmeteo(MeteoData::TSS) = IOUtils::nodata; //if tss has been set, then it needs to be reset manually
	*/
}

double GSNIO::olwr_to_tss(const double& olwr) {
	const double ea = 1.;
	if (olwr == IOUtils::nodata)
		return IOUtils::nodata;

	return pow( olwr / ( ea * Cst::stefan_boltzmann ), 0.25);
}

void GSNIO::readAssimilationData(const Date&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::write2DGrid(const Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GSNIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
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

	double& hs = meteo(MeteoData::HS);
	if (hs != IOUtils::nodata)
		hs /= 100.;
}

size_t GSNIO::data_write(void* buf, size_t size, size_t nmemb, void* userp)
{
	if(userp)
	{
		std::ostream& os = *static_cast<std::ostream*>(userp);
		std::streamsize len = size * nmemb;
		if(os.write(static_cast<char*>(buf), len))
			return len;
	}

	return 0;
}

CURLcode GSNIO::curl_read(const std::string& url_query, std::ostream& os)
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	const string url = endpoint + url_query;
	cout << url << endl;

	if (curl) {
		if(CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &data_write))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FILE, &os))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, GSNIO::http_timeout))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_URL, url.c_str())))
		{
			code = curl_easy_perform(curl);
		}
		curl_easy_cleanup(curl);
	}

	return code;
}

// The GARAGE:
void GSNIO::listSensors(std::vector<std::string>& /*vec_names*/)
{
	/**
	 * Parse through the io.ini file and copy the desired station names STATION#
	 * into vecStationName, if no stations are configured explicitly simply take
	 * all stations that are available in the current GSN instance
	 */

	/*
	_ns1__listVirtualSensorNamesResponse sensor_names;
	_ns1__listVirtualSensorNames sensor_req;

	if (gsn.listVirtualSensorNames(&sensor_req, &sensor_names) == SOAP_OK){
		cerr << "[I] No STATIONS specified... Using all " << sensor_names.virtualSensorName.size() << " sensors available through GSN\n";
		for (size_t ii=0; ii<sensor_names.virtualSensorName.size(); ii++){
			cerr << "\tSTATION" << ii+1 << " = " << sensor_names.virtualSensorName[ii] << "\n";
			vec_names.push_back(sensor_names.virtualSensorName[ii]);
		}
	} else {
		soap_print_fault(&gsn, stderr);
		throw IOException("Error in communication with GSN",AT);
	}
	*/
}

} //namespace




	/*
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
			const size_t nr_meta = metadata.virtualSensorDetails.size();

			if (nr_meta == 0)
				throw IOException("No meta data for sensor " + vecStationName[ii] + ". Is it a valid sensor?", AT);

			for (size_t jj=0; jj<nr_meta; jj++){
				const size_t predicates = metadata.virtualSensorDetails[jj]->addressing->predicates.size();
				for (size_t kk=0; kk<predicates; kk++){
					const string field_name = IOUtils::strToUpper( metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->name );
					const string field_val  = IOUtils::strToUpper( metadata.virtualSensorDetails[jj]->addressing->predicates[kk]->__item );

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
							std::string tmp2;
							IOUtils::convertString(tmp2, field_val);
							if (IOUtils::isNumeric(tmp2)) IOUtils::convertString(slope_azi, field_val);
							else slope_azi=IOUtils::bearing(tmp2);
							info_complete |= 16;
						}
					}
				}

				if (info_complete != 31){
					;//throw IOException("Incomplete meta data (location info) for sensor " + vecStationName[ii], AT);
				}
			}
		} else {
			soap_print_fault(&gsn, stderr);
			throw IOException("Error in communication with GSN while retrieving virtual sensor meta data, see soap error message", AT);
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
	*/



	/*
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
		const string sid = data.queryResult.at(0)->sid;
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
	*/
