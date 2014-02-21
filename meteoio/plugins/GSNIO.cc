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

#include <algorithm>
#include <sstream>
#include <iostream>

#include <curl/curl.h>

using namespace std;

namespace mio {
/**
 * @page gsn GSN
 * @section gsn_format Format
 * This plugin reads meteorological data from GSN (Global Sensor Network, see <a href="http://sourceforge.net/apps/trac/gsn/"> GSN home page</a>)
 * via the RESTful web service. To compile the plugin you need to have the <a href="http://curl.haxx.se/">CURL library</a> with its headers present.
 * @subsection gsn_fields Field mapping
 * The following GSN fields are read from GSN and mapped to MeteoData attributes:
 * <center><table border="0">
 * <tr><td>
 * <table border="1">
 * <tr><th>GSN attribute</th><th>MeteoData field</th></tr>
 * <tr><td>RELATIVE_HUMIDITY or AIR_HUMID</td><td>MeteoData::RH</td></tr>
 * <tr><td>AIR_TEMPERATURE or AIR_TEMP</td><td>MeteoData::TA</td></tr>
 * <tr><td>WIND_DIRECTION</td><td>MeteoData::DW</td></tr>
 * <tr><td>WIND_SPEED_MAX</td><td>MeteoData::VW_MAX</td></tr>
 * <tr><td>WIND_SPEED_SCALAR_AV or WIND_SPEED</td><td>MeteoData::VW</td></tr>
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
 * Please keep in mind that the names in GSN have currently not been standardized. This means that any sensor that does
 * not use the above names will not be properly supported (fields will not be missing but might appear under a different name)!
 *
 * @section gsn_units Units
 * The units of measurements are sometimes listed in the response headers, they are then parsed by the plugin and if known,
 * like <b>°C</b> or <b>\%</b>, offsets and multipliers are set to convert the data to MKSA
 *
 * Otherwise the units are assumed to be the following:
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
 * - GSN_URL: The URL of the RESTful web service e.g. http://planetdata.epfl.ch:22001/rest
 * - GSN_USER: The username to access the service
 * - GSN_PASS: The password to authenticate the USER
 * - STATION#: station code for the given number #, e. g. la_fouly_1034 (case sensitive!)
 *
 * If no STATION keys are given, the full list of ALL stations available to the user in GSN will be used!
 * This may result in a long download.
 *
 * @code
 * METEO	= GSN
 * GSN_URL	= http://montblanc.slf.ch:22001/rest
 * GSN_USER	= mylogin
 * GSN_PASS	= mypasswd
 * STATION1	= wind_tunnel_meteo
 * @endcode
 *
 */

const int GSNIO::http_timeout = 60; // seconds until connect time out for libcurl
const std::string GSNIO::sensors_endpoint = "sensors";
const std::string GSNIO::null_string = "null";

GSNIO::GSNIO(const std::string& configfile)
      : cfg(configfile), vecStationName(), multiplier(), offset(), vecMeta(), vecAllMeta(), coordin(),
        coordinparam(), coordout(), coordoutparam(), endpoint(), userid(), passwd(), default_timezone(1.)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::GSNIO(const Config& cfgreader)
      : cfg(cfgreader), vecStationName(), multiplier(), offset(), vecMeta(), vecAllMeta(), coordin(),
        coordinparam(), coordout(), coordoutparam(), endpoint(), userid(), passwd(), default_timezone(1.)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();
}

GSNIO::~GSNIO() throw(){}

void GSNIO::initGSNConnection() {
	curl_global_init(CURL_GLOBAL_ALL);

	default_timezone = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", default_timezone, IOUtils::nothrow);

	cfg.getValue("GSN_URL", "Input", endpoint, IOUtils::nothrow);
	if (!endpoint.empty()){
		if (*endpoint.rbegin() != '/') endpoint += "/";
		cerr << "[i] Using GSN Endpoint: " << endpoint << endl;
	}

	cfg.getValue("GSN_USER", "Input", userid);
	cfg.getValue("GSN_PASS", "Input", passwd);
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
		cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName

	//Get Meta Data for all stations first
	getAllStations();

	if (!vecStationName.empty()) { //if the user has specified a subset of stations
		for (size_t ii=0; ii<vecStationName.size(); ii++) {
			for (size_t jj=0; jj<vecAllMeta.size(); jj++) {
				if (vecAllMeta[jj].stationID == vecStationName[ii]) {
					vecMeta.push_back(vecAllMeta[jj]);
				}
			}

			if (vecMeta.size() != (ii+1)) { // could not find station in list of available stations
				throw NoAvailableDataException("Could not retrieve meta data for station " + vecStationName[ii], AT);
			}
		}
	} else { //otherwise use all available stations
		vecMeta = vecAllMeta;
	}
}

void GSNIO::save_station(const std::string& id, const std::string& name, const double& lat, const double& lon,
                         const double& alt, const double& slope_angle, const double& slope_azi)
{
	Coords current_coord(coordin, coordinparam);
	current_coord.setLatLon(lat, lon, alt);
	StationData sd(current_coord, id, name);

	if (slope_angle != IOUtils::nodata) {
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
	 * and which are accessible for the current user (see Input::GSN_USER)
	 */
	const string vsname_str("# vsname:");
	const string altitude_str("# altitude:");
	const string longitude_str("# longitude:");
	const string latitude_str("# latitude:");
	const string slope_str("# slope:");
	const string exposition_str("# exposition:");
	const string name_str("# name:");

	stringstream ss;
	string line("");

	vecAllMeta.clear();

	if (curl_read(sensors_endpoint + "?username=" + userid + "&password=" + passwd, ss)) {
		string name(""), id(""), azi("");
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
	const string fields_str("# fields:");
	const string units_str("# units:");

	const string request = sensors_endpoint + "/" + vecMeta[stationindex].stationID + "?from=" + dateStart.toString(Date::ISO) + ":00"
	                 + "&to=" + dateEnd.toString(Date::ISO) + ":00" + "&username=" + userid + "&password=" + passwd;

	stringstream ss;

	if (curl_read(request, ss)) {
		vector<size_t> index;
		bool olwr_present = false;

		MeteoData tmpmeteo;
		tmpmeteo.meta = vecMeta.at(stationindex);

		string line(""), fields(""), units("");
		while (getline(ss, line)) { //parse header section

			if (line.size() && (line[0] != '#')) break;

			if (!line.compare(0, fields_str.size(), fields_str)) {
				fields = line.substr(fields_str.size());
			} else if (!line.compare(0, units_str.size(), units_str)) {
				units = line.substr(units_str.size()) + " "; // the extra space is important if no units are specified
			}
		}

		if (units.empty() || fields.empty()) {
			throw InvalidFormatException("Invalid header for station " + tmpmeteo.meta.stationID, AT);
		}

		map_parameters(fields, units, tmpmeteo, index);
		olwr_present = tmpmeteo.param_exists("OLWR");

		do { //parse data section, the first line should already be buffered
			parse_streamElement(line, index, olwr_present, vecMeteo, tmpmeteo);
		} while (getline(ss, line));
	} else {
		throw IOException("Could not retrieve data for station " + vecMeta[stationindex].stationID, AT);
	}
}

void GSNIO::map_parameters(const std::string& fields, const std::string& units, MeteoData& md, std::vector<size_t>& index)
{
	vector<string> field, unit;
	size_t timestamp_field = IOUtils::npos;
	multiplier.clear();
	offset.clear();

	IOUtils::readLineToVec(fields, field, ',');
	IOUtils::readLineToVec(units, unit, ',');

	if ((field.size() != unit.size()) || (field.size() < 2)) {
		throw InvalidFormatException("Fields and units are inconsistent for station " + md.meta.stationID, AT);
	}

	for (size_t ii=0; ii<field.size(); ii++) {
		const string field_name( IOUtils::strToUpper(field[ii]) );

		if (field_name == "RELATIVE_HUMIDITY" || field_name == "RH" || field_name == "AIR_HUMID" || field_name == "REL_HUMIDITY") {
			index.push_back(MeteoData::RH);
		} else if (field_name == "AIR_TEMPERATURE" || field_name == "TA" || field_name == "AIR_TEMP") {
			index.push_back(MeteoData::TA);
		} else if (field_name == "WIND_DIRECTION" || field_name == "DW") {
			index.push_back(MeteoData::DW);
		} else if (field_name == "WIND_SPEED_MAX" || field_name == "VW_MAX") {
			index.push_back(MeteoData::VW_MAX);
		} else if (field_name == "WIND_SPEED_SCALAR_AV" || field_name == "VW" || field_name == "WIND_SPEED") {
			index.push_back(MeteoData::VW);
		} else if (field_name == "INCOMING_SHORTWAVE_RADIATION" || field_name == "ISWR" || field_name == "SOLAR_RAD") {
			index.push_back(MeteoData::ISWR);
		} else if (field_name == "INCOMING_LONGWAVE_RADIATION" || field_name == "ILWR") {
			index.push_back(MeteoData::ILWR);
		} else if (field_name == "OUTGOING_SHORTWAVE_RADIATION" || field_name == "RSWR") {
			index.push_back(MeteoData::RSWR);
		} else if (field_name == "OUTGOING_LONGWAVE_RADIATION" || field_name == "RLWR") { //is used to calculate TSS
			md.addParameter("OLWR");
			index.push_back(md.getParameterIndex("OLWR"));
		} else if (field_name == "SNOW_HEIGHT" || field_name == "HS1") {
			index.push_back(MeteoData::HS);
		} else if (field_name == "RAIN_METER" || field_name == "PINT") {
			index.push_back(MeteoData::HNW);
		} else if (field_name == "SURFACE_TEMP" || field_name == "TSS" || field_name == "SNOW_SURFACE_TEMPERATURE") {
			index.push_back(MeteoData::TSS);
		} else if (field_name == "ATM_PRESSURE" || field_name == "P") {
			index.push_back(MeteoData::P);
		} else if (field_name == "TIMESTAMP") {
			timestamp_field = ii;
			index.push_back(IOUtils::npos);
		} else if (field_name == "TIME") {
			index.push_back(IOUtils::npos);
		} else { //this is an extra parameter
			md.addParameter(field_name);
			const size_t parindex = md.getParameterIndex(field_name);
			index.push_back(parindex);

			//For the parameters unknown to MeteoIO we can store the units qualification °C, %, etc
			//and make it possible for the values to be converted to MKSA in the convertUnits procedure
			string name( unit[ii] );
			IOUtils::trim(name);

			if (name == "%") {
				multiplier[parindex] = 0.01;
			} else if (name.size() == 2 && (int)((unsigned char)name[0]) == 176 && name[1] == 'C') { //in °C, UTF8
				offset[parindex] = 273.15;
			}
		}
	}

	if (timestamp_field != IOUtils::npos) { //store timestamp index at index[0]
		index[0] = timestamp_field;
	} else {
		throw InvalidFormatException("No timestamp field for station " + md.meta.stationID, AT);
	}
}

void GSNIO::parse_streamElement(const std::string& line, const std::vector<size_t>& index, const bool& olwr_present, std::vector<MeteoData>& vecMeteo, MeteoData& tmpmeteo) const
{
	static vector<string> data;
	static double timestamp;
	static const size_t timestamp_index = index[0];

	const size_t size = IOUtils::readLineToVec(line, data, ',');
	if (size < 2) return; // Malformed for sure, retire gracefully, no exception thrown

	//The timestamp index is stored in index[0]
	IOUtils::convertString(timestamp, data[timestamp_index]);
	tmpmeteo.date.setUnixDate((time_t)(floor(timestamp/1000.0)));
	tmpmeteo.date.setTimeZone(default_timezone);

	for (size_t jj=2; jj<size; jj++) {
		const string& value = data[jj];
		if (value != GSNIO::null_string){
			IOUtils::convertString(tmpmeteo(index[jj]), value);
		} else {
			tmpmeteo(index[jj]) = IOUtils::nodata;
		}
	}

	convertUnits(tmpmeteo);
	if ((olwr_present) && (tmpmeteo(MeteoData::TSS) == IOUtils::nodata))
		tmpmeteo(MeteoData::TSS) = olwr_to_tss(tmpmeteo("OLWR"));

	vecMeteo.push_back(tmpmeteo);
	tmpmeteo(MeteoData::TSS) = IOUtils::nodata; //if tss has been set, then it needs to be reset manually
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

void GSNIO::convertUnits(MeteoData& meteo) const
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

	// For all parameters that have either an offset or an multiplier to bring to MKSA
	map<size_t, double>::const_iterator it;
	for (it = multiplier.begin(); it != multiplier.end(); it++) {
		double& tmp = meteo(it->first);
		if (tmp != IOUtils::nodata) tmp *= it->second;
	}

	for (it = offset.begin(); it != offset.end(); it++) {
		double& tmp = meteo(it->first);
		if (tmp != IOUtils::nodata) tmp += it->second;
	}
}

size_t GSNIO::data_write(void* buf, size_t size, size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os = *static_cast<ostream*>(userp);
		const streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool GSNIO::curl_read(const std::string& url_query, std::ostream& os)
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	const string url = endpoint + url_query;

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

	if(code!=CURLE_OK)
		std::cout << "[E] " << curl_easy_strerror(code) << "\n";

	return (code==CURLE_OK);
}

} //namespace
