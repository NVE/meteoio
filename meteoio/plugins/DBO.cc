/***********************************************************************************/
/*  Copyright 2017 SLF                                                                                                                                */
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
#include <meteoio/plugins/DBO.h>

#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/meteoLaws/Meteoconst.h>

#include <algorithm>
#include <sstream>
#include <iostream>

#include <curl/curl.h>
#include <meteoio/plugins/picojson.h>

using namespace std;

namespace mio {
/**
 * @page dbo DBO
 * @section dbo_format Format
 * This plugin reads meteorological data from DBO
 * via the RESTful web service. To compile the plugin you need to have the <a href="http://curl.haxx.se/">CURL library</a> with its headers present.
 *
 * @section dbo_keywords Keywords
 * This plugin uses the following keywords:
 * - DBO_URL: The URL of the RESTful web service e.g.http://developwis.wsl.ch:8730/osper-api
 * - DBO_USER: The username to access the service (optional)
 * - DBO_PASS: The password to authenticate the USER (optional)
 * - STATION#: station code for the given station, prefixed by the network it belongs ot (for example: IMIS::SLF2)
 * - DBO_TIMEOUT: timeout (in seconds) for the connection to the server (default: 60s)
 * - DBO_DEBUG: print the full requests/answers from the server when something does not work as expected
 *
 * @code
 * METEO	= DBO
 * DBO_URL	= http://developwis.wsl.ch:8730/osper-api
 * DBO_USER	= mylogin
 * DBO_PASS	= mypasswd
 * STATION1	= wind_tunnel_meteo
 * @endcode
 *
 * @section dbo_dependencies Picojson
 * This plugin relies on <A HREF="https://github.com/kazuho/picojson/">picojson</A> for reading and parsing
 * <A HREF="https://en.wikipedia.org/wiki/JSON">JSON</A> data. Picojson is released under a
 * <A HREF="https://opensource.org/licenses/BSD-2-Clause">2-Clause BSD License</A>. Please find here below
 * the full license agreement for picojson:
 *
 * @code
 * Copyright 2009-2010 Cybozu Labs, Inc.
 * Copyright 2011-2014 Kazuho Oku
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * @endcode
 */

//we keep it as a simple function in order to avoid exposing picojson stuff in the header
void printJSON(const picojson::value& v, const unsigned int& depth)
{
	if (v.is<picojson::null>()) {
		for (unsigned int jj=0; jj<depth; jj++) std::cout << "\t";
		std::cout << "NULL\n";
		return;
	}

	if (v.is<picojson::object>()) {
		const picojson::value::object& obj = v.get<picojson::object>();
		for (picojson::value::object::const_iterator ii = obj.begin(); ii != obj.end(); ++ii) {
			for (unsigned int jj=0; jj<depth; jj++) std::cout << "\t";
			std::cout << ii->first << "\n";
			printJSON(ii->second, depth+1);
		}
	} else if (v.is<std::string>()){
		for (unsigned int jj=0; jj<depth; jj++) std::cout << "\t";
		std::cout << v.get<std::string>() << "\n";
	} else if (v.is<double>()){
		for (unsigned int jj=0; jj<depth; jj++) std::cout << "\t";
		std::cout << v.get<double>() << "\n";
	} else if (v.is<bool>()){
		for (unsigned int jj=0; jj<depth; jj++) std::cout << "\t";
		std::cout << std::boolalpha << v.get<bool>() << "\n";
	} else if (v.is<picojson::array>()){ //ie vector<picojson::value>
		for (unsigned int jj=0; jj<depth; jj++) std::cout << "\t";
		const picojson::array& array = v.get<picojson::array>();
		std::cout << "array " << array.size() << "\n";
		for (size_t jj=0; jj<array.size(); jj++)
			printJSON(array[jj], depth+1);
	}
}

picojson::value goToJSONPath(const std::string& path, picojson::value& v)
{
	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	if (v.is<picojson::object>()) {
		picojson::value::object& obj = v.get<picojson::object>();
		for (std::map<std::string,picojson::value>::iterator it = obj.begin(); it != obj.end(); ++it) {
			if (it->first==local_path) {
				if (!remaining_path.empty())
					goToJSONPath(remaining_path, it->second);
				else
					return it->second;
			}
		}
	}

	return picojson::value();
}

void JSONQuery(const std::string& path, picojson::value& v, std::vector<picojson::value>& results)
{
	if (v.is<picojson::null>()) return;

	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	if (v.is<picojson::object>()) {
		picojson::value::object& obj = v.get<picojson::object>();
		for (std::map<std::string,picojson::value>::iterator it = obj.begin(); it != obj.end(); ++it) {
			if (it->first==local_path) {
				if (!remaining_path.empty()) {
					 if (it->second.is<picojson::array>()){ //ie vector<picojson::value>
						picojson::array& array = it->second.get<picojson::array>();
						for (size_t jj=0; jj<array.size(); jj++)
							JSONQuery(remaining_path, array[jj], results);
					} else
						JSONQuery(remaining_path, it->second, results);
				} else {
					results.push_back( it->second );
				}
			}
		}
	}
}

std::string getString(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);
	if (!results.empty()) {
		if (! results.front().is<picojson::null>() &&  results.front().is<std::string>()) return  results.front().get<std::string>();
	}

	return std::string();
}

std::vector<std::string> getStrings(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);

	std::vector<std::string> vecString;
	for (size_t ii=0; ii<results.size(); ii++) {
		 if (results[ii].is<picojson::array>()){ //ie vector<picojson::value>
			const picojson::array& array = results[ii].get<picojson::array>();
			for (size_t jj=0; jj<array.size(); jj++) {
				if (! array[jj].is<picojson::null>() &&  array[jj].is<std::string>()) vecString.push_back( array[jj].get<std::string>() );
			}
		} else
			if (! results[ii].is<picojson::null>() &&  results[ii].is<std::string>()) vecString.push_back( results[ii].get<std::string>() );
	}

	return vecString;
}

double getDouble(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);
	if (!results.empty()) {
		if (! results.front().is<picojson::null>() &&  results.front().is<double>()) return  results.front().get<double>();
	}

	return IOUtils::nodata;
}

std::vector<double> getDoubles(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);

	std::vector<double> vecDouble;
	for (size_t ii=0; ii<results.size(); ii++) {
		 if (results[ii].is<picojson::array>()){ //ie vector<picojson::value>
			const picojson::array& array = results[ii].get<picojson::array>();
			for (size_t jj=0; jj<array.size(); jj++) {
				if (! array[jj].is<picojson::null>() &&  array[jj].is<double>()) vecDouble.push_back( array[jj].get<double>() );
			}
		} else
			if (! results[ii].is<picojson::null>() &&  results[ii].is<double>()) vecDouble.push_back( results[ii].get<double>() );
	}

	return vecDouble;
}

//converts C to Kelvin, converts RH to [0,1], HS to m
double convertUnits(const MeteoData::Parameters& param, const double& value)
{
	switch (param) {
		case MeteoData::TA: case MeteoData::TSG: case MeteoData::TSS:
			return IOUtils::C_TO_K(value);
		case MeteoData::RH: case MeteoData::HS:
			return value/100.;
		default:
			return value;
	}
}

bool parseTsPoint(const picojson::value& v, Date& datum, double& value)
{
	if (!v.is<picojson::array>()) return false;

	const picojson::array& array = v.get<picojson::array>();
	if (array.size()!=2) return false;

	if (array[0].is<std::string>())
		IOUtils::convertString(datum, array[0].get<std::string>(), 0.); //HACK: hard-coding GMT
	else
		return false;

	if (array[1].is<double>())
		value = array[1].get<double>();
	else
		return false;

	return true;
}

void parseTimeSerie(const std::string& tsID, const MeteoData::Parameters& param, const StationData& sd, picojson::value& v, std::vector<MeteoData>& vecMeteo)
{
	picojson::value ts( goToJSONPath("$.measurements", v) );
	if (!ts.is<picojson::array>())
		throw InvalidFormatException("Could not parse timeserie "+tsID, AT);
	const picojson::array& vecRaw = ts.get<picojson::array>();

	if (vecMeteo.empty()) {
		vecMeteo.reserve(vecRaw.size());
		for (size_t ii=0; ii<vecRaw.size(); ii++) {
			Date datum;
			double value;
			if (!parseTsPoint(vecRaw[ii], datum, value)) {
				std::ostringstream ss; ss << "Error parsing element " << ii << " of timeserie " << tsID;
				throw InvalidFormatException(ss.str(), AT);
			}

			MeteoData md(datum, sd);
			md(param) = convertUnits(param, value);
			vecMeteo.push_back( md );
		}
	} else {
		for (size_t ii=0; ii<vecRaw.size(); ii++) {
			Date datum;
			double value;
			if (!parseTsPoint(vecRaw[ii], datum, value)) {
				std::ostringstream ss; ss << "Error parsing element " << ii << " of timeserie " << tsID;
				throw InvalidFormatException(ss.str(), AT);
			}

			/*MeteoData md(datum, sd);
			md(param) = convertUnits(param, value);
			vecMeteo.push_back( md );*/
		}
	}
}

unsigned int parseInterval(const std::string& interval_str)
{
	unsigned int hour, minute, second;
	if (sscanf(interval_str.c_str(), "%u:%u:%u", &hour, &minute, &second) == 3) {
		return (hour*3600 + minute*60 + second);
	} else
		throw ConversionFailedException("Could not read aggregation interval '"+interval_str+"'", AT);
}

std::map<std::string, std::vector<DBO::tsMeta> > getTsProperties(picojson::value& v)
{
	std::map<std::string, std::vector<DBO::tsMeta> > tsMap;
	std::vector<picojson::value> results;
	JSONQuery("$.properties.timeseries", v, results);

	for (size_t ii=0; ii<results.size(); ii++) {
		 if (results[ii].is<picojson::array>()){
			const picojson::array& array = results[ii].get<picojson::array>();
			for (size_t jj=0; jj<array.size(); jj++) {
				if (! array[jj].is<picojson::null>()) {
					std::string code, device_code, agg_type;
					double ts_id;
					unsigned int interval;
					Date since, until;

					const picojson::value::object& obj = array[jj].get<picojson::object>();
					for (picojson::value::object::const_iterator it = obj.begin(); it != obj.end(); ++it) {
						if (it->first=="code" && it->second.is<std::string>()) code = it->second.get<std::string>();
						if (it->first=="deviceCode" && it->second.is<std::string>()) device_code = it->second.get<std::string>();
						if (it->first=="id" && it->second.is<double>()) ts_id = it->second.get<double>();
						if (it->first=="since" && it->second.is<std::string>()) IOUtils::convertString(since, it->second.get<std::string>(), 0.);
						if (it->first=="until" && it->second.is<std::string>()) IOUtils::convertString(until, it->second.get<std::string>(), 0.);
						if (it->first=="aggregationType" && it->second.is<std::string>()) agg_type = it->second.get<std::string>();
						if (it->first=="aggregationInterval" && it->second.is<std::string>()) interval = parseInterval(it->second.get<std::string>());
					}

					if (device_code=="BATTERY" || device_code=="LOGGER") break;
					if (agg_type=="SD") break; //we don't care about standard deviation anyway

					std::string param_str( code.substr(0, code.find('_')) );
					IOUtils::toUpper( param_str );
					tsMap[param_str].push_back( DBO::tsMeta(since, until, agg_type, ts_id, interval) );
				}
			}
		}
	}

	return tsMap;
}

/*************************************************************************************************/
const int DBO::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string DBO::metadata_endpoint = "/osper-api/osper/stations/";
const std::string DBO::data_endpoint = "/osper-api/osper/timeseries/";
const std::string DBO::null_string = "null";

DBO::DBO(const std::string& configfile)
      : cfg(configfile), vecStationName(), vecMeta(), vecTsMeta(),
        coordin(), coordinparam(), coordout(), coordoutparam(),
        endpoint(), default_timezone(1.),
        http_timeout(http_timeout_dflt), dbo_debug(false)
{
	initDBOConnection();
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

DBO::DBO(const Config& cfgreader)
      : cfg(cfgreader), vecStationName(), vecMeta(), vecTsMeta(),
        coordin(), coordinparam(), coordout(), coordoutparam(),
        endpoint(), default_timezone(1.),
        http_timeout(http_timeout_dflt), dbo_debug(false)
{
	initDBOConnection();
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

void DBO::initDBOConnection() {
	curl_global_init(CURL_GLOBAL_ALL);

	cfg.getValue("DBO_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Input", default_timezone, IOUtils::nothrow);

	cfg.getValue("DBO_URL", "Input", endpoint);
	if (*endpoint.rbegin() != '/') endpoint += "/";
	cerr << "[i] Using DBO URL: " << endpoint << endl;

	cfg.getValue("DBO_DEBUG", "INPUT", dbo_debug, IOUtils::nothrow);
}

void DBO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();
	if (vecMeta.empty()) fillStationMeta();
	vecStation = vecMeta;
}

void DBO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                          std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	if (vecMeta.empty()) fillStationMeta();

	vecMeteo.resize(vecMeta.size());
	for(size_t ii=0; ii<vecMeta.size(); ii++)
		readData(dateStart, dateEnd, vecMeteo[ii], ii);
}

void DBO::fillStationMeta()
{
	vecMeta.clear();
	vecTsMeta.resize( vecStationName.size() );

	for(size_t ii=0; ii<vecStationName.size(); ii++) {
		std::string station_id( vecStationName[ii] );
		if (station_id.find(':')==std::string::npos) station_id = "IMIS::" + station_id;
		const std::string request( metadata_endpoint + "name=" + IOUtils::strToLower( station_id ) );

		std::stringstream ss;
		if (curl_read(request, ss)) {
			if (ss.str().empty()) throw UnknownValueException("Station not found: '"+station_id+"'", AT);

			picojson::value v;
			const std::string err( picojson::parse(v, ss.str()) );
			if (!err.empty()) throw IOException("Error while parsing JSON: "+err, AT);

			const std::vector<double> coordinates( getDoubles("$.geometry.coordinates", v) );
			if (coordinates.size()!=3) throw InvalidFormatException("Wrong coordinates specification!", AT);

			Coords position(coordin, coordinparam);
			position.setLatLon(coordinates[1], coordinates[0], coordinates[2]);
			const StationData sd(position, getString("$.properties.name", v), getString("$.properties.locationName", v));
			vecMeta.push_back( sd );

			//select proper time series
			vecTsMeta[ii] = getTsProperties(v);
		} else {
			if (dbo_debug)
				std::cout << "****\nRequest: " << request << "\n****\n";
			throw IOException("Could not retrieve data for station " + station_id, AT);
		}
	}

	//debug info
	/*for(size_t ii=0; ii<vecStationName.size(); ii++) {
		for (std::map<std::string, std::vector<DBO::tsMeta> >::iterator it = vecTsMeta[ii].begin(); it != vecTsMeta[ii].end(); ++it) {
			for(size_t jj=0; jj<it->second.size(); jj++)
				std::cout << it->first << " " << it->second[jj].toString() << "\n";
		}
	}*/
}

//read all data for the given station
void DBO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	//TODO: for station stationindex, loop over the timeseries that cover [dateStart, dateEnd] for the current station
	//vecTsMeta[ stationindex ]

	const unsigned int ts_id = 39;
	readTimeSerie(ts_id, dateStart, dateEnd, vecMeta[stationindex], vecMeteo);
}

void DBO::readTimeSerie(const unsigned int& ts_id, const Date& dateStart, const Date& dateEnd, const StationData& sd, std::vector<MeteoData>& vecMeteo)
{
	std::ostringstream ss_ID; ss_ID << ts_id;
	const std::string base_url( data_endpoint + ss_ID.str() );
	//const std::string period( "?from=" + dateStart.toString(Date::ISO_TZ) + "&until=" + dateEnd.toString(Date::ISO_TZ) );
	const std::string period( string("?from=2016-11-17T13:00Z") + "&until=" + "2017-01-05T13:00Z" );
	const std::string request( base_url + period );

	std::stringstream ss;
	if (curl_read(request, ss)) {
		if (ss.str().empty()) throw UnknownValueException("Timeserie not found: '"+ss_ID.str()+"'", AT);
		picojson::value v;
		const std::string err( picojson::parse(v, ss.str()) );
		if (!err.empty()) throw IOException("Error while parsing JSON: "+err, AT);

		parseTimeSerie(ss_ID.str(), MeteoData::RH, sd, v, vecMeteo);
	} else {
		if (dbo_debug)
			std::cout << "****\nRequest: " << request << "\n****\n";
		throw IOException("Could not retrieve data for timeseries " + ss_ID.str(), AT);
	}
}

size_t DBO::data_write(void* buf, size_t size, size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os = *static_cast<ostream*>(userp);
		const std::streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool DBO::curl_read(const std::string& url_query, std::ostream& os)
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	const std::string url( endpoint + url_query );

	if (curl) {
		if (CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &data_write))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FILE, &os))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, DBO::http_timeout))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_URL, url.c_str())))
		{
			code = curl_easy_perform(curl);
		}
		curl_easy_cleanup(curl);
	}

	if (code!=CURLE_OK) {
		if (dbo_debug)
			std::cout << "****\nRequest: " << url_query << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}

	return (code==CURLE_OK);
}

} //namespace
