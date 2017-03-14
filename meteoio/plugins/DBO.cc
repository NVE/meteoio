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
			//results.reserve( results.size()+array.size() ); //most of the time, we will come here with an empty vector
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
	if (value==IOUtils::nodata) return value;

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
		IOUtils::convertString(datum, array[0].get<std::string>(), 0.); //hard-coding GMT
	else
		return false;

	if (array[1].is<double>())
		value = array[1].get<double>();
	else {
		if (!array[1].is<picojson::null>()) return false;
		value = IOUtils::nodata;
		return true;
	}

	return true;
}

const std::vector<DBO::tsData> parseTimeSerie(const std::string& tsID, picojson::value& v)
{
	picojson::value ts( goToJSONPath("$.measurements", v) );
	if (!ts.is<picojson::array>())
		throw InvalidFormatException("Could not parse timeserie "+tsID, AT);
	const picojson::array& vecRaw = ts.get<picojson::array>();

	if (vecRaw.empty()) return std::vector<DBO::tsData>();

	std::vector<DBO::tsData> vecData( vecRaw.size() );
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		double value;
		Date datum;
		if (!parseTsPoint(vecRaw[ii], datum, value))  {
			printJSON(vecRaw[ii], 0);
			std::ostringstream ss; ss << "Error parsing element " << ii << " of timeserie " << tsID;
			throw InvalidFormatException(ss.str(), AT);
		}
		vecData[ii] = DBO::tsData(datum, value);
	}

	return vecData;
}

unsigned int parseInterval(const std::string& interval_str)
{
	unsigned int hour, minute, second;
	if (sscanf(interval_str.c_str(), "%uMIN", &minute) == 1) {
		return (minute*60);
	} if (sscanf(interval_str.c_str(), "%u:%u:%u", &hour, &minute, &second) == 3) {
		return (hour*3600 + minute*60 + second);
	} else
		throw ConversionFailedException("Could not read measure interval '"+interval_str+"'", AT);
}

std::map<size_t, DBO::tsMeta> getTsProperties(picojson::value& v)
{
	std::map<size_t, DBO::tsMeta> tsMap;
	std::vector<picojson::value> results;
	JSONQuery("$.properties.timeseries", v, results);

	for (size_t ii=0; ii<results.size(); ii++) {
		 if (results[ii].is<picojson::array>()){
			const picojson::array& array = results[ii].get<picojson::array>();
			for (size_t jj=0; jj<array.size(); jj++) {
				if (! array[jj].is<picojson::null>()) {
					std::string code, device_code, agg_type;
					double id = -1.;
					unsigned int interval = 0;
					Date since, until;

					const picojson::value::object& obj = array[jj].get<picojson::object>();
					for (picojson::value::object::const_iterator it = obj.begin(); it != obj.end(); ++it) {
						if (it->first=="code" && it->second.is<std::string>()) code = it->second.get<std::string>();
						if (it->first=="deviceCode" && it->second.is<std::string>()) device_code = it->second.get<std::string>();
						if (it->first=="id" && it->second.is<double>()) id = it->second.get<double>();
						if (it->first=="since" && it->second.is<std::string>()) IOUtils::convertString(since, it->second.get<std::string>(), 0.);
						if (it->first=="until" && it->second.is<std::string>()) IOUtils::convertString(until, it->second.get<std::string>(), 0.);
						if (it->first=="aggregationType" && it->second.is<std::string>()) agg_type = it->second.get<std::string>();
						if (it->first=="measureInterval" && it->second.is<std::string>()) interval = parseInterval(it->second.get<std::string>());
					}

					if (device_code=="BATTERY" || device_code=="LOGGER") break;
					if (agg_type=="SD") break; //we don't care about standard deviation anyway
					if (id==-1.) break; //no id was provided

					const std::string param( IOUtils::strToUpper( code.substr(0, code.find('_')) ) );
					tsMap[ static_cast<size_t>(id) ] = DBO::tsMeta(param, since, until, agg_type, interval);
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

void DBO::initDBOConnection()
{
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

			//parse and store the time series belonging to this station
			vecTsMeta[ii] = getTsProperties(v);
		} else {
			if (dbo_debug)
				std::cout << "****\nRequest: " << request << "\n****\n";
			throw IOException("Could not retrieve data for station " + station_id, AT);
		}
	}
}

bool DBO::getParameter(const std::string& param_str, const std::string& agg_type, MeteoData::Parameters &param)
{
	if (agg_type=="SD") return false; //we don't care about standard deviation anyway

	if (param_str=="P") param = MeteoData::P;
	else if (param_str=="TA") param = MeteoData::TA;
	else if (param_str=="RH") param = MeteoData::RH;
	else if (param_str=="TSG") param = MeteoData::TSG;
	else if (param_str=="TSS") param = MeteoData::TSS;
	else if (param_str=="HS") param = MeteoData::HS;
	else if (param_str=="VW" && agg_type=="MAX") param = MeteoData::VW_MAX;
	else if (param_str=="VW") param = MeteoData::VW;
	else if (param_str=="DW") param = MeteoData::DW;
	else if (param_str=="RSWR") param = MeteoData::RSWR;
	else if (param_str=="ISWR") param = MeteoData::ISWR;
	else if (param_str=="ILWR") param = MeteoData::ILWR;
	else if (param_str=="RRI") param = MeteoData::PSUM;
	else return false;

	return true;
}

//read all data for the given station
void DBO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	const std::string Start( dateStart.toString(Date::ISO_Z) );
	const std::string End( dateEnd.toString(Date::ISO_Z) );

	//for each parameter, a vector of tsID that should be used
	std::map<MeteoData::Parameters, std::vector<size_t> > mapParams;

	//for station stationindex, loop over the timeseries that cover [Start, End] for the current station
	for (std::map<size_t, DBO::tsMeta>::const_iterator it = vecTsMeta[stationindex].begin(); it != vecTsMeta[stationindex].end(); ++it) {
		std::cerr << it->first << " " << it->second.toString() << "\n";
		MeteoData::Parameters param;
		if (getParameter(it->second.param, it->second.agg_type, param)==false) continue; //unrecognized parameter

		if (it->second.interval>3600 && it->second.interval<600) continue; //HACK for now, to keep things simpler for IMIS
		if (!it->second.since.isUndef() && it->second.since>dateEnd) continue; //this TS does not contain our period of interest
		if (!it->second.until.isUndef() && it->second.until<dateStart) continue; //this TS does not contain our period of interest

		mapParams[param].push_back( it->first );
	}

	//second pass: we try a finer selection when multiple TS cover the same parameter and get the suitable TS
	for (std::map<MeteoData::Parameters, std::vector<size_t> >::const_iterator it = mapParams.begin(); it != mapParams.end(); ++it) {
		const MeteoData::Parameters param = it->first;
		for(size_t ii=0; ii<it->second.size(); ii++) {
			//std::cerr << "Reading tsID " << it->second[ii] << " for " << MeteoData::getParameterName(param) << "\n";
			readTimeSerie(it->second[ii], param, Start, End, vecMeta[stationindex], vecMeteo);
		}
	}
}

//dateStart and dateEnd should already be GMT
void DBO::readTimeSerie(const size_t& ts_id, const MeteoData::Parameters& param, const std::string& Start, const std::string& End, const StationData& sd, std::vector<MeteoData>& vecMeteo)
{
	std::ostringstream ss_ID; ss_ID << ts_id;
	const std::string base_url( data_endpoint + ss_ID.str() );
	const std::string request( base_url + "?from=" + Start + "&until=" + End );

	std::stringstream ss;
	if (curl_read(request, ss)) {
		if (ss.str().empty()) throw UnknownValueException("Timeseries not found: '"+ss_ID.str()+"'", AT);
		picojson::value v;
		const std::string err( picojson::parse(v, ss.str()) );
		if (!err.empty()) throw IOException("Error while parsing JSON: "+ss.str(), AT);

		const std::vector<DBO::tsData> vecData( parseTimeSerie(ss_ID.str(), v) );
		mergeTimeSeries(param, vecData, sd, vecMeteo);
	} else {
		if (dbo_debug)
			std::cout << "****\nRequest: " << request << "\n****\n";
		throw IOException("Could not retrieve data for timeseries " + ss_ID.str(), AT);
	}
}

void DBO::mergeTimeSeries(const MeteoData::Parameters& param, const std::vector<DBO::tsData>& vecData, const StationData& sd, std::vector<MeteoData>& vecMeteo)
{
	if (vecData.empty()) return;

	if (vecMeteo.empty()) { //easy case: the initial vector is empty
		vecMeteo.reserve(vecData.size());
		for (size_t ii=0; ii<vecData.size(); ii++) {
			MeteoData md(vecData[ii].date, sd);
			md(param) = convertUnits(param, vecData[ii].val);
			vecMeteo.push_back( md );
		}
	} else {
		size_t vecM_start = 0; //the index in vecRaw that matches the original start of vecMeteo
		size_t vecM_end = 0; //the index in vecRaw that matches the original end of vecMeteo
		const MeteoData md_pattern(Date(), sd); //This assumes that the station is not moving!

		//filling data before vecMeteo
		if (vecData.front().date<vecMeteo.front().date) {
			const Date start_date( vecMeteo.front().date );
			vecM_start = vecData.size(); //if no overlap is found, take all vecData
			for(size_t ii=0; ii<vecData.size(); ii++) { //find the range of elements to add
				if (vecData[ii].date>=start_date) {
					vecM_start = ii;
					break;
				}
			}

			vecMeteo.insert(vecMeteo.begin(), vecM_start, md_pattern);
			for (size_t ii=0; ii<vecM_start; ii++) {
				vecMeteo[ii].date = vecData[ii].date;
				vecMeteo[ii](param) = convertUnits(param, vecData[ii].val);
			}
		}

		//general case: merge one timestamp at a time
		std::vector<MeteoData> tmp;
		tmp.reserve( vecMeteo.size() + (vecData.size() - vecM_start)); //"worst case" scenario: all elements will be added

		size_t idx2 = vecM_start; //all previous elements were handled before
		size_t last_vM = vecM_start; //last element from vecMeteo that will have to be invalidated
		for(size_t ii=vecM_start; ii<vecMeteo.size(); ii++) {
			const Date curr_date( vecMeteo[ii].date );
			while ((idx2<vecData.size()) && (curr_date>vecData[idx2].date)) {
				tmp.push_back( md_pattern );
				tmp.back().date = vecData[idx2].date;
				tmp.back()(param) = convertUnits(param, vecData[idx2].val);
				idx2++;
			}
			if (idx2==vecData.size())  break; //nothing left to merge

			if (curr_date==vecData[idx2].date) {
				vecMeteo[ii](param) = convertUnits(param, vecData[idx2].val);
				idx2++;
			}
			tmp.push_back( vecMeteo[ii] );
			last_vM = ii;
		}

		const size_t new_count = last_vM - vecM_start + 1;
		if (new_count<tmp.size())
			vecMeteo.insert( vecMeteo.begin() + vecM_start, tmp.size()-new_count, tmp.front()); //so room for the extra params is allocated

		for(size_t ii=0; ii<tmp.size(); ii++)
			vecMeteo[vecM_start+ii] = tmp[ii];

		vecM_end = idx2;

		//filling data after vecMeteo
		if (vecMeteo.back().date<vecData.back().date) {
			if (vecM_end!=vecData.size()) {
				for (size_t ii=vecM_end; ii<vecData.size(); ii++) {
					vecMeteo.push_back( md_pattern );
					vecMeteo.back().date = vecData[ii].date;
					vecMeteo.back()(param) = convertUnits(param, vecData[ii].val);
				}
			}
		}
	}
}

size_t DBO::data_write(void* buf, const size_t size, const size_t nmemb, void* userp)
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
