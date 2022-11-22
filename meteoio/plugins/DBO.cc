// SPDX-License-Identifier: LGPL-3.0-or-later
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

#include <sstream>
#include <regex>

#include <curl/curl.h>
#include <meteoio/thirdParty/picojson.h>

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
 * - DBO_URL: The URL of the RESTful web service (default: https://pgdata.int.slf.ch)
 * - STATION#: station code for the given station, prefixed by the network it belongs ot (for example: IMIS::SLF2, by default the network is assumed to be IMIS)
 * - DBO_TIMEOUT: timeout (in seconds) for the connection to the server (default: 60s)
 * - DBO_DEBUG: print the full requests/answers from the server when something does not work as expected (default: false)
 *
 * @code
 * METEO	= DBO
 * STATION1	= WFJ2
 * STATION2	= SMN::*WFJ1
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

static const double dbo_tz = 0.; //assuming GMT

std::string indent(const unsigned int& depth) 
{
	std::ostringstream ss;
	for (unsigned int jj=0; jj<depth; jj++) 
		ss << "\t";
	return ss.str();
}

//we keep these below as a simple functions in order to avoid exposing picojson stuff in the header
void printJSON(const picojson::value& v, const unsigned int& depth)
{
	if (v.is<picojson::null>()) {
		std::cout << indent(depth) << "NULL\n";
		return;
	}

	if (v.is<picojson::object>()) {
		//NOTE: v.get() will only be called once, see
		//https://stackoverflow.com/questions/15766020/does-a-c11-range-based-for-loop-condition-get-evaluated-every-cycle
		for (const auto& it : v.get<picojson::object>()) { //std::map<std::string,picojson::value>
			std::cout << indent(depth) << it.first << "\n";
			printJSON(it.second, depth+1);
		}
	} else if (v.is<std::string>()){
		std::cout << indent(depth) << v.get<std::string>() << "\n";
	} else if (v.is<double>()){
		std::cout << indent(depth) << v.get<double>() << "\n";
	} else if (v.is<bool>()){
		std::cout << indent(depth) << std::boolalpha << v.get<bool>() << "\n";
	} else if (v.is<picojson::array>()){ //ie vector<picojson::value>
		const auto& array = v.get<picojson::array>();
		std::cout << indent(depth) << "array " << array.size() << "\n";
		for (const auto& vec_elem : array)
			printJSON(vec_elem, depth+1);
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
		for (auto& keyvalue : v.get<picojson::object>()) { //std::map<std::string,picojson::value>
			if (keyvalue.first==local_path) {
				if (!remaining_path.empty())
					goToJSONPath(remaining_path, keyvalue.second);
				else
					return keyvalue.second;
			}
		}
	}

	return picojson::value();
}

void JSONQuery(const std::string& path, const picojson::value& v, std::vector<picojson::value>& results)
{
	if (v.is<picojson::null>()) return;

	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	if (v.is<picojson::object>()) {
		for (const auto& keyvalue : v.get<picojson::object>()) { //std::map<std::string,picojson::value>
			if (keyvalue.first==local_path) {
				if (!remaining_path.empty()) {
					 if (keyvalue.second.is<picojson::array>()){ //ie vector<picojson::value>
						for (const auto& vec_elem : keyvalue.second.get<picojson::array>())
							JSONQuery(remaining_path, vec_elem, results);
					} else
						JSONQuery(remaining_path, keyvalue.second, results);
				} else {
					results.push_back( keyvalue.second );
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
		if (! results.front().is<picojson::null>() && results.front().is<std::string>()) return  results.front().get<std::string>();
	}

	return std::string();
}

std::vector<std::string> getStrings(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);

	std::vector<std::string> vecString;
	for (const auto& result : results) {
		 if (result.is<picojson::array>()){
			for (const auto& vec_elem : result.get<picojson::array>()) { //loop over vector<picojson::value>
				if (! vec_elem.is<picojson::null>() && vec_elem.is<std::string>()) vecString.push_back( vec_elem.get<std::string>() );
			}
		} else
			if (! result.is<picojson::null>() && result.is<std::string>()) vecString.push_back( result.get<std::string>() );
	}

	return vecString;
}

double getDouble(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);
	if (!results.empty()) {
		if (! results.front().is<picojson::null>() && results.front().is<double>()) return  results.front().get<double>();
	}

	return IOUtils::nodata;
}

std::vector<double> getDoubles(const std::string& path, picojson::value& v)
{
	std::vector<picojson::value> results;
	JSONQuery(path, v, results);

	std::vector<double> vecDouble;
	for (const auto& result : results) {
		 if (result.is<picojson::array>()) {
			for (const auto& vec_elem : result.get<picojson::array>()) {//loop over vector<picojson::value>
				if (! vec_elem.is<picojson::null>() && vec_elem.is<double>()) vecDouble.push_back( vec_elem.get<double>() );
			}
		} else
			if (! result.is<picojson::null>() && result.is<double>()) vecDouble.push_back( result.get<double>() );
	}

	return vecDouble;
}

bool parseTsPoint(const picojson::value& v, Date& datum, double& value)
{
	if (!v.is<picojson::array>()) return false;

	const picojson::array& array = v.get<picojson::array>();
	if (array.size()!=2) return false;

	if (array[0].is<std::string>())
		IOUtils::convertString(datum, array[0].get<std::string>(), dbo_tz);
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

const std::vector<DBO::tsData> parseTimeSerie(const size_t& tsID, const double& factor, const double& offset, picojson::value& v)
{
	const picojson::value ts( goToJSONPath("$.measurements", v) );
	if (!ts.is<picojson::array>())
		throw InvalidFormatException("Could not parse timeseries "+IOUtils::toString(tsID), AT);

	const picojson::array& vecRaw = ts.get<picojson::array>();
	if (vecRaw.empty()) return std::vector<DBO::tsData>();

	std::vector<DBO::tsData> vecData( vecRaw.size() );
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		double value;
		Date datum;
		if (!parseTsPoint(vecRaw[ii], datum, value))  {
			printJSON(vecRaw[ii], 0);
			throw InvalidFormatException("Error parsing element "+IOUtils::toString(ii)+" of timeseries "+IOUtils::toString(tsID), AT);
		}

		if (value!=IOUtils::nodata) value = value * factor + offset;
		vecData[ii] = DBO::tsData(datum, value);
	}

	return vecData;
}

//get the properties of all timeseries belonging to the current station
std::vector<DBO::tsMeta> getTsProperties(picojson::value& v)
{
	std::vector<DBO::tsMeta> tsVec;
	std::vector<picojson::value> results;
	JSONQuery("$.timeseries", v, results);

	for (size_t ii=0; ii<results.size(); ii++) {
		if (!results[ii].is<picojson::array>()) continue;

		for (const auto& ts_obj : results[ii].get<picojson::array>()) { //loop over all provided timeseries
			if (ts_obj.is<picojson::null>()) continue;
			
			std::string code, device_code, agg_type;
			double id = -1., seqNr = 1;
			double interval=0, ts_offset=0;
			Date since, until;

			for (const auto& keyValue : ts_obj.get<picojson::object>()) { //loop over all key/values of a given timeseries
				const std::string key( keyValue.first );

				if (key=="sequenceNumber" && keyValue.second.is<double>()) seqNr = keyValue.second.get<double>();
				if (key=="id" && keyValue.second.is<double>()) id = keyValue.second.get<double>();
				if (key=="measurandCode" && keyValue.second.is<std::string>()) code = keyValue.second.get<std::string>();
				if (key=="deviceCode" && keyValue.second.is<std::string>()) device_code = keyValue.second.get<std::string>();
				if (key=="since" && keyValue.second.is<std::string>()) IOUtils::convertString(since, keyValue.second.get<std::string>(), 0.);
				if (key=="until" && keyValue.second.is<std::string>()) IOUtils::convertString(until, keyValue.second.get<std::string>(), 0.);
				if (key=="aggregationType" && keyValue.second.is<std::string>()) agg_type = keyValue.second.get<std::string>();
				if (key=="measureIntervalInMinutes" && keyValue.second.is<double>()) interval = keyValue.second.get<double>(); //TODO check that it can be cast
				if (key=="measureIntervalOffsetInMinutes" && keyValue.second.is<double>()) ts_offset = keyValue.second.get<double>();
			}

			//reject some timeseries
			if (seqNr!=1) continue; //HACK per WIS, only consider seq number 1
			if (device_code=="BATTERY" || device_code=="LOGGER") continue;
			if (device_code=="MODEL_SNOWPACK") continue; //TODO add an option to use Snowpack computed parameters
			if (agg_type=="SD") continue; //we don't care about standard deviation anyway
			if (id==-1.) continue; //no id was provided

			const std::string param_dbo( IOUtils::strToUpper( code.substr(0, code.find('_')) ) );
			const std::string parname( DBO::getParameter(param_dbo, agg_type) );
			if (parname.empty()) continue;
			
			DBO::tsMeta tmp(param_dbo, since, until, agg_type, static_cast<size_t>(id), static_cast<unsigned int>(interval*60.), static_cast<unsigned int>(ts_offset*60.), static_cast<unsigned int>(seqNr));
			tmp.parname = parname;
			DBO::setUnitsConversion(tmp);
			tsVec.push_back( tmp );
		}
	}

	return tsVec;
}

/*************************************************************************************************/
//example metadata query: https://pgdata.int.slf.ch/data/stations/IMIS/WFJ2
const int DBO::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string DBO::endpoint_default = "https://pgdata.int.slf.ch";
const std::string DBO::metadata_api = "/data/stations/";
const std::string DBO::data_api = "/data/timeseries/";

DBO::DBO(const std::string& configfile)
      : vecStationName(), vecMeta(), vecTsMeta(),
        coordin(), coordinparam(),
        endpoint(endpoint_default),
        http_timeout(http_timeout_dflt), dbo_debug(false)
{
	const Config cfg( configfile );
	initDBOConnection(cfg);
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

DBO::DBO(const Config& cfgreader)
      : vecStationName(), vecMeta(), vecTsMeta(),
        coordin(), coordinparam(),
        endpoint(endpoint_default),
        http_timeout(http_timeout_dflt), dbo_debug(false)
{
	initDBOConnection(cfgreader);
	IOUtils::getProjectionParameters(cfgreader, coordin, coordinparam);
	cfgreader.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

void DBO::initDBOConnection(const Config& cfg)
{
	curl_global_init(CURL_GLOBAL_ALL);

	cfg.getValue("DBO_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	cfg.getValue("DBO_URL", "Input", endpoint, IOUtils::nothrow);
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

/**
* @brief Read and cache the stations' metadata
*/
void DBO::fillStationMeta()
{
	static const std::regex stat_id_regex("([^:]+)::([^:]+)");
	std::smatch stat_id_matches;

	vecMeta.clear();
	vecMeta.resize( vecStationName.size() );
	vecTsMeta.resize( vecStationName.size() );

	for(size_t ii=0; ii<vecStationName.size(); ii++) {
		std::string station_id( IOUtils::strToUpper(vecStationName[ii]) );
		std::string network = "IMIS";

		if (std::regex_match(station_id, stat_id_matches, stat_id_regex)) {
			network = stat_id_matches.str(1);
			station_id = stat_id_matches.str(2);
		}

		const std::string request( metadata_api + network + "/" + station_id );

		std::stringstream ss;
		if (curl_read(request, ss)) {
			picojson::value v;
			
			//handling possible errors
			const std::string err( picojson::parse(v, ss.str()) );
			if (!err.empty()) throw IOException("Error while parsing JSON: "+err, AT);
			std::string error_msg = getString("$.error", v);
			if (!error_msg.empty()) throw UnknownValueException("Error with station '"+network+"::"+station_id+"': "+error_msg, AT);

			Coords position(coordin, coordinparam);
			position.setLatLon(getDouble("$.location.lat", v), getDouble("$.location.lon", v), getDouble("$.location.elevation", v));
			const StationData sd(position, getString("$.code", v), getString("$.label", v));
			vecMeta[ii] = sd;

			//parse and store the time series belonging to this station
			vecTsMeta[ii] = getTsProperties(v);
		} else {
			if (dbo_debug)
				std::cout << "****\nRequest: " << request << "\n****\n";
			throw IOException("Could not retrieve data for station " + station_id, AT);
		}
		
		if (dbo_debug) {
			std::cout << "<Station " << station_id << ">\n";
			for (const auto& ts : vecTsMeta[ii]) std::cout << ts.toString() << "\n";
			std::cout << "</Station " << station_id << ">\n";
		}
	}
}

/**
* @brief Identify the relevant MeteoData::Parameters from DBO provided information
* @param[in] param_str DBO string representation of the meteo parameter
* @param[in] agg_type DBO aggregation type
* @return standardized parameter name or empty string
*/
std::string DBO::getParameter(const std::string& param_str, const std::string& agg_type)
{
	if (param_str=="P") return "P";
	else if (param_str=="TA") return "TA";
	else if (param_str=="RH") return "RH";
	else if (param_str=="TS0") return "TSG";
	else if (param_str=="TSS") return "TSS";
	else if (param_str=="HS") return "HS";
	else if (param_str=="VW" && agg_type=="MAX") return "VW_MAX";
	else if (param_str=="VW") return "VW";
	else if (param_str=="DW") return "DW";
	else if (param_str=="RSWR") return "RSWR";
	else if (param_str=="ISWR") return "ISWR";
	else if (param_str=="ILWR") return "ILWR";
	else if (param_str=="RRI") return "PSUM";
	else if (param_str=="TS25") return "TS1";
	else if (param_str=="TS50") return "TS2";
	else if (param_str=="TS100") return "TS3";
	else if (param_str=="TG10") return "TSOIL10";
	else if (param_str=="TG30") return "TSOIL30";
	else if (param_str=="TG50") return "TSOIL50";
	else return "";
}

/**
* @brief Provide the way to convert the DBO units into standardized units (SI).
* It is assume that we can first multiply by a factor, then add an offset.
* @param[in] ts DBO timeseries properties
*/
void DBO::setUnitsConversion(DBO::tsMeta& ts)
{
	//compute the conversion parameters (C to K, cm to m, % to [0-1], PINT to PSUM
	if (ts.parname=="TA" || ts.parname=="TSG" || ts.parname=="TSS") {
		ts.units_offset = Cst::t_water_freezing_pt;
	} else if(ts.parname=="RH" || ts.parname=="HS") {
		ts.units_factor = 0.01;
	} else if(ts.parname=="PSUM") {
		ts.units_factor = 3600. / ts.interval;
	} else if(ts.parname=="TS1" || ts.parname=="TS2" || ts.parname=="TS3") {
		ts.units_offset = Cst::t_water_freezing_pt;
	} else if(ts.parname=="TSOIL10" || ts.parname=="TSOIL30" || ts.parname=="TSOIL50") {
		ts.units_offset = Cst::t_water_freezing_pt;
	}
	
	return;
}

//read all data for the given station
void DBO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	const std::string Start( dateStart.toString(Date::ISO_Z) );
	const std::string End( dateEnd.toString(Date::ISO_Z) );
	const StationData &sd = vecMeta[stationindex];

	//now get the data
	for (const DBO::tsMeta &ts : vecTsMeta[stationindex]) {
		
		//check if the current ts contains our period of interest
		const Date tsStart(ts.since), tsEnd(ts.until);
		if ((!tsStart.isUndef() && tsStart>dateEnd) || (!tsEnd.isUndef() && tsEnd<dateStart)) continue;

		std::ostringstream ss_ID; ss_ID << ts.id;
		const std::string request( data_api + ss_ID.str() + "?from=" + Start + "&until=" + End );

		std::stringstream ss;
		if (curl_read(request, ss)) { //retrieve the page from the formed URL
			if (ss.str().empty()) throw UnknownValueException("Timeseries not found: '"+ss_ID.str()+"'", AT);
			picojson::value v;
			const std::string err( picojson::parse(v, ss.str()) );
			if (!err.empty()) throw IOException("Error while parsing JSON: "+ss.str(), AT);

			const std::vector<DBO::tsData> vecData( parseTimeSerie(ts.id, ts.units_factor, ts.units_offset, v) );
			if (vecData.empty()) {
				if (dbo_debug) 
					std::cout << vecMeta[stationindex].getStationID() << " has no data for " << ts.parname << " in the requested period\n";
				continue;
			}

			MeteoData md_pattern = (vecMeteo.empty())? MeteoData(Date(), sd) : vecMeteo.front(); //This assumes that the station is not moving!
			if (!md_pattern.param_exists(ts.parname)) {
				md_pattern.addParameter( ts.parname );
				for (size_t jj=0; jj<vecMeteo.size(); jj++) vecMeteo[jj].addParameter( ts.parname ); //TODO rewrite addParameter to create and attribute a value
			}
			
			const size_t parindex = md_pattern.getParameterIndex( ts.parname );
			mergeTimeSeries(md_pattern, parindex, vecData, vecMeteo);
		} else {
			if (dbo_debug)
				std::cout << "****\nRequest: " << request << "\n****\n";
			throw IOException("Could not retrieve data for timeseries " + ss_ID.str(), AT);
		}
	}
}

/**
* @brief Merge a newly read timeseries into vecMeteo
* @param[in] md_pattern pattern MeteoData to be used to insert new elements
* @param[in] param index of the current meteo parameter
* @param[in] vecData the raw (but parsed) data
* @param vecMeteo the vector that will receive the new values (as well as inserts if necessary)
*/
void DBO::mergeTimeSeries(const MeteoData& md_pattern, const size_t& param, const std::vector<DBO::tsData>& vecData, std::vector<MeteoData>& vecMeteo) const
{
	if (vecData.empty()) return;

	if (vecMeteo.empty()) { //easy case: the initial vector is empty
		vecMeteo.resize( vecData.size() );
		for (size_t ii=0; ii<vecData.size(); ii++) {
			MeteoData md( md_pattern );
			md.date = vecData[ii].date;
			md(param) = vecData[ii].val;
			vecMeteo[ii] = md;
		}
	} else {
		size_t vecM_start = 0; //the index in vecRaw that matches the original start of vecMeteo
		size_t vecM_end = 0; //the index in vecRaw that matches the original end of vecMeteo

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
				vecMeteo[ii](param) = vecData[ii].val;
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
				tmp.back()(param) = vecData[idx2].val;
				idx2++;
			}
			if (idx2==vecData.size())  break; //nothing left to merge

			if (curr_date==vecData[idx2].date) {
				vecMeteo[ii](param) = vecData[idx2].val;
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
					vecMeteo.back()(param) = vecData[ii].val;
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

bool DBO::curl_read(const std::string& url_query, std::ostream& os) const
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
