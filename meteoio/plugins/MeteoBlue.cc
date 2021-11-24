// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 SLF                                                                                                                                */
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
#include <meteoio/plugins/MeteoBlue.h>

#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/meteoLaws/Meteoconst.h>

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cstring>

#include <curl/curl.h>
#include <meteoio/thirdParty/picojson.h>

using namespace std;

namespace mio {
	
//we keep these below as a simple functions in order to avoid exposing picojson stuff in the header
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

picojson::value goToJSONPath(const std::string& path, const picojson::value& v)
{
	size_t start_pos = 0;
	if (path[0]=='$') start_pos++;
	if (path[1]=='.') start_pos++;

	const size_t end_pos = path.find(".", start_pos);
	const std::string local_path = (end_pos!=std::string::npos)? path.substr(start_pos, end_pos-start_pos) : path.substr(start_pos);
	const std::string remaining_path = (end_pos!=std::string::npos)? path.substr(end_pos+1) : "";

	if (v.is<picojson::object>()) {
		const picojson::value::object& obj = v.get<picojson::object>();
		for (std::map<std::string,picojson::value>::const_iterator it = obj.begin(); it != obj.end(); ++it) {
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

/*************************************************************************************************/
const int MeteoBlue::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string MeteoBlue::dflt_endpoint = "http://my.meteoblue.com/";
std::map< std::string, MeteoGrids::Parameters > MeteoBlue::params_map;
const bool MeteoBlue::__init = MeteoBlue::initStaticData();

bool MeteoBlue::initStaticData()
{
	params_map[ "precipitation" ] = MeteoGrids::PSUM;
	params_map[ "snowfraction" ] = MeteoGrids::PSUM_PH;	//WARNING the opposite of us!!
	params_map[ "temperature" ] = MeteoGrids::TA;
	params_map[ "relativehumidity" ] = MeteoGrids::RH;
	params_map[ "windspeed" ] = MeteoGrids::VW;
	params_map[ "winddirection" ] = MeteoGrids::DW;
	params_map[ "sealevelpressure" ] = MeteoGrids::P_SEA;
	params_map[ "skintemperature" ] = MeteoGrids::TSS;
	params_map[ "ghi_total" ] = MeteoGrids::ISWR;
	params_map[ "dif_total" ] = MeteoGrids::ISWR_DIFF;
	params_map[ "surfaceairpressure" ] = MeteoGrids::P;
	params_map[ "gust" ] = MeteoGrids::VW_MAX;
	
	return true;
}

MeteoBlue::MeteoBlue(const std::string& configfile)
      : cfg(configfile), vecMeta(),
        coordin(), coordinparam(),
        endpoint(dflt_endpoint), apikey(), packages(),
        http_timeout(http_timeout_dflt), debug(false)
{
	
	init();
}

MeteoBlue::MeteoBlue(const Config& cfgreader)
      : cfg(cfgreader), vecMeta(),
        coordin(), coordinparam(),
        endpoint(dflt_endpoint), apikey(), packages(),
        http_timeout(http_timeout_dflt), debug(false)
{
	init();
}

void MeteoBlue::init()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);

	cfg.getValue("METEOBLUE_DEBUG", "INPUT", debug, IOUtils::nothrow);
	cfg.getValue("METEOBLUE_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	cfg.getValue("METEOBLUE_APIKEY", "Input", apikey);
	cfg.getValue("METEOBLUE_URL", "Input", endpoint, IOUtils::nothrow);
	if (*endpoint.rbegin() != '/') endpoint += "/";
	
	//building the packages arguments
	std::vector<std::string> vecStr;
	cfg.getValue("METEOBLUE_PACKAGES", "Input", vecStr);
	packages = vecStr[0];
	for (size_t ii=1; ii<vecStr.size(); ii++) packages.append( "_"+vecStr[ii] );
	
	//reading the stations' coordinates to retrieve
	const std::vector< std::pair<std::string, std::string> > vecStationSpecs( cfg.getValues("STATION", "Input") );
	for (size_t ii=0; ii<vecStationSpecs.size(); ii++) {
		if (vecStationSpecs[ii].first.find('_') != std::string::npos) continue; //so we skip the other station#_xxx parameters
		
		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		const Coords curr_point(coordin, coordinparam, vecStationSpecs[ii].second);
		if (curr_point.isNodata()) continue;
		
		const std::string id_num( vecStationSpecs[ii].first.substr(std::string("STATION").length()) );
		const bool has_id = cfg.keyExists("STATION"+id_num+"_ID", "Input");
		const std::string stat_id = (has_id)? cfg.get("STATION"+id_num+"_ID", "Input"): "STAT"+id_num;
		const bool has_name = cfg.keyExists("STATION"+id_num+"_NAME", "Input");
		const std::string stat_name = (has_name)? cfg.get("STATION"+id_num+"_NAME", "Input"): "STATION"+id_num;
		
		vecMeta.push_back( StationData(curr_point, stat_id, stat_name) );
	}
	if (vecMeta.empty())
		throw InvalidArgumentException("Please provide stations' coordinates for the METEOBLUE plugin!", AT);
	
	curl_global_init(CURL_GLOBAL_ALL);
}

void MeteoBlue::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation = vecMeta;
}

void MeteoBlue::readMeteoData(const Date& dateStart, const Date& dateEnd,
                          std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	
	vecMeteo.resize(vecMeta.size());
	for(size_t ii=0; ii<vecMeta.size(); ii++)
		readData(dateStart, dateEnd, vecMeta[ii], vecMeteo[ii]);
}

void MeteoBlue::readTime(const picojson::value &v, const StationData& sd, std::vector<MeteoData> &vecMeteo) const
{
	picojson::value time( goToJSONPath("$.time", v) );
	if (!time.is<picojson::array>()) {
		throw InvalidFormatException("Could not parse the data section '"+time.to_str()+"'", AT);
	}
	const picojson::array vecRaw = time.get<picojson::array>();
	
	vecMeteo.reserve( vecRaw.size() );
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		if (!vecRaw[ii].is<std::string>()) {
			std::ostringstream ss; ss << "Could not parse date '" << vecRaw[ii].to_str() << "'";
			throw InvalidFormatException(ss.str(), AT);
		}
		Date datum;
		IOUtils::convertString(datum, vecRaw[ii].get<std::string>(), 0.); //hard-coding GMT as this is what is delivered in ISO8601
		vecMeteo.push_back( MeteoData(datum, sd) );
	}
}

void MeteoBlue::readParameter(const picojson::value &v, const std::string& paramname, std::vector<MeteoData> &vecMeteo) const
{
	const std::map< std::string, MeteoGrids::Parameters >::const_iterator it = params_map.find( paramname );
	if (it==params_map.end()) return;
	const std::string mio_parname( MeteoGrids::getParameterName( it->second ) );
	
	picojson::value param( goToJSONPath("$."+paramname, v) );
	if (!param.is<picojson::array>()) {
		throw InvalidFormatException("Could not parse the data section '"+param.to_str()+"'", AT);
	}
	const picojson::array vecRaw = param.get<picojson::array>();
	if (vecMeteo.size()!=vecRaw.size()) 
		throw InvalidFormatException("Trying to insert "+IOUtils::toString(vecRaw.size())+" values into vecMeteo of size "+IOUtils::toString(vecMeteo.size()), AT);
	
	for (size_t ii=0; ii<vecRaw.size(); ii++) {
		if (vecRaw[ii].is<picojson::null>()) continue;
		if (!vecRaw[ii].is<double>())
			throw InvalidFormatException("Could not parse '"+vecRaw[ii].to_str()+"' as double", AT);
		
		const size_t parindex = vecMeteo[ii].addParameter( mio_parname );
		vecMeteo[ii]( parindex ) = vecRaw[ii].get<double>();
	}
}

void MeteoBlue::readData(const Date& dateStart, const Date& dateEnd, const StationData& sd, std::vector<MeteoData> &vecMeteo)
{
	//build the query URL
	std::ostringstream url;
	/*url << endpoint << "packages/" << packages << "?lat=" << sd.position.getLat() << "&lon=" << sd.position.getLon();
	const double altitude( sd.position.getAltitude() );
	if (altitude!=IOUtils::nodata) url << "&asl=" << altitude;
	url << "&timeformat=iso8601" << "&history_days=3" << "&apikey=" << apikey;*/
	url << "http://my.meteoblue.com/packages/basic-1h_agro-1h?lat=47.56&lon=7.57&apikey=DEMOKEY&sig=e85c990f1d5d476b29eddd989ca56859";
	
std::cout << "Using the following URL: " << url.str() << "\n";
	
	std::stringstream ss;
	if (curl_read(url.str(), ss)) { //retrieve the page from the formed URL
		if (ss.str().empty()) throw UnknownValueException("No data returned for query '"+url.str()+"'", AT);
		picojson::value v;
		const std::string err( picojson::parse(v, ss.str()) );
		if (!err.empty()) throw IOException("Error while parsing JSON: "+ss.str(), AT);
		
		//go to dataset, loop over the parameters
		picojson::value data1h( goToJSONPath("$.data_1h", v) );
		readTime(data1h, sd, vecMeteo);
		
		const picojson::value::object& obj = data1h.get<picojson::object>();
		for (picojson::value::object::const_iterator it = obj.begin(); it != obj.end(); ++it) {
			std::cout << it->first << ": " << it->second.to_str() << std::endl;
			const std::string paramname( it->first );
			
			if (paramname!="time") readParameter(data1h, paramname, vecMeteo);
		}
	}
	
	std::cout << "Read:\n";
	for (size_t ii=0; ii<vecMeteo.size(); ii++) std::cout << vecMeteo[ii].toString();
	std::cout << "\n";
	
}

size_t MeteoBlue::data_write(void* buf, const size_t size, const size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os = *static_cast<ostream*>(userp);
		const std::streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool MeteoBlue::curl_read(const std::string& url, std::ostream& os) const
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	if (curl) {
		if (CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &data_write))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FILE, &os))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, http_timeout))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_URL, url.c_str())))
		{
			code = curl_easy_perform(curl);
		}
		curl_easy_cleanup(curl);
	}

	if (code!=CURLE_OK) {
		if (debug)
			std::cout << "****\nRequest: " << url << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}

	return (code==CURLE_OK);
}

} //namespace
