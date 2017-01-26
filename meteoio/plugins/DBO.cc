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
 * - STATION#: station code for the given station, e. g. la_fouly_1034 (case sensitive!)
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
 */

const int DBO::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string DBO::sensors_endpoint = "sensors";
const std::string DBO::sensors_format = "format=csv";
const std::string DBO::null_string = "null";

DBO::DBO(const std::string& configfile)
      : cfg(configfile), vecStationName(),
        endpoint(), userid(), passwd(), default_timezone(1.),
        http_timeout(http_timeout_dflt), dbo_debug(false)
{
	initDBOConnection();
	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

DBO::DBO(const Config& cfgreader)
      : cfg(cfgreader), vecStationName(),
        endpoint(), userid(), passwd(), default_timezone(1.),
        http_timeout(http_timeout_dflt), dbo_debug(false)
{
	initDBOConnection();
	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

void DBO::initDBOConnection() {
	curl_global_init(CURL_GLOBAL_ALL);

	cfg.getValue("DBO_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Input", default_timezone, IOUtils::nothrow);

	cfg.getValue("DBO_URL", "Input", endpoint);
	if (*endpoint.rbegin() != '/') endpoint += "/";
	cerr << "[i] Using DBO URL: " << endpoint << endl;

	cfg.getValue("DBO_USER", "Input", userid, IOUtils::nothrow);
	cfg.getValue("DBO_PASS", "Input", passwd, IOUtils::nothrow);
	cfg.getValue("DBO_DEBUG", "INPUT", dbo_debug, IOUtils::nothrow);
}

void DBO::readStationData(const Date& /*date*/, std::vector<StationData>& vecStation)
{
	vecStation.clear();

}

void DBO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                          std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
}

//this method is called on each station in order to parse the header and set the metadata
bool DBO::parseMetadata(std::stringstream& /*ss*/, StationData &/*sd*/, std::string &/*fields*/, std::string &/*units*/) const
{

	return true;
}

void DBO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	const std::string station_id = vecStationName[stationindex];
	const string anon_request = sensors_endpoint + "/" + IOUtils::strToLower( station_id ) + "?" + sensors_format + "&from=" + dateStart.toString(Date::ISO)
	                            + "&to=" + dateEnd.toString(Date::ISO);
	const string auth = "&username=" + userid + "&password=" + passwd;
	const string request = (!userid.empty())? anon_request+auth : anon_request;

	stringstream ss;
	if (curl_read(request, ss)) {
		vector<size_t> index;
		string line, fields, units;

		MeteoData tmpmeteo;
		const bool meta_status = parseMetadata(ss, tmpmeteo.meta, fields, units); //read just one station

		if (units.empty() || fields.empty() || meta_status==false) {
			//when printing out a DBO error message, the # and ' ' have to be stripped from the begining -> substr(2)
			if (ss.str().find("doesn't exist in DBO!") != std::string::npos)
				throw NotFoundException(ss.str().substr(2), AT);
			if (ss.str().find("doesn't have access to the sensor") != std::string::npos)
				throw AccessException(ss.str().substr(2), AT);
			if (ss.str().find("There is no user with the provided") != std::string::npos)
				throw AccessException(ss.str().substr(2), AT);
			if (ss.str().find("The query consumed too many server resources!") != std::string::npos)
				throw IOException(ss.str().substr(2), AT);
			if (dbo_debug) {
				std::cout << "****\nRequest: " << request << "\n";
				std::cout << "Reply: " << ss.str() << "\n****\nPlease check the station name!\n";
			}
			throw InvalidFormatException("Invalid header for station " + station_id, AT);
		}
		//map_parameters(fields, units, tmpmeteo, index);

		do { //parse data section, the first line should already be buffered
			if (line.empty() || (line[0] == '#') || !isdigit(line[0])) continue; //skip empty lines
			//parse_streamElement(line, index, tmpmeteo);
			vecMeteo.push_back( tmpmeteo );
		} while (getline(ss, line));
		if (vecMeteo.empty()) //ie the data section was empty
			vecMeteo.push_back( tmpmeteo );
	} else {
		if (dbo_debug)
			std::cout << "****\nRequest: " << request << "\n****\n";
		throw IOException("Could not retrieve data for station " + station_id, AT);
	}
}

void DBO::convertUnits(MeteoData& meteo) const
{
	//converts C to Kelvin, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	ta = IOUtils::C_TO_K(ta);

	double& tsg = meteo(MeteoData::TSG);
	tsg = IOUtils::C_TO_K(tsg);

	double& tss = meteo(MeteoData::TSS);
	tss = IOUtils::C_TO_K(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
		rh /= 100.;

	double& hs = meteo(MeteoData::HS);
	if (hs != IOUtils::nodata)
		hs /= 100.;
}

size_t DBO::data_write(void* buf, size_t size, size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os = *static_cast<ostream*>(userp);
		const streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool DBO::curl_read(const std::string& url_query, std::ostream& os)
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	const string url = endpoint + url_query;

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
