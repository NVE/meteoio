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
#ifndef DBO_H
#define DBO_H

#include <meteoio/IOInterface.h>

#include <string>
#include <vector>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

/**
 * @class DBO
 * @brief This class enables the access to the DBO RESTful web service
 *
 * @ingroup plugins
 * @date   2017-01-26
 */

class DBO : public IOInterface {
	public:
		DBO(const std::string& configfile);
		DBO(const DBO&);
		DBO(const Config&);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void fillStationMeta();
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		void convertUnits(MeteoData& meteo) const;

		void initDBOConnection();
		static size_t data_write(void* buf, size_t size, size_t nmemb, void* userp);
		bool curl_read(const std::string& url, std::ostream& os);

		const Config cfg;
		std::vector<std::string> vecStationName;
		std::vector<StationData> vecMeta;
		std::string coordin, coordinparam, coordout, coordoutparam; ///< projection parameters
		std::string endpoint, userid, passwd; ///< Variables for endpoint configuration
		double default_timezone;
		int http_timeout; //time out for http connections
		bool dbo_debug;

		static const int http_timeout_dflt;
		static const std::string sensors_endpoint, null_string;
};

} //end namespace mio

#endif
