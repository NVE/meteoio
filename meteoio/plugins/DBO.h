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

		typedef struct ts_Meta {
			ts_Meta(const Date& i_since, const Date& i_until, const std::string& i_agg_type, const double& i_ts_id, const unsigned int& i_interval)
			                : since(i_since), until(i_until), agg_type(i_agg_type), ts_id(i_ts_id), interval(i_interval) {}

			std::string toString() {
				std::ostringstream os;
				os << ts_id << " [";
				os << ((since.isUndef())? "-∞" : since.toString(Date::ISO)) << " - ";
				os << ((until.isUndef())? "∞" : until.toString(Date::ISO)) << "] ";
				os << agg_type << " - " << interval << " s";
				return os.str();
			}

			Date since, until;
			std::string agg_type;
			double ts_id;
			unsigned int interval;
		} tsMeta;
	private:
		void fillStationMeta();
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		void readTimeSerie(const unsigned int& ts_id, const MeteoData::Parameters& param, const Date& dateStart, const Date& dateEnd, const StationData& sd, std::vector<MeteoData>& vecMeteo);

		void initDBOConnection();
		static size_t data_write(void* buf, size_t size, size_t nmemb, void* userp);
		bool curl_read(const std::string& url, std::ostream& os);

		const Config cfg;
		std::vector<std::string> vecStationName;
		std::vector<StationData> vecMeta;
		std::vector< std::map<std::string, std::vector<DBO::tsMeta> > > vecTsMeta; ///< for every station, a map that contains for each parameter the relevant timeseries properties
		std::string coordin, coordinparam, coordout, coordoutparam; ///< projection parameters
		std::string endpoint; ///< Variables for endpoint configuration
		double default_timezone;
		int http_timeout; //time out for http connections
		bool dbo_debug;

		static const int http_timeout_dflt;
		static const std::string metadata_endpoint, data_endpoint, null_string;
};

} //end namespace mio

#endif
