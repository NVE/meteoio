// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef WWCSIO_H
#define WWCSIO_H

#include <meteoio/IOInterface.h>

#include <string>

namespace mio {

/**
 * @class WWCSIO_H
 * @brief This is the plugin required to get meteorological data from the WWCS database.
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2022-02-24
 */
class WWCSIO : public IOInterface {
	public:
		WWCSIO(const std::string& configfile);
		WWCSIO(const WWCSIO&);
		WWCSIO(const Config& cfgreader);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		typedef struct DB_FIELD {
			DB_FIELD(const std::string& i_param) 
			        : param( i_param ), processing( 0 ), isDate(param=="DATETIME") {}
			DB_FIELD(const std::string& i_param, const unsigned int &i_processing) 
			        : param( i_param ), processing( i_processing ), isDate(false) {}
			const std::string param;
			const unsigned int processing;
			const bool isDate;
		} db_field;
		
		void readConfig();
		std::vector<std::string> readStationIDs() const;
		void readStationMetaData();
		void readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
		              const size_t& stationindex) const;

		const Config cfg;
		std::vector<std::string> vecStationIDs;
		std::vector<StationData> vecStationMetaData;
		std::string mysqlhost, mysqldb, mysqluser, mysqlpass;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		double in_dflt_TZ, out_dflt_TZ;
		static const size_t nrMeteoFields;
		unsigned int mysql_options;

		static const std::string MySQLQueryStationMetaData;
		static const std::string MySQLQueryMeteoData;
		static const std::vector< db_field > meteoFields;
};

} //namespace
#endif

