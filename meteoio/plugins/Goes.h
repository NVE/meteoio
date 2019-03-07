/***********************************************************************************/
/*  Copyright 2019 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef GOES_H
#define GOES_H

#include <meteoio/IOInterface.h>

namespace mio {

/**
 * @class GoesIO
 * @brief This plugin deals with data that has been transmitted through the
 * <A HREF="https://en.wikipedia.org/wiki/Geostationary_Operational_Environmental_Satellite">GOES</A> satellites
 * (see also https://www.goes-r.gov/resources/docs.html and https://www.rtl-sdr.com/tag/goes/). In order to reduce data
 * transfers, no headers are provided with the data and therefore all the metadata will have to be provided to the plugin.
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2019-03-05
 */
class GoesIO : public IOInterface {
	public:
		GoesIO(const std::string& configfile);
		GoesIO(const GoesIO&);
		GoesIO(const Config& cfgreader);
		
		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);

		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void parseInputOutputSection(const Config& cfgreader);
		void readRawGoes(const std::string& file_and_path, const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo);
		Date parseDate(const std::vector<float>& raw_data) const;
		MeteoData parseDataLine(const std::string& goesID, const Date& dt, const std::vector<float>& raw_data) const;
		void parseFieldsSpecs(const std::vector<std::string>& fieldsNames, MeteoData &meteo_template, std::vector<size_t> &idx);
		bool addStation(const std::string& goesID);

		std::vector<std::string> vecFilenames;
		std::map<std::string, std::vector<size_t> > fields_idx; ///< for each stationID, vector of MeteoData field index
		std::map<std::string, std::vector<double> > units_offset, units_multiplier;
		std::map<std::string, MeteoData> md_template;
		std::map<std::string, size_t> stationsIdx; ///< the station index within vecMeteo
		Config metaCfg;
		std::string meteopath;
		std::string coordin, coordinparam; //projection parameters
		double in_TZ, in_nodata;
		size_t stationID_idx, year_idx, jdn_idx, hour_idx;
		bool debug;
};

} //namespace
#endif
