/***********************************************************************************/
/*  Copyright 2016 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef OSHDIO_H
#define OSHDIO_H

#include <meteoio/IOInterface.h>

#include <string>

namespace mio {

/**
 * @class OshdIO
 * @brief This plugin reads Matlab binary files, relying on the MatIO library
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2016-03-24
 */
class OshdIO : public IOInterface {
	public:
		OshdIO(const std::string& configfile);
		OshdIO(const OshdIO&);
		OshdIO(const Config& cfgreader);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void parseInputOutputSection();
		void readSWRad(const Date& station_date, const std::string& file_suffix, const size_t& nrIDs, std::vector< std::vector<MeteoData> >& vecMeteo) const;
		void readFromFile(const std::string& filename, const MeteoData::Parameters& param, const Date& in_timestep, std::vector<double> &vecData) const;
		void buildVecIdx(const std::vector<std::string>& vecAcro);
		void fillStationMeta();
		
		size_t getFileIdx(const Date& start_date) const;
		static void scanMeteoPath(const std::string& meteopath_in,  std::vector< std::pair<mio::Date,std::string> > &meteo_files);
		static void checkFieldType(const MeteoData::Parameters& param, const std::string& type);
		static double convertUnits(const double& val, const std::string& units, const MeteoData::Parameters& param);
		
		const Config cfg;
		std::vector< std::pair<Date,std::string> > cache_meteo_files; //cache of meteo files in METEOPATH
		std::vector<StationData> vecMeta;
		std::vector<std::string> vecIDs; ///< IDs of the stations that have to be read
		std::vector< std::pair<MeteoData::Parameters, std::string> > params_map; ///< parameters to extract from the files
		std::vector<size_t> vecIdx; ///< index of each ID that should be read within the 'acro', 'names' and 'data' vectors
		std::string in_meteopath, in_metafile;
		bool debug; ///< write out extra information to help understand what is being read
		//std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		
		static const char* meteo_ext; //for the file naming scheme
		static const double in_dflt_TZ;     //default time zones
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
};

} //namespace
#endif
