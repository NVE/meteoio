// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
// #include <meteoio/plugins/bufrio.h>
#include <meteoio/plugins/BUFRIO.h>
#include <meteoio/plugins/plugin_utils.h>
#include <unordered_set>
using namespace std;

namespace mio {
/**
 * @page bufrio BUFRIO
 * @section bufrio_format Format
 * *Put here the information about the standard format that is implemented*
 *
 * @section bufrio_units Units
 *
 *
 * @section bufrio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */
using namespace PLUGIN;
const double BUFRIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const std::string dflt_extension_BUFR = ".bufr";
const std::string BUFRIO::template_filename = "MeteoIO.bufr";

BUFRIO::BUFRIO(const std::string& configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), station_files(), additional_params()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	
	parseInputSection();
}

BUFRIO::BUFRIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), station_files(), additional_params()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputSection();
}

void BUFRIO::parseInputSection() {
	const std::string in_meteo = IOUtils::strToUpper(cfg.get("METEO", "Input", ""));
	if (in_meteo == "BUFR") { // keep it synchronized with IOHandler.cc for plugin mapping!!
        const std::string inpath = cfg.get("METEOPATH", "Input");
        std::vector<std::string> vecFilenames;
        cfg.getValues("STATION", "INPUT", vecFilenames);

		std::string bufr_ext = dflt_extension_BUFR;
		cfg.getValue("BUFREXT", "Input", bufr_ext, IOUtils::nothrow);
		if (bufr_ext == "none")
			bufr_ext.clear();

		cfg.getValue("ADDITIONAL_PARAMS", "INPUT", additional_params, IOUtils::nothrow);

		double fallback_tz = 0;
		cfg.getValue("FALLBACKTZ", "Input", fallback_tz, IOUtils::nothrow);

		bool verbose = false;
		cfg.getValue("VERBOSE", "Input", verbose, IOUtils::nothrow);

        if (vecFilenames.empty())
            scanMeteoPath(cfg, inpath, vecFilenames, bufr_ext);

        const std::vector<std::string> all_files_and_paths = getFilesWithPaths(vecFilenames, inpath, dflt_extension_BUFR);
		for (const auto &filename : all_files_and_paths) {
			station_files.push_back(BUFRFile(filename, coordin, verbose, fallback_tz));
		}
	}
}


void BUFRIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData>>& vecvecMeteo)
{
	vecvecMeteo.clear();
	std::map<std::string, size_t> station_ids;
	for (auto &station_file : station_files) {
		station_file.readData(vecvecMeteo, station_ids, additional_params);
	}
	// sort the data, as multiple bufr files can have information about the same station, but order is not guaranteed
	for (auto &vecMeteo : vecvecMeteo) {
		std::sort(vecMeteo.begin(), vecMeteo.end());
	}
}

void BUFRIO::readStationData(const Date &/* date */, std::vector<StationData> &vecStation) {
	vecStation.clear();
	std::unordered_set<std::string> station_id_set;
	for (const auto &station_file : station_files) {
		for (const auto &station : station_file.getMetadata()) {
			auto success = station_id_set.insert(station.second.getStationID());
			if (success.second) {
				vecStation.push_back(station.second);
			}
		}
	}
}

void BUFRIO::writeMeteoData(const std::vector<std::vector<MeteoData>> &/* vecMeteo */, const std::string &/* name */) {
	throw IOException("Writing BUFR files is not implemented", AT);
}


} //namespace
