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
 * BUFR is a binary data format, that was introduced by the WMO to universally store meteorological observation data. It is very flexible and powerful, 
 * but can be quite complex to handle. An introduction can be found at 	https://en.wikipedia.org/wiki/BUFR. 
 * 
 * A BUFR file can contain multiple messages, each of which can contain multiple subsets. Where the data is stored with its metadata and can be read by 
 * its respective key. Here we use keys from the WMO element table 40 (https://confluence.ecmwf.int/display/ECC/WMO%3D40+element+table) to read the data.
 * 
 * To inspect your BUFR files we recommend using bufr tools: (https://confluence.ecmwf.int/display/ECC/BUFR+tools) or the 
 * metview software (https://metview.readthedocs.io/en/latest/). These are provided by ECMWF and therefore use the same BUFR reading library (ecCodes) as MeteoIO.
 * 
 * Prior to reading the data a scan of the file is done, to find all stations in the file, together with their meta data, and if it is valid (we require position (lat,lon,alt) and station ID) 
 * If no station ID but a station Name is found, we use this as station ID. Additionally we require unique timepoints per station.
 * 
 * So far only reading parameters known to MeteoIO is supported. And the BUFR file can not contain multiple datapoints for the same station and time. Such as
 * for height profiles. I.e., we require only static station observations.
 * @note If you require these functionality please reach out to us, so we can improve on this plugin.
 * 
 * It is not yet possible to write BUFR files with MeteoIO, but will soon be possible.
 * 
 * @section bufrio_keywords Keywords
 * This plugin uses the following keywords:
 * - METEOPATH: Path to the BUFR files
 * - BUFREXT: Extension of the BUFR files (default: .bufr)
 * - STATION#: File name per station
 * - ADDITIONAL_PARAMS: Additional parameters to read from the BUFR file, as a comma separated list of WMO element table 40 keys
 * - FALLBACKTZ: Fallback timezone in hours, if no timezone information is found in the BUFR file (default: GMT+0)
 * - VERBOSE: Print additional information during reading (default: false)
 * 
 * @section bufrio_example Example
 * @code
 * METEO = BUFR
 * METEOPATH = path/to/bufr/files
 * BUFREXT = .bfr
 * ADDITIONAL_PARAMS = airTemperature, dewPointTemperature
 * STATION1 = example.bfr
 * VERBOSE = TRUE
 * @endcode
 * 
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
