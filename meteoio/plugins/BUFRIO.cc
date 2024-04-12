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
#include "BUFRIO.h"

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

const double BUFRIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

BUFRIO::BUFRIO(const std::string& configfile) : cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	
	/* Example: how to read keys from the Config object*/
	/*const double factor = cfg.get("PLUGIN_FACTOR", "Input"); //if the key PLUGIN_FACTOR is not found in the [Input] section, an exception will be thrown
	 * 
	 * bool enable_feature = false;
	 * cfg.getValue("ENABLE_FEATURE", "Input", enable_feature, IOUtils::nothrow); //if the key is not found, it simply keeps its previous value
	 * 
	 * int parameter = 0;
	 * cfg.getValue("PLUGIN_NR_PARAMS", "Output", parameter); //if the key is not found, an exception will be thrown
	 * 
	 * //it is also possible to get all the keys starting with a given pattern at once and then loop through them:
	 * std::vector<std::string> vecFilenames;
	* cfg.getValues("STATION", "INPUT", vecFilenames);
	 */
}

BUFRIO::BUFRIO(const Config& cfgreader) : cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}


void BUFRIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
	
}

void BUFRIO::readStationData(const Date &/* date */, std::vector<StationData> &vecStation) {
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}


} //namespace
