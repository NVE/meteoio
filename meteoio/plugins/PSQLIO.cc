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
#include "PSQLIO.h"

using namespace std;
using namespace pqxx;

namespace mio {
/**
 * @page template PSQLIO
 * @section psql_format Format
 * 
 *
 * @section psql_units Units
 *
 *
 * @section psql_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

const double PSQLIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

PSQLIO::PSQLIO(const std::string& configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PSQLIO::PSQLIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PSQLIO::~PSQLIO() throw()
{

}

void PSQLIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*name_in*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::read2DGrid(Grid2DObject& /*grid_out*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                             const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                              const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::write2DGrid(const Grid2DObject& /*grid_in*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PSQLIO::cleanup() throw()
{

}

} //namespace
