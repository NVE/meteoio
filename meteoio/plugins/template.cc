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
#include "template.h"

using namespace std;

namespace mio {
/**
 * @page template TEMPLATE
 * @section template_format Format
 * *Put here the informations about the standard format that is implemented*
 *
 * @section template_units Units
 *
 *
 * @section template_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

const double TEMPLATE::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

TEMPLATE::TEMPLATE(void (*delObj)(void*), const std::string& filename) : IOInterface(delObj), cfg(filename)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

TEMPLATE::TEMPLATE(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

TEMPLATE::TEMPLATE(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

TEMPLATE::~TEMPLATE() throw()
{
	
}

void TEMPLATE::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
					 std::vector< std::vector<MeteoData> >& /*vecMeteo*/, 
					 std::vector< std::vector<StationData> >& /*vecStation*/,
					 const unsigned int&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                          const std::vector< std::vector<StationData> >& /*vecStation*/,
                          const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void TEMPLATE::cleanup() throw()
{
	
}

#ifndef _METEOIO_JNI
extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const string& classname, const string& filename) {
		if(classname == "TEMPLATE") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new TEMPLATE(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
