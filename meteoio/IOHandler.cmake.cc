/***********************************************************************************/
/*  Copyright 2009-2012 WSL Institute for Snow and Avalanche Research  SLF-DAVOS   */
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

#include <meteoio/IOHandler.h>

#cmakedefine PLUGIN_ARCIO
#cmakedefine PLUGIN_A3DIO
#cmakedefine PLUGIN_ARPSIO
#cmakedefine PLUGIN_GRASSIO
#cmakedefine PLUGIN_GEOTOPIO
#cmakedefine PLUGIN_SMETIO
#cmakedefine PLUGIN_SNIO
#cmakedefine PLUGIN_PGMIO
#cmakedefine PLUGIN_IMISIO
#cmakedefine PLUGIN_GRIBIO
#cmakedefine PLUGIN_PNGIO
#cmakedefine PLUGIN_BORMAIO
#cmakedefine PLUGIN_COSMOXMLIO
#cmakedefine PLUGIN_GSNIO

#include <meteoio/plugins/ARCIO.h>
#include <meteoio/plugins/A3DIO.h>
#include <meteoio/plugins/ARPSIO.h>
#include <meteoio/plugins/GrassIO.h>
#include <meteoio/plugins/GeotopIO.h>
#include <meteoio/plugins/PGMIO.h>
#include <meteoio/plugins/SMETIO.h>
#include <meteoio/plugins/SNIO.h>

#ifdef PLUGIN_BORMAIO
#include <meteoio/plugins/BormaIO.h>
#endif

#ifdef PLUGIN_COSMOXMLIO
#include <meteoio/plugins/CosmoXMLIO.h>
#endif

#ifdef PLUGIN_IMISIO
#include <meteoio/plugins/ImisIO.h>
#endif

#ifdef PLUGIN_GRIBIO
#include <meteoio/plugins/GRIBIO.h>
#endif

#ifdef PLUGIN_GSNIO
#include <meteoio/plugins/GSNIO.h>
#endif

#ifdef PLUGIN_PNGIO
#include <meteoio/plugins/PNGIO.h>
#endif

using namespace std;

namespace mio {
 /**
 * @page plugins Plugins overview
 * The data access is handled by a system of plugins. They all offer the same interface, meaning that a plugin can transparently be replaced by another one. Since they might rely on third party libraries for accessing the data, they have been created as plugins, that is they are loaded on demand (and also compiled only if requested at compile time). A plugin can therefore fail to load (for example if it does not exist) at run time.
 *
 * @section available_categories Data sources categories
 * Several data sources categories have been defined that can be provided by a different plugin. Each data source category is defined by a specific key in the configuration file (usually, io.ini):
 * - METEO, for meteorological time series
 * - DEM, for Digital Elevation Maps
 * - LANDUSE, for land cover information
 * - GRID2D, for generic 2D grids (they can contain meteo fields and be recognized as such or arbitrary gridded data)
 * - SPECIALPTS, for a list of points that can be used for providing extra information at some specific location (extracting time series at a few selected points, etc)
 *
 * A plugin is "connected" to a given data source category simply by giving its keyword as value for the data source key:
 * @code
 * METEO = SMET
 * DEM = ARC
 * @endcode
 * Each plugin might have its own specific options, meaning that it might require its own keywords. Please check in each plugin documentation the supported options and keys (see links below).
 * Moreover, a given plugin might only support a given category for read or write (for example, PNG: there is no easy and safe way to interpret a given color as a given numeric value without knowing its color scale, so reading a png has been disabled).
 * Finally, the plugins usually don't implement all these categories (for example, ArcGIS file format only describes 2D grids, so the ARC plugin will only deal with 2D grids), so please check what a given plugin implements before connecting it to a specific data source category.
 *
 * @section available_plugins Available plugins
 * So far the following plugins have been implemented (by keyword for the io.ini key/value config file). Please read the documentation for each plugin in order to know the plugin-specific keywords:
 * <center><table border="1">
 * <tr><th>Plugin keyword</th><th>Provides</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><td>\subpage a3d "A3D"</td><td>meteo, specialpts</td><td>original Alpine3D meteo files</td><td></td></tr>
 * <tr><td>\subpage arc "ARC"</td><td>dem, landuse, grid2d</td><td>ESRI/ARC ascii grid files</td><td></td></tr>
 * <tr><td>\subpage arps "ARPS"</td><td>dem, grid2d</td><td>ARPS ascii formatted grids</td><td></td></tr>
 * <tr><td>\subpage borma "BORMA"</td><td>meteo</td><td>Borma xml meteo files</td><td><A HREF="http://libxmlplusplus.sourceforge.net/">libxml++</A></td></tr>
 * <tr><td>\subpage cosmoxml "COSMOXML"</td><td>meteo</td><td>MeteoSwiss COSMO's postprocessing XML format</td><td><A HREF="http://libxmlplusplus.sourceforge.net/">libxml++</A></td></tr>
 * <tr><td>\subpage geotop "GEOTOP"</td><td>meteo</td><td>GeoTop meteo files</td><td></td></tr>
 * <tr><td>\subpage grass "GRASS"</td><td>dem, landuse, grid2d</td><td>Grass grid files</td><td></td></tr>
 * <tr><td>\subpage gribio "GRIB"</td><td>meteo, dem, grid2d</td><td>GRIB meteo grid files</td><td><A HREF="http://www.ecmwf.int/products/data/software/grib_api.html">grib-api</A></td></tr>
 * <tr><td>\subpage gsn "GSN"</td><td>meteo</td><td>connects to the Global Sensor Network web service interface</td><td><SUP>1</SUP></td></tr>
 * <tr><td>\subpage imis "IMIS"</td><td>meteo</td><td>connects to the IMIS database</td><td><A HREF="http://docs.oracle.com/cd/B12037_01/appdev.101/b10778/introduction.htm">Oracle's OCCI library</A></td></tr>
 * <tr><td>\subpage pgmio "PGM"</td><td>dem, grid2d</td><td>PGM grid files</td><td></td></tr>
 * <tr><td>\subpage pngio "PNG"</td><td>dem, grid2d</td><td>PNG grid files</td><td><A HREF="http://www.libpng.org/pub/png/libpng.html">libpng</A></td></tr>
 * <tr><td>\subpage smetio "SMET"</td><td>meteo, specialpts</td><td>SMET data files</td><td></td></tr>
 * <tr><td>\subpage snowpack "SNOWPACK"</td><td>meteo</td><td>original SNOWPACK meteo files</td><td></td></tr>
 * </table></center>
 * <i><SUP>1</SUP>In order to rebuild the soap bindings for GSN, <A HREF="http://gsoap2.sourceforge.net/">gsoap</A> is required. This is only relevant to the plugin developers.</i>
 *
 * @section data_generators Data generators
 * It is also possible to duplicate a meteorological parameter as another meteorological parameter. This is done by specifying a COPY key, following the syntax
 * COPY::new_name = existing_parameter. For example:
 * @code
 * COPY::VW_avg = VW
 * @endcode
 * This creates a new parameter VW_avg that starts as an exact copy of the raw data of VW, for each station. This newly created parameter is
 * then processed as any other meteorological parameter (thus going through filtering, generic processing, spatial interpolations). This only current
 * limitation is that the parameter providing the raw data must be defined for all stations (even if filled with nodata, this is good enough).
 */

void IOHandler::registerPlugins()
{
	//mapPlugins[io.ini KEY]= IOPlugin(library file name, class name, NULL, NULL);
#ifdef PLUGIN_ARCIO
	mapPlugins["ARC"]       = IOPlugin("ARCIO", NULL, &IOPlugin::createInstance<ARCIO>);
#endif
#ifdef PLUGIN_A3DIO
	mapPlugins["A3D"]       = IOPlugin("A3DIO", NULL, &IOPlugin::createInstance<A3DIO>);
#endif
#ifdef PLUGIN_ARPSIO
	mapPlugins["ARPS"]      = IOPlugin("ARPSIO", NULL, &IOPlugin::createInstance<ARPSIO>);
#endif
#ifdef PLUGIN_GRASSIO
	mapPlugins["GRASS"]     = IOPlugin("GrassIO", NULL, &IOPlugin::createInstance<GrassIO>);
#endif
#ifdef PLUGIN_GEOTOPIO
	mapPlugins["GEOTOP"]    = IOPlugin("GeotopIO", NULL, &IOPlugin::createInstance<GeotopIO>);
#endif
#ifdef PLUGIN_SMETIO
	mapPlugins["SMET"]      = IOPlugin("SMETIO", NULL, &IOPlugin::createInstance<SMETIO>);
#endif
#ifdef PLUGIN_SNIO
	mapPlugins["SNOWPACK"]  = IOPlugin("SNIO", NULL, &IOPlugin::createInstance<SNIO>);
#endif
#ifdef PLUGIN_PGMIO
	mapPlugins["PGM"]       = IOPlugin("PGMIO", NULL, &IOPlugin::createInstance<PGMIO>);
#endif
#ifdef PLUGIN_IMISIO
	mapPlugins["IMIS"]      = IOPlugin("ImisIO", NULL, &IOPlugin::createInstance<ImisIO>);
#endif
#ifdef PLUGIN_GRIBIO
	mapPlugins["GRIB"]      = IOPlugin("GRIBIO", NULL, &IOPlugin::createInstance<GRIBIO>);
#endif
#ifdef PLUGIN_PNGIO
	mapPlugins["PNG"]       = IOPlugin("PNGIO", NULL, &IOPlugin::createInstance<PNGIO>);
#endif
#ifdef PLUGIN_BORMAIO
	mapPlugins["BORMA"]     = IOPlugin("BormaIO", NULL, &IOPlugin::createInstance<BormaIO>);
#endif
#ifdef PLUGIN_COSMOXMLIO
	mapPlugins["COSMOXML"]  = IOPlugin("CosmoXMLIO", NULL, &IOPlugin::createInstance<CosmoXMLIO>);
#endif
#ifdef PLUGIN_GSNIO
	mapPlugins["GSN"]       = IOPlugin("GSNIO", NULL, &IOPlugin::createInstance<GSNIO>);
#endif
}

//Copy constructor
#ifdef _POPC_
IOHandler::IOHandler(const IOHandler& aio)
           : cfg(aio.cfg), mapPlugins(), copy_parameter(aio.copy_parameter), copy_name(aio.copy_name), enable_copying(aio.enable_copying)
{
	//We might be on a different machine, even with different architecture -> rebuild plugins map
	registerPlugins();
}
#else
IOHandler::IOHandler(const IOHandler& aio)
           : IOInterface(), cfg(aio.cfg), mapPlugins(aio.mapPlugins), copy_parameter(aio.copy_parameter),
             copy_name(aio.copy_name), enable_copying(aio.enable_copying)
{

}
#endif

#ifdef _POPC_
IOHandler::IOHandler(const Config& cfgreader)
           : cfg(cfgreader), mapPlugins(), copy_parameter(), copy_name(), enable_copying(false)
#else
IOHandler::IOHandler(const Config& cfgreader)
           : IOInterface(), cfg(cfgreader), mapPlugins(), copy_parameter(), copy_name(), enable_copying(false)
#endif
{
	registerPlugins();
	parse_copy_config();
}

#ifdef _POPC_
IOHandler::~IOHandler(){
#else
IOHandler::~IOHandler() throw(){
#endif
	// Get rid of the objects
	std::map<std::string, IOPlugin>::iterator mapit;
	for (mapit = mapPlugins.begin(); mapit!=mapPlugins.end(); ++mapit){
		IOInterface*& io = (mapit->second).io;
		if (io != NULL) {
			delete io;
			io = NULL;
		}
	}
}

IOHandler& IOHandler::operator=(const IOHandler& source) {
	if(this != &source) {
		mapPlugins = source.mapPlugins;
		copy_parameter = source.copy_parameter;
		copy_name = source.copy_name;
		enable_copying = source.enable_copying;
	}
	return *this;
}

IOInterface* IOHandler::getPlugin(const std::string& cfgkey, const std::string& cfgsection)
{
	std::string op_src;
	cfg.getValue(cfgkey, cfgsection, op_src);

	const std::map<std::string, IOPlugin>::iterator mapit = mapPlugins.find(op_src);
	if (mapit == mapPlugins.end())
		throw IOException("Cannot find plugin " + op_src + " as requested in file " + cfg.getSourceName() + ". Has it been activated through ccmake? Is it registered in IOHandler::registerPlugins?", AT);
	if ((mapit->second).io == NULL){
		(mapit->second).io = (mapit->second).creator_func(cfg);
	}

	//Now check if it is correctly loaded
	if ((mapit->second).io == NULL) {
		throw IOException("Requesting to read/write data with plugin '" + op_src + "', but plugin is not loaded", AT);
	}

	return (mapit->second).io;
}

void IOHandler::read2DGrid(Grid2DObject& grid_out, const std::string& i_filename)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(grid_out, i_filename);
}

void IOHandler::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(grid_out, parameter, date);
}

void IOHandler::readDEM(DEMObject& dem_out)
{
	IOInterface *plugin = getPlugin("DEM", "Input");
	plugin->readDEM(dem_out);
	dem_out.update();
}

void IOHandler::readLanduse(Grid2DObject& landuse_out)
{
	IOInterface *plugin = getPlugin("LANDUSE", "Input");
	plugin->readLanduse(landuse_out);
}

void IOHandler::readStationData(const Date& date, STATIONS_SET& vecStation)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readStationData(date, vecStation);
}

void IOHandler::readMeteoData(const Date& date, METEO_SET& vecMeteo)
{
	std::vector< std::vector<MeteoData> > meteoTmpBuffer;
	readMeteoData(date, date, meteoTmpBuffer);

	vecMeteo.clear();
	vecMeteo.reserve( meteoTmpBuffer.size() );

	for (size_t ii=0; ii<meteoTmpBuffer.size(); ++ii) {//stations
		if (!meteoTmpBuffer[ii].empty())
			vecMeteo.push_back( meteoTmpBuffer[ii].front() );
	}
}

void IOHandler::checkTimestamps(const std::vector<METEO_SET>& vecVecMeteo) const
{
	for(size_t stat_idx=0; stat_idx<vecVecMeteo.size(); ++stat_idx) {
		const size_t nr_timestamps = vecVecMeteo[stat_idx].size();
		if(nr_timestamps==0) continue;

		Date previous_date( vecVecMeteo[stat_idx][0].date );
		for(size_t ii=1; ii<nr_timestamps; ++ii) {
			const StationData& station = vecVecMeteo[stat_idx][ii].meta;
			const Date& current_date = vecVecMeteo[stat_idx][ii].date;
			if(current_date<=previous_date) {
				throw IOException("Error at time "+current_date.toString(Date::ISO)+" for station \""+station.stationName+"\" ("+station.stationID+") : timestamps must be in increasing order and unique!", AT);
			}
			previous_date = current_date;
		}
	}
}

void IOHandler::readMeteoData(const Date& dateStart, const Date& dateEnd,
                              std::vector<METEO_SET>& vecMeteo,
                              const size_t& stationindex)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readMeteoData(dateStart, dateEnd, vecMeteo, stationindex);
	checkTimestamps(vecMeteo);

	copy_parameters(stationindex, vecMeteo);
}
#ifdef _POPC_
void IOHandler::writeMeteoData(std::vector<METEO_SET>& vecMeteo,
                               const std::string& name)
#else
void IOHandler::writeMeteoData(const std::vector<METEO_SET>& vecMeteo,
                               const std::string& name)
#endif
{
	IOInterface *plugin = getPlugin("METEO", "Output");
	plugin->writeMeteoData(vecMeteo, name);
}

void IOHandler::readAssimilationData(const Date& date_in, Grid2DObject& da_out)
{
	IOInterface *plugin = getPlugin("DA", "Input");
	plugin->readAssimilationData(date_in, da_out);
}

void IOHandler::readSpecialPoints(std::vector<Coords>& pts) {
	IOInterface *plugin = getPlugin("SPECIALPTS", "Input");
	plugin->readSpecialPoints(pts);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	IOInterface *plugin = getPlugin("GRID2D", "Output");
	plugin->write2DGrid(grid_in, name);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID2D", "Output");
	plugin->write2DGrid(grid_in, parameter, date);
}

void IOHandler::parse_copy_config()
{
	/**
	 * Parse [Input] section for potential parameters that the user wants
	 * duplicated (starting with 'COPY::')
	 */
	vector<string> copy_keys;
	const size_t nrOfMatches = cfg.findKeys(copy_keys, "COPY::", "Input");

	for (size_t ii=0; ii<nrOfMatches; ++ii) {
		const string name_of_copy = copy_keys[ii].substr(6);
		string initial_name;
		cfg.getValue(copy_keys[ii], "Input", initial_name);

		if ((name_of_copy.length() > 0) && (initial_name.length() > 0)){
			copy_parameter.push_back(initial_name);
			copy_name.push_back(name_of_copy);
			enable_copying = true;
		}
	}
}

/**
* This procedure runs through the MeteoData objects in vecMeteo and according to user
* configuration copies a certain present meteo parameter to another one, named by the
* user in the [Input] section of the io.ini, e.g.
* [Input]
* COPY::TA2 = TA
* means that TA2 will be the name of a new parameter in MeteoData with the copied value
* of the meteo parameter MeteoData::TA
*/
void IOHandler::copy_parameters(const size_t& stationindex, std::vector< METEO_SET >& vecMeteo) const
{
	if (!enable_copying) return; //Nothing configured

	size_t station_start=0, station_end=vecMeteo.size();
	if (stationindex != IOUtils::npos) {
		if (stationindex < vecMeteo.size()) {
			station_start = stationindex;
			station_end   = stationindex+1;
		} else {
			throw IndexOutOfBoundsException("Accessing stationindex in readMeteoData that is out of bounds", AT);
		}
	}

	const size_t nr_of_params = copy_parameter.size();
	vector<size_t> indices; //will hold the indices of the parameters to be copied

	for (size_t ii=station_start; ii<station_end; ++ii) { //for each station
		for (size_t jj=0; jj<vecMeteo[ii].size(); ++jj) { //for each MeteoData object of one station

			if (jj==0) { //buffer the index numbers
				for (size_t kk=0; kk<nr_of_params; ++kk) {
					const size_t param_index = vecMeteo[ii][jj].getParameterIndex(copy_parameter[kk]);
					if (param_index == IOUtils::npos) {
						std::ostringstream ss;
						ss << "At " << vecMeteo[ii][jj].date.toString(Date::ISO) << ", station " << vecMeteo[ii][jj].meta.stationID;
						ss << " has no parameter \"" << copy_parameter[kk] << "\" to copy!\n";
						throw InvalidArgumentException(ss.str(), AT);
					}

					indices.push_back(param_index);
				}
			}

			for (size_t kk=0; kk<nr_of_params; ++kk) {
				const size_t newparam = vecMeteo[ii][jj].addParameter(copy_name[kk]);
				vecMeteo[ii][jj](newparam) = vecMeteo[ii][jj](indices[kk]);
			}
		}
		indices.clear(); //may change for every station
	}
}


#ifdef _POPC_
const std::string IOHandler::toString()
#else
const std::string IOHandler::toString() const
#endif
{
	std::ostringstream os;
	os << "<IOHandler>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";

	os << "<mapPlugins>\n";
	os << setw(10) << "Keyword" << " = " << IOPlugin::header << "\n";
	std::map<std::string, IOPlugin>::const_iterator it1;
	for (it1=mapPlugins.begin(); it1 != mapPlugins.end(); ++it1){
		os << setw(10) << it1->first << " = " <<  it1->second.toString();
	}
	os << "</mapPlugins>\n";
	os << "</IOHandler>\n";
	return os.str();
}

} //end namespace
