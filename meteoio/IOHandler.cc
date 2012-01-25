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

#include <meteoio/IOHandler.h>

using namespace std;

namespace mio {
 /**
 * @page plugins Plugins overview
 * The data access is handled by a system of plugins. They all offer the same interface, meaning that a plugin can transparently be replaced by another one. Since they might rely on third party libraries for accessing the data, they have been created as plugins, that is they are loaded on demand (and also compiled only if requested at compile time). A plugin can therefore fail to load (for example if it does not exist) at run time.
 *
 * @section available_plugins Available plugins
 * So far the following children have been implemented (by keyword for the io.ini key/value config file). Please read the documentation for each plugin in order to know the plugin-specific keywords:
 * - \subpage a3d "A3D" for reading original Alpine3D meteo files (no extra requirements)
 * - \subpage borma "BORMA" for reading Borma xml meteo files (requires libxml)
 * - \subpage cosmoxml "COSMO" for treading XML files produced by MeteoSwiss COSMO's postprocessing (requires libxml)
 * - \subpage imis "IMIS" for reading meteo data out of the IMIS database (requires Oracle's OCCI library)
 * - \subpage geotop "GEOTOP" for reading original GeoTop meteo files (no extra requirements)
 * - \subpage snowpack "SNOWPACK" for reading original SNOWPACK meteo files (no extra requirements)
 * - \subpage gsn "GSN" for reading meteo data out of the Global Sensor Network web service interface (requires GSoap)
 * - \subpage arc "ARC" for reading ESRI/ARC DEM files (no extra requirements)
 * - \subpage grass "GRASS" for reading Grass DEM files (no extra requirements)
 * - \subpage arps "ARPSIO" for reading ARPS formatted DEM (no extra requirements)
 * - \subpage pgmio "PGMIO" for reading PGM grid files (no extra requirements)
 * - \subpage pngio "PNGIO" for reading PNG grid files (requires libpng)
 * - \subpage smetio "SMETIO" for reading SMET meteo data files (no extra requirements)
 * - \subpage gribio "GribIO" for reading GRIB meteo grid files (requires libgrib)
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
#if defined(_WIN32)
	const std::string libsuffix = ".dll";
#elif defined(APPLE)
	const std::string libsuffix = ".dylib";
#else
	const std::string libsuffix = ".so";
#endif
#ifdef _POPC_
	const std::string popc_extra = "popc";
#else
	const std::string popc_extra = "";
#endif
	//mapPlugins[io.ini KEY]= IOPlugin(library file name, class name, NULL, NULL);
	mapPlugins["A3D"]       = IOPlugin("", "A3DIO", &fileio, NULL);
	mapPlugins["BORMA"]     = IOPlugin("libbormaio"+popc_extra+libsuffix, "BormaIO", NULL, NULL);
	mapPlugins["IMIS"]      = IOPlugin("libimisio"+popc_extra+libsuffix, "ImisIO", NULL, NULL);
	mapPlugins["GEOTOP"]    = IOPlugin("libgeotopio"+popc_extra+libsuffix, "GeotopIO", NULL, NULL);
	mapPlugins["SNOWPACK"]  = IOPlugin("libsnio"+popc_extra+libsuffix, "SNIO", NULL, NULL);
	mapPlugins["GSN"]       = IOPlugin("libgsnio"+popc_extra+libsuffix, "GSNIO", NULL, NULL);
	mapPlugins["ARC"]       = IOPlugin("libarcio"+popc_extra+libsuffix, "ARCIO", NULL, NULL);
	mapPlugins["GRASS"]     = IOPlugin("libgrassio"+popc_extra+libsuffix, "GrassIO", NULL, NULL);
	mapPlugins["ARPS"]      = IOPlugin("libarpsio"+popc_extra+libsuffix, "ARPSIO", NULL, NULL);
	mapPlugins["PGM"]       = IOPlugin("libpgmio"+popc_extra+libsuffix, "PGMIO", NULL, NULL);
	mapPlugins["PNG"]       = IOPlugin("libpngio"+popc_extra+libsuffix, "PNGIO", NULL, NULL);
	mapPlugins["SMET"]      = IOPlugin("libsmetio"+popc_extra+libsuffix, "SMETIO", NULL, NULL);
	mapPlugins["COSMOXML"]  = IOPlugin("libcosmoxmlio"+popc_extra+libsuffix, "CosmoXMLIO", NULL, NULL);
	mapPlugins["GRIB"]      = IOPlugin("libgribio"+popc_extra+libsuffix, "GRIBIO", NULL, NULL);
}

//Copy constructor
#ifdef _POPC_
//IOHandler::IOHandler(const IOHandler& aio) : cfg(aio.cfg), fileio(aio.cfg), bormaio(aio.cfg), imisio(aio.cfg){
	//Nothing else so far //HACK for POPC
//}
#else
IOHandler::IOHandler(const IOHandler& aio) : IOInterface(NULL), cfg(aio.cfg), fileio(aio.cfg), enable_copying(false)
{
	//Nothing else so far
	//TODO: Deal with the IOInterface* pointers, e.g. bormaio
}
#endif

#ifdef _POPC_
IOHandler::IOHandler(const Config& cfgreader) : cfg(cfgreader), fileio(cfgreader), enable_copying(false)
#else
IOHandler::IOHandler(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader), fileio(cfgreader), enable_copying(false)
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
	for (mapit = mapPlugins.begin(); mapit!=mapPlugins.end(); mapit++){
		deletePlugin((mapit->second).dynLibrary, (mapit->second).io);
	}
}
#ifdef _POPC_
void IOHandler::deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io)
#else
void IOHandler::deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io) throw()
#endif
{
	if (dynLibrary != NULL) {
		if (io != NULL) {
			io->deleteSelf();
			io = NULL;
		}

		// Close the dynamic library
#ifndef _POPC_ //HACK: this line causes a segfault in the parallel version for unknown reasons, lwk
		delete dynLibrary;
#endif
	}
}

void IOHandler::loadPlugin(const std::string& libname, const std::string& classname, DynamicLibrary*& dynLibrary, IOInterface*& io)
{
	std::string pluginpath = "";

	try {
		cfg.getValue("PLUGINPATH", pluginpath, Config::nothrow);
		if (pluginpath != "")
			pluginpath += "/";

		//Which dynamic library needs to be loaded
		std::string filename = pluginpath + libname;
		dynLibrary = DynamicLoader::loadObjectFile(filename);

		if(dynLibrary == NULL) {
			cout << AT << ": [E] Failed loading dynamic plugin " << classname << " from " << filename << std::endl;
			cout << "\t" << DynamicLoader::getErrorMessage() << std::endl;
			cout << "Please check your PLUGINPATH in your configuration file!" << std::endl;
		} else {
			io = dynamic_cast<IOInterface*>((dynLibrary)->newObject(classname, cfg));
			if(io == NULL) {
				cout << AT << ": [E] Failed loading dynamic plugin " << classname << " from " << filename << "(NULL pointer to plugin's class)" << std::endl;
				//delete dynLibrary; This causes a segfault !!
			} else {
				cout << "[i] Success loading dynamic plugin " << classname << " from " << filename << std::endl;
			}
		}
	} catch (const std::exception& e) {
	#ifndef _POPC_ //HACK: this line causes a segfault in the parallel version for unknown reasons, lwk
		if (dynLibrary != NULL)
			delete dynLibrary;
	#endif
		std::cerr << AT << ": [E] failed while loading plugin with error: \n" << e.what() << std::endl;
	}
}

IOInterface* IOHandler::getPlugin(const std::string& cfgkey, const std::string& cfgsection)
{
	std::string op_src="";
	cfg.getValue(cfgkey, cfgsection, op_src);

	std::map<std::string, IOPlugin>::iterator mapit = mapPlugins.find(op_src);
	if (mapit == mapPlugins.end())
		throw IOException("Can not find plugin " + op_src + " as requested in file " + cfg.getSourceName() + ". Has its developer declared it in IOHandler::registerPlugins?", AT);
	if ((mapit->second).io == NULL){
		loadPlugin((mapit->second).libname, (mapit->second).classname, (mapit->second).dynLibrary, (mapit->second).io);
	}
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

void IOHandler::readStationData(const Date& date, STATION_TIMESERIE& vecStation)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readStationData(date, vecStation);
}

void IOHandler::readMeteoData(const Date& date, METEO_TIMESERIE& vecMeteo)
{
	vecMeteo.clear();

	std::vector< std::vector<MeteoData> > meteoTmpBuffer;
	readMeteoData(date, date, meteoTmpBuffer);

	unsigned int emptycounter = 0;
	for (unsigned int ii=0; ii<meteoTmpBuffer.size(); ii++){//stations
		if (meteoTmpBuffer[ii].size() > 0){
			vecMeteo.push_back(meteoTmpBuffer[ii][0]);
		} else {
			emptycounter++;
		}
	}

	if (emptycounter == meteoTmpBuffer.size()){
		vecMeteo.clear();
	}
}

void IOHandler::readMeteoData(const Date& dateStart, const Date& dateEnd,
                              std::vector<METEO_TIMESERIE>& vecMeteo,
                              const size_t& stationindex)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readMeteoData(dateStart, dateEnd, vecMeteo, stationindex);

	copy_parameters(stationindex, vecMeteo);
}
#ifdef _POPC_
void IOHandler::writeMeteoData(std::vector<METEO_TIMESERIE>& vecMeteo,
                               const std::string& name)
#else
void IOHandler::writeMeteoData(const std::vector<METEO_TIMESERIE>& vecMeteo,
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
	size_t nrOfMatches = cfg.findKeys(copy_keys, "COPY::", "Input");

	for (size_t ii=0; ii<nrOfMatches; ii++) {
		string initial_name = "";
		string name_of_copy = copy_keys[ii].substr(6);
		cfg.getValue(copy_keys[ii], "Input", initial_name);

		if ((name_of_copy.length() > 0) && (initial_name.length() > 0)){
			copy_parameter.push_back(initial_name);
			copy_name.push_back(name_of_copy);
			//cout << "Param: " << initial_name << " copied to name: " << name_of_copy << endl;
			enable_copying = true;
		}
	}
}

void IOHandler::copy_parameters(const size_t& stationindex, std::vector< METEO_TIMESERIE >& vecMeteo) const
{
	/**
	 * This procedure runs through the MeteoData objects in vecMeteo and according to user
	 * configuration copies a certain present meteo parameter to another one, named by the
	 * user in the [Input] section of the io.ini, e.g.
	 * [Input]
	 * COPY::TA2 = TA
	 * means that TA2 will be the name of a new parameter in MeteoData with the copied value
	 * of the meteo parameter MeteoData::TA
	 */
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

	size_t nr_of_params = copy_parameter.size();
	vector<size_t> indices; //will hold the indices of the parameters to be copied

	for (size_t ii=station_start; ii<station_end; ii++) { //for each station
		for (size_t jj=0; jj<vecMeteo[ii].size(); jj++) { //for each MeteoData object of one station

			if (jj==0) { //buffer the index numbers
				for (size_t kk=0; kk<nr_of_params; kk++) {
					size_t param_index = vecMeteo[ii][jj].getParameterIndex(copy_parameter[kk]);
					if (param_index == IOUtils::npos) {
						std::stringstream ss;
						ss << "At " << vecMeteo[ii][jj].date.toString(Date::ISO) << ", station " << vecMeteo[ii][jj].meta.stationID;
						ss << " has no parameter \"" << copy_parameter[kk] << "\" to copy!\n";
						throw InvalidArgumentException(ss.str(), AT);
					}

					indices.push_back(param_index);
				}
			}

			for (size_t kk=0; kk<nr_of_params; kk++) {
				size_t newindex = vecMeteo[ii][jj].addParameter(copy_name[kk]);
				vecMeteo[ii][jj](newindex) = vecMeteo[ii][jj](indices[kk]);
			}
		}
		indices.clear(); //may change for every station
	}
}

#ifndef _POPC_
std::ostream& operator<<(std::ostream& os, const IOHandler& data)
{
	os << data.toString();
	return os;
}
#endif

#ifdef _POPC_
std::string IOHandler::toString()
#else
std::string IOHandler::toString() const
#endif
{
	std::stringstream os;
	os << "<IOHandler>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";

	os << "<mapPlugins>\n";
	os << setw(10) << "Keyword" << " = " << IOPlugin::header << "\n";
	std::map<std::string, IOPlugin>::const_iterator it1;
	for (it1=mapPlugins.begin(); it1 != mapPlugins.end(); it1++){
		os << setw(10) << it1->first << " = " <<  it1->second;
	}
	os << "</mapPlugins>\n";
	os << "</IOHandler>\n";
	return os.str();
}

} //end namespace
