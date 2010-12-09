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
#ifdef _POPC_
	#include <meteoio/IOHandler.ph>
	#ifdef DEBUG_ARITHM
		#ifndef _GNU_SOURCE
			#define _GNU_SOURCE
		#endif
		#include <fenv.h>
	#endif
#else
	#include <meteoio/IOHandler.h>
#endif

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
 * - \subpage imis "IMIS" for reading meteo data out of the IMIS database (requires Oracle's OCCI library)
 * - \subpage geotop "GEOTOP" for reading original GeoTop meteo files (no extra requirements)
 * - \subpage snowpack "SNOWPACK" for reading original SNOWPACK meteo files (no extra requirements)
 * - \subpage gsn "GSN" for reading meteo data out of the Global Sensor Network web service interface (requires GSoap)
 * - \subpage arc "ARC" for reading ESRI/ARC DEM files (no extra requirements)
 * - \subpage grass "GRASS" for reading Grass DEM files (no extra requirements)
 * - \subpage arps "ARPSIO" for reading ARPS formatted DEM (no extra requirements)
 * - \subpage pgmio "PGMIO" for reading PGM grid files (no extra requirements)
 * - \subpage smetio "SMETIO" for reading SMET meteo data files (no extra requirements)
 *
 */

void IOHandler::registerPlugins()
{
#if defined(WIN32)
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
	mapPlugins["SMET"]      = IOPlugin("libsmetio"+popc_extra+libsuffix, "SMETIO", NULL, NULL);
}

#ifdef _POPC_
IOHandler::IOHandler(const std::string& configfile) :  cfg(configfile), fileio(configfile){
#else
IOHandler::IOHandler(const std::string& configfile) : IOInterface(NULL), cfg(configfile), fileio(configfile)
{
#endif
	registerPlugins();
}

//Copy constructor
#ifdef _POPC_
//IOHandler::IOHandler(const IOHandler& aio) : cfg(aio.cfg), fileio(aio.cfg), bormaio(aio.cfg), imisio(aio.cfg){
	//Nothing else so far //HACK for POPC
//}
#else
IOHandler::IOHandler(const IOHandler& aio) : IOInterface(NULL), cfg(aio.cfg), fileio(aio.cfg)
{
	//Nothing else so far
	//TODO: Deal with the IOInterface* pointers, e.g. bormaio
}
#endif

#ifdef _POPC_
IOHandler::IOHandler(const Config& cfgreader) : cfg(cfgreader), fileio(cfgreader)
#else
IOHandler::IOHandler(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader), fileio(cfgreader)
#endif
{
#if defined(_POPC_) && defined(DEBUG_ARITHM)
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
#endif
	registerPlugins();
}

#ifdef _POPC_
IOHandler::~IOHandler(){
#else
IOHandler::~IOHandler() throw(){
#endif
	// Get rid of the objects
	std::map<std::string, IOPlugin::IOPlugin>::iterator mapit;
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
	std::cout << "[i] " << AT << ": Loading dynamic plugin: " << libname << std::endl;
	std::string pluginpath = "";

	try {
		cfg.getValue("PLUGINPATH", pluginpath, Config::nothrow);
		if (pluginpath != "")
			pluginpath += "/";

		//Which dynamic library needs to be loaded
		std::cout << "\t" << "Trying to load " << libname << " ... ";
		std::string filename = pluginpath + libname;
		dynLibrary = DynamicLoader::loadObjectFile(filename);

		if(dynLibrary == NULL) {
			std::cout << "failed\n\tCouldn't load the dynamic library " << filename << "\n\t" << DynamicLoader::getErrorMessage() << std::endl;
		} else {
			io = dynamic_cast<IOInterface*>((dynLibrary)->newObject(classname, cfg));

			if(io == NULL) {
				std::cout << "failed" << std::endl;
				//delete dynLibrary; This causes a segfault !!
			} else {
				std::cout << "success" << std::endl;
			}
		}
	} catch (std::exception& e) {
		if (dynLibrary != NULL)
			delete dynLibrary;
		std::cerr << "\t" << e.what() << std::endl;
	}
}

IOInterface* IOHandler::getPlugin(const std::string& cfgkey, const std::string& cfgsection)
{
	std::string op_src="";
	cfg.getValue(cfgkey, cfgsection, op_src);

	std::map<std::string, IOPlugin::IOPlugin>::iterator mapit = mapPlugins.find(op_src);
	if (mapit == mapPlugins.end())
		throw IOException("Can not find plugin " + op_src + " as requested in file " + cfg.getSourceName() + ". Has its developer declared it in IOHandler::registerPlugins?", AT);

	if ((mapit->second).io == NULL){
		loadPlugin((mapit->second).libname, (mapit->second).classname, (mapit->second).dynLibrary, (mapit->second).io);
	}

	if ((mapit->second).io == NULL) {
		throw IOException("Requesting to read/write data with plugin for " + cfgkey + ", but plugin is not loaded", AT);
	}

	return (mapit->second).io;
}

void IOHandler::read2DGrid(Grid2DObject& out_grid, const std::string& _filename)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(out_grid, _filename);
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

void IOHandler::readStationData(const Date& date, STATION_DATASET& vecStation)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readStationData(date, vecStation);
}

void IOHandler::readMeteoData(const Date& date, METEO_DATASET& vecMeteo)
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
						std::vector<METEO_DATASET>& vecMeteo,
						const unsigned& stationindex)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readMeteoData(dateStart, dateEnd, vecMeteo, stationindex);
}
#ifdef _POPC_
void IOHandler::writeMeteoData(std::vector<METEO_DATASET>& vecMeteo,
                               const std::string& name)
#else
void IOHandler::writeMeteoData(const std::vector<METEO_DATASET>& vecMeteo,
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
	os << cfg;

	os << "<mapPlugins>\n";
	os << setw(10) << "Keyword" << " = " << IOPlugin::header << "\n";
	std::map<std::string, IOPlugin::IOPlugin>::const_iterator it1;
	for (it1=mapPlugins.begin(); it1 != mapPlugins.end(); it1++){
		os << setw(10) << it1->first << " = " <<  it1->second;
	}
	os << "</mapPlugins>\n";
	os << "</IOHandler>\n";
	return os.str();
}

} //end namespace
