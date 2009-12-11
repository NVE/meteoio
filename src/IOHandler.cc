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
	#include "IOHandler.ph"
	#ifdef DEBUG_ARITHM
		#ifndef _GNU_SOURCE
			#define _GNU_SOURCE
		#endif
		#include <fenv.h>
	#endif
#else
	#include "IOHandler.h"
#endif

 /**
 * @page plugins Plugins overview
 * The data access is handled by a system of plugins. They all offer the same interface, meaning that a plugin can transparently be replaced by another one. Since they might rely on third party libraries for accessing the data, they have been created as plugins, that is they are loaded on demand (and also compiled only if requested at compile time). A plugin can therefore fail to load (for example if it does not exist) at run time. 
 *
 * @section available_plugins Available plugins
 * So far the following children have been implemented (by keyword for the io.ini key/value config file). Please read the documentation for each plugin in order to know the plugin-specific keywords:
 * - \subpage a3d "A3D" for reading original Alpine3D meteo files (no extra requirements)
 * - \subpage boschung "BOSCHUNG" for reading Boshung xml meteo files (requires libxml)
 * - \subpage imis "IMIS" for reading meteo data out of the IMIS database (requires Oracle's OCCI library)
 * - \subpage geotop "GEOTOP" for reading original GeoTop meteo files (no extra requirements)
 * - \subpage gsn "GSN" for reading meteo data out of the Global Sensor Network web service interface (requires GSoap)
 * - \subpage arc "ARC" for reading ESRI/ARC DEM files (no extra requirements)
 * - \subpage grass "GRASS" for reading Grass DEM files (no extra requirements)
 * - \subpage arps "ARPS" for reading ARPS grid files (no extra requirements)
 * 
 */

void IOHandler::registerPlugins()
{
#if defined(WIN32)
	const std::string libsuffix = ".dll";
#elif defined(APPLE)
	const std::string libsuffix = ".dynlib";
#else
	const std::string libsuffix = ".so";
#endif
#ifdef _POPC_
	const std::string popc_extra = "popc";
#else
	const std::string popc_extra = "";
#endif
	mapPlugins["A3D"]	= IOPlugin("", "A3DIO", &fileio, NULL);
	mapPlugins["BOSCHUNG"]	= IOPlugin("libboschungio"+popc_extra+libsuffix, "BoschungIO", NULL, NULL);
	mapPlugins["IMIS"]	= IOPlugin("libimisio"+popc_extra+libsuffix, "ImisIO", NULL, NULL);
	mapPlugins["GEOTOP"]	= IOPlugin("libgeotopio"+popc_extra+libsuffix, "GeotopIO", NULL, NULL);
	mapPlugins["GSN"]	= IOPlugin("libgsnio"+popc_extra+libsuffix, "GSNIO", NULL, NULL);
	mapPlugins["ARC"]	= IOPlugin("libarcio"+popc_extra+libsuffix, "ARCIO", NULL, NULL);
	mapPlugins["GRASS"]	= IOPlugin("libgrassio"+popc_extra+libsuffix, "GrassIO", NULL, NULL);
}

#ifdef _POPC_
IOHandler::IOHandler(const std::string& configfile) :  cfg(configfile), fileio(cfg){
#else
IOHandler::IOHandler(const std::string& configfile) : IOInterface(NULL), cfg(configfile), fileio(cfg)
{
#endif
	registerPlugins();
}

//Copy constructor
#ifdef _POPC_
//IOHandler::IOHandler(const IOHandler& aio) : cfg(aio.cfg), fileio(cfg), boschungio(cfg), imisio(cfg){
	//Nothing else so far
//}
#else
IOHandler::IOHandler(const IOHandler& aio) : IOInterface(NULL), cfg(aio.cfg), fileio(cfg)
{
	//Nothing else so far
	//TODO: Deal with the IOInterface* pointers, e.g. boschungio
}
#endif

#ifdef _POPC_
IOHandler::IOHandler(const ConfigReader& cfgreader) : cfg(cfgreader), fileio(cfg)
#else
IOHandler::IOHandler(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader), fileio(cfg)
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
	cout << "[i] " << AT << ": Loading dynamic plugin: " << libname << endl;
	string pluginpath = "";

	try {
		cfg.getValue("PLUGINPATH", pluginpath, ConfigReader::nothrow);
		if (pluginpath != "")
			pluginpath += "/";

		//Which dynamic library needs to be loaded
		cout << "\t" << "Trying to load " << libname << " ... ";
		std::string filename = pluginpath + libname;
		dynLibrary = DynamicLoader::loadObjectFile(filename);
		
		if(dynLibrary == NULL) {
			cout << "failed\n\tCouldn't load the dynamic library " << filename << "\n\t" << DynamicLoader::getErrorMessage() << endl;
		} else {
			io = dynamic_cast<IOInterface*>((dynLibrary)->newObject(classname, cfg.getFileName()));

			if(io == NULL) {
				cout << "failed" << endl;
				//delete dynLibrary; This causes a segfault !!
			} else {
				cout << "success" << endl;
			}
		}
	} catch (exception& e) {
		if (dynLibrary != NULL)
		delete dynLibrary;
		cerr << "\t" << e.what() << endl;
	}
}

IOInterface* IOHandler::getPlugin(const std::string& cfgvalue)
{
	std::string op_src="";
	cfg.getValue(cfgvalue, op_src);

	mapit = mapPlugins.find(op_src);
	if (mapit == mapPlugins.end())
		throw IOException(cfgvalue + " does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
	
	if ((mapit->second).io == NULL){
		loadPlugin((mapit->second).libname, (mapit->second).classname, (mapit->second).dynLibrary, (mapit->second).io);
	}
	
	if ((mapit->second).io == NULL) {
		throw IOException("Requesting to read/write data with plugin for " + cfgvalue + ", but plugin is not loaded", AT);
	}

	return (mapit->second).io;
}

void IOHandler::read2DGrid(Grid2DObject& _grid, const std::string& _filename)
{
	IOInterface *plugin = getPlugin("GRID2DSRC");
	plugin->read2DGrid(_grid, _filename);
}

void IOHandler::readDEM(DEMObject& dem_out)
{
	IOInterface *plugin = getPlugin("DEMSRC");
	plugin->readDEM(dem_out);
	dem_out.update();
}

void IOHandler::readLanduse(Grid2DObject& landuse_out)
{
	IOInterface *plugin = getPlugin("LANDUSESRC");
	plugin->readLanduse(landuse_out);
}

void IOHandler::readMeteoData(const Date_IO& date, METEO_DATASET& vecMeteo, STATION_DATASET& vecStation)
{
	vecMeteo.clear();
	vecStation.clear();
	
	std::vector< std::vector<MeteoData> > meteoTmpBuffer;
	std::vector< std::vector<StationData> > stationTmpBuffer;
	readMeteoData(date, date, meteoTmpBuffer, stationTmpBuffer);

	unsigned int emptycounter = 0;
	for (unsigned int ii=0; ii<meteoTmpBuffer.size(); ii++){//stations
		if ((meteoTmpBuffer[ii].size() > 0) && (stationTmpBuffer[ii].size() > 0)){
			vecMeteo.push_back(meteoTmpBuffer[ii][0]);
			vecStation.push_back(stationTmpBuffer[ii][0]);
		} else {
			emptycounter++;
		}
	}

	if (emptycounter == meteoTmpBuffer.size()){
		vecMeteo.clear();
		vecStation.clear();
	}
}

void IOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
						std::vector<METEO_DATASET>& vecMeteo, 
						std::vector<STATION_DATASET>& vecStation, 
						const unsigned& stationindex)
{
	IOInterface *plugin = getPlugin("METEOSRC");
	plugin->readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
}

void IOHandler::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	IOInterface *plugin = getPlugin("DASRC");
	plugin->readAssimilationData(date_in, da_out);
}

void IOHandler::readSpecialPoints(CSpecialPTSArray& pts) {
	IOInterface *plugin = getPlugin("SPECIALPTSSRC");
	plugin->readSpecialPoints(pts);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	IOInterface *plugin = getPlugin("OUTPUT");
	plugin->write2DGrid(grid_in, name);
}


