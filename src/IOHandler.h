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
#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#ifdef _POPC_
#error
#endif

#include <map>
#include <string>

#include "IOInterface.h"
#include "A3DIO.h"
#include "IOExceptions.h"
#include "IOPlugin.h"

typedef std::map<std::string, IOPlugin::IOPlugin>::iterator PLUGIN_ITERATOR;

/**
* @file IOHandler.h
* The is the class implementing the interface as defined by the IOInterface class.
* This class is responsible for loading the necessary plugins and getting the data through them.
*/

class IOHandler : public IOInterface {
	public:
		IOHandler(const std::string& configfile);
		IOHandler(const IOHandler&);
		IOHandler(const ConfigReader&);
		~IOHandler() throw();

		//methods defined in the IOInterface class
		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readStationData(const Date_IO& date, std::vector<StationData>& vecStation);

		virtual void writeMeteoData(const std::vector<METEO_DATASET>& vecMeteo,
							   const std::vector<STATION_DATASET>& vecStation,
							   const std::string& name="");
		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd,
						std::vector<METEO_DATASET>& vecMeteo, 
						std::vector<STATION_DATASET>& vecStation,
						const unsigned& stationindex=IOUtils::npos);
		void readMeteoData(const Date_IO& date, METEO_DATASET& vecMeteo, STATION_DATASET& vecStation);
		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);

	private:
		void loadDynamicPlugins();
		void loadPlugin(const std::string& libname, const std::string& classname, 
					 DynamicLibrary*& dynLibrary, IOInterface*& io);
		void deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io) throw();
		void registerPlugins();
		IOInterface *getPlugin(const std::string& cfgkey, const std::string& cfgsection="GENERAL");

		ConfigReader cfg;
		std::map<std::string, IOPlugin::IOPlugin> mapPlugins;
		PLUGIN_ITERATOR mapit;
		A3DIO fileio;
};

#endif
