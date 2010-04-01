/***********************************************************************************/
/*  Copyright 2009 HES-SO Fribourg                                                 */
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
#include <map>
#include <string>

#include "IOInterface.h"
#include "A3DIO.h"
#include "IOExceptions.h"
#include "IOPlugin.h"
#include "marshal_meteoio.h"

typedef std::map<std::string, IOPlugin::IOPlugin>::iterator PLUGIN_ITERATOR;

/**
* @file IOHandler.ph
* The is the parclass implementing the interface as defined by the IOInterface class.
* This class is responsible for loading the necessary plugins and getting the data through them.
*/

parclass IOHandler {
// Note : No heritage here for POPC++ : a parclass cannot herit from a class
		classuid(1003);
	public:
		IOHandler(const std::string& configfile) @{od.url("localhost");}; // @{ power=100 ?: 50; };
		//IOHandler(const IOHandler&) @{od.url("localhost");};
		IOHandler(const ConfigReader&) @{od.url("localhost");}; //@{ power=100 ?: 50; };
		~IOHandler();

		//methods defined in the IOInterface class
		virtual void read2DGrid([out]Grid2DObject& dem_out, const std::string& parameter="");
		virtual void readDEM([out]DEMObject& dem_out);
		virtual void readLanduse([out]Grid2DObject& landuse_out);
		virtual void readStationData([in]const Date_IO& date, 
			     	[proc=marshal_STATION_DATASET] STATION_DATASET& vecStation);
		virtual void writeMeteoData([in,proc=marshal_vector_METEO_DATASET] std::vector<METEO_DATASET>& vecMeteo,
			     [in,proc=marshal_vector_STATION_DATASET] std::vector<STATION_DATASET>& vecStation,
			     [in]const std::string& name);
		virtual void readMeteoData([in]const Date_IO& dateStart, [in]const Date_IO& dateEnd,
			[proc=marshal_vector_METEO_DATASET] std::vector<METEO_DATASET>& vecMeteo,
			[proc=marshal_vector_STATION_DATASET] std::vector<STATION_DATASET>& vecStation,
				const unsigned& stationindex=IOUtils::npos);
		void readMeteoData([in]const Date_IO& date, [proc=marshal_METEO_DATASET] METEO_DATASET& vecMeteo, [proc=marshal_STATION_DATASET] STATION_DATASET& vecStation);
		virtual void readAssimilationData([in] const Date_IO&,[out] Grid2DObject& da_out);
		virtual void readSpecialPoints([out,proc=marshal_vec_coords]std::vector<Coords>& pts);
		virtual void write2DGrid([in]const Grid2DObject& grid_in, [in]const std::string& name);

	private:
		void loadDynamicPlugins();
		void loadPlugin(const std::string& libname, const std::string& classname, 
					 DynamicLibrary*& dynLibrary, IOInterface*& io);
		void deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io);
		void registerPlugins();
		IOInterface *getPlugin(const std::string& cfgkey, const std::string& cfgsection="GENERAL");

		ConfigReader cfg;
		std::map<std::string, IOPlugin::IOPlugin> mapPlugins;
		PLUGIN_ITERATOR mapit;
		A3DIO fileio;
};

#endif
