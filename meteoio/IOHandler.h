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

/*#ifdef _POPC_
#error
#endif*/

#include <meteoio/IOInterface.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOPlugin.h>

#include <map>
#include <string>

namespace mio {

/**
* @file IOHandler.h
* @class IOHandler
* @brief This class is the class to use for raw I/O operations. It is responsible for transparently loading the plugins
* and it follows the interface defined by the IOInterface class with the addition of a few convenience methods.
*/
#ifdef _POPC_
class IOHandler {
#else
class IOHandler : public IOInterface {
#endif
	public:
	#ifndef _POPC_
		IOHandler(const IOHandler&);
	#endif
		IOHandler(const Config&);
	#ifdef _POPC_
		virtual ~IOHandler();
	#else
		virtual ~IOHandler() throw();
	#endif

		//methods defined in the IOInterface class
		virtual void read2DGrid(Grid2DObject& out_grid, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readStationData(const Date& date,
		                             STATION_TIMESERIE& vecStation);
	#ifdef _POPC_
		virtual void writeMeteoData(std::vector<METEO_TIMESERIE>& vecMeteo,
		                            const std::string& name="");
	#else
		virtual void writeMeteoData(const std::vector<METEO_TIMESERIE>& vecMeteo,
		                            const std::string& name="");
	#endif
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector<METEO_TIMESERIE>& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);
		void readMeteoData(const Date& date, METEO_TIMESERIE& vecMeteo);

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	#ifndef _POPC_
		friend std::ostream& operator<<(std::ostream& os, const IOHandler& data);
	#endif
	#ifdef _POPC_
		std::string toString();
	#else
		std::string toString() const;
	#endif

	private:
		void registerPlugins();
		IOInterface *getPlugin(const std::string& cfgkey, const std::string& cfgsection="GENERAL");

		void parse_copy_config();
		void copy_parameters(const size_t& stationindex, std::vector< METEO_TIMESERIE >& vecMeteo) const;
		const Config& cfg;
		std::map<std::string, IOPlugin> mapPlugins;

		std::vector<std::string> copy_parameter, copy_name;
		bool enable_copying;
};

} //namespace

#endif
