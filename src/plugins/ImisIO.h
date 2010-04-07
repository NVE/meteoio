/***********************************************************************************/
/*  Copyright 2009, 2010 WSL Institute for Snow and Avalanche Research   SLF-DAVOS */
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
#ifndef __IMISIO_H__
#define __IMISIO_H__


#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <occi.h>

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"
#include "DynamicLibrary.h"
#include "MeteoData.h"
#include "StationData.h"
#include "IOExceptions.h"

#define IMIS_BUFF_SIZE 5000

class ImisIO : public IOInterface {
	public:
		ImisIO(void (*delObj)(void*), const std::string& filename);

		ImisIO(const std::string& configfile);
		ImisIO(const ImisIO&);
		ImisIO(const ConfigReader&);
		~ImisIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date_IO& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           std::vector< std::vector<StationData> >& vecStation,
		                           const unsigned int& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, 
		                            const std::vector< std::vector<StationData> >& vecStation,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);

	private:
		void cleanup() throw();
		void getProjectionParameters();
		void getStation2Data(const std::string stat_abk, unsigned int stao_nr, std::vector<std::string>& data2S);
		void getImisData(const std::string &stat_abk, const unsigned int &stao_nr, 
		                 const std::vector<int>& datestart, const std::vector<int>& dateend,
		                 std::vector< std::vector<std::string> >& dataImis);
		void parseDataSet(const std::vector<std::string>& meteo_in, MeteoData& md);		
		void readData(const Date_IO& dateStart, const Date_IO& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, 
		              std::vector< std::vector<StationData> >& vecStation, const unsigned int& stationindex);
		void readStationNames(std::vector<std::string>& vecStationName);
		void parseStationName(const std::string& stationName, std::string& stName, unsigned int& stNumber);
		void readStationMetaData();
		void convertUnits(MeteoData& meteo);

		ConfigReader cfg;
		std::string coordsys, coordparam; //projection parameters
		std::vector<StationData> vecMyStation;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
};

#endif

