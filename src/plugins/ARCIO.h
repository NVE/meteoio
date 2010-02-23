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
#ifndef __ARCIO_H__
#define __ARCIO_H__

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"
#include "Coords.h"
#include "IOExceptions.h"
#include "DynamicLibrary.h"

#include <string>
#include <sstream>
#include <iostream>


/**
 * @class ARCIO
 * @brief This class enables the access to 2D grids stored in ESRI ASCII (ARCGIS) format
 *
 * @author Thomas Egger
 * @date   2009-08-28
 */
class ARCIO : public IOInterface {
	public:
		ARCIO(void (*delObj)(void*), const std::string& filename);

		ARCIO(const std::string& configfile);
		ARCIO(const ARCIO&);
		ARCIO(const ConfigReader&);
		~ARCIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);

	private:
		void cleanup() throw();
		void getProjectionParameters();
		ConfigReader cfg;
		std::ifstream fin; //Input file streams
		std::ofstream fout;//Output file streams$
		std::string coordsys, coordparam; //projection parameters
};

#endif
