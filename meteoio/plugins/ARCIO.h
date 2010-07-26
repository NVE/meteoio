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

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/DynamicLibrary.h>

#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
 * @class ARCIO
 * @brief This class enables the access to 2D grids stored in ESRI ASCII (ARCGIS) format
 *
 * @author Thomas Egger
 * @date   2009-08-28
 */

class ARCIO : public IOInterface {
	public:
		ARCIO(void (*delObj)(void*), const Config& i_cfg);

		ARCIO(const std::string& configfile);
		ARCIO(const ARCIO&);
		ARCIO(const Config&);
		~ARCIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);
		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, 
							   const std::vector< std::vector<StationData> >& vecStation,
							   const std::string& name="");
		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);

	private:
		void cleanup() throw();
		Config cfg;
		std::ifstream fin; //Input file streams
		std::ofstream fout;//Output file streams
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //end namespace mio

#endif
