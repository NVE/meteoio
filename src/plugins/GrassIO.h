/***********************************************************************************/
/*  Copyright 2009 EPFL                                                            */
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
#ifndef __GRASSIO_H__
#define __GRASSIO_H__

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"
#include "IOExceptions.h"
#include "Coords.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


/**
 * @class GrassIO
 * @brief This class enables the access to 2D grids stored in GRASS ASCII (e.g. JGrass) format
 *
 * @author Thomas Egger
 * @date   2008-08-03
 */

namespace mio {

class GrassIO : public IOInterface {
	public:
		GrassIO(void (*delObj)(void*), const std::string& filename);

		GrassIO(const std::string& configfile);
		GrassIO(const GrassIO&);
		GrassIO(const ConfigReader&);
		~GrassIO() throw();

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

		ConfigReader cfg;
		std::ifstream fin; //Input file streams
		std::ofstream fout;//Output file streams
		static const double plugin_nodata;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //end namespace mio

#endif
