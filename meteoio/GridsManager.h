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
#ifndef __GRIDSMANAGER_H__
#define __GRIDSMANAGER_H__

#include <meteoio/DataGenerator.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOHandler.h>
#include <meteoio/Config.h>

namespace mio {

class GridsManager {
	public:
		GridsManager(IOHandler& in_iohandler, const Config& in_cfg);

		//Legacy support to support functionality of the IOInterface superclass:
		void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		void readDEM(DEMObject& dem_out);
		void readAssimilationData(const Date& date_in, Grid2DObject& da_out);
		void readLanduse(Grid2DObject& landuse_out);
		void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");
		void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);
		//end legacy support

		void setProcessingLevel(const unsigned int& i_level);
		void clear_cache();

		const std::string toString() const;

	private:
		void addToBuffer(const Grid2DObject& in_grid2Dobj, const std::string& grid_hash);
		bool getFromBuffer(const std::string& grid_hash, Grid2DObject& grid) const;

		IOHandler& iohandler;
		const Config& cfg;

		std::map<std::string, Grid2DObject> mapBufferedGrids;
		std::vector<DEMObject> dem_buffer;
		std::vector<std::string> IndexBufferedGrids; // this is required in order to know which grid is the oldest one
		size_t max_grids; ///< How many grids to buffer (grids, dem, landuse and assimilation grids together)
		unsigned int processing_level;
};
} //end namespace
#endif