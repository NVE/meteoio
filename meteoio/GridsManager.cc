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

#include <meteoio/GridsManager.h>

using namespace std;

namespace mio {

GridsManager::GridsManager(IOHandler& in_iohandler, const Config& in_cfg)
             : iohandler(in_iohandler), cfg(in_cfg), mapBufferedGrids(), dem_buffer(), IndexBufferedGrids(), max_grids(10),
               processing_level(IOUtils::filtered | IOUtils::resampled | IOUtils::generated)
{
	cfg.getValue("BUFF_GRIDS", "General", max_grids, IOUtils::nothrow);
}

/**
* @brief Set the desired ProcessingLevel
*        The processing level affects the way meteo data is read and processed
*        Three values are possible:
*        - IOUtils::raw data shall be read directly from the buffer
*        - IOUtils::filtered data shall be filtered before returned to the user
*        - IOUtils::resampled data shall be resampled before returned to the user
*          this only affects the function getMeteoData(const Date&, METEO_DATASET&);
*
*        The three values can be combined: e.g. IOUtils::filtered | IOUtils:resampled
* @param i_level The ProcessingLevel values that shall be used to process data
*/
void GridsManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOUtils::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

	if (((i_level & IOUtils::raw) == IOUtils::raw)
	    && ((i_level & IOUtils::filtered) == IOUtils::filtered))
		throw InvalidArgumentException("The processing level is invalid (raw and filtered at the same time)", AT);

	processing_level = i_level;
}

void GridsManager::clear_cache()
{
	mapBufferedGrids.clear();
	IndexBufferedGrids.clear();
}

void GridsManager::read2DGrid(Grid2DObject& grid2D, const std::string& filename)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, filename);
	} else {
		if (getFromBuffer(filename, grid2D))
			return;

		iohandler.read2DGrid(grid2D, filename);
		addToBuffer(grid2D, filename);
	}
}

void GridsManager::read2DGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, parameter, date);
	} else {
		if (max_grids>0) {
			const string grid_hash = date.toString(Date::ISO)+"::"+MeteoGrids::getParameterName(parameter);
			if (getFromBuffer(grid_hash, grid2D))
				return;

			iohandler.read2DGrid(grid2D, parameter, date);
			addToBuffer(grid2D, grid_hash); //the STL containers make a copy
		} else {
			iohandler.read2DGrid(grid2D, parameter, date);
		}
	}
}

void GridsManager::readDEM(DEMObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readDEM(grid2D);
	} else {
		if (max_grids>0) {
			if (dem_buffer.size() == 1) { //HACK: properly manage buffering multiple dems!
				//already in buffer. If the update properties have changed,
				//we copy the ones given in input and force the update of the object
				const DEMObject::update_type in_ppt = (DEMObject::update_type)grid2D.getUpdatePpt();
				const DEMObject::slope_type in_slope_alg = (DEMObject::slope_type)grid2D.getDefaultAlgorithm();

				grid2D = dem_buffer[0];
				const DEMObject::update_type buff_ppt = (DEMObject::update_type)grid2D.getUpdatePpt();
				const DEMObject::slope_type buff_slope_alg = (DEMObject::slope_type)grid2D.getDefaultAlgorithm();

				if (in_ppt!=buff_ppt || in_slope_alg!=buff_slope_alg) {
					grid2D.setDefaultAlgorithm(in_slope_alg);
					grid2D.setUpdatePpt(in_ppt);
					grid2D.update();
				}

				return;
			}

			iohandler.readDEM(grid2D);
			dem_buffer.push_back(grid2D); //the STL containers make a copy
		} else {
			iohandler.readDEM(grid2D);
		}
	}
}

void GridsManager::readLanduse(Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readLanduse(grid2D);
	} else {
		if (getFromBuffer("/:LANDUSE", grid2D))
			return;

		iohandler.readLanduse(grid2D);
		addToBuffer(grid2D, "/:LANDUSE");
	}
}

void GridsManager::readAssimilationData(const Date& date, Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readAssimilationData(date, grid2D);
	} else {
		if(max_grids>0) {
			const string grid_hash = "/:ASSIMILATIONDATA"+date.toString(Date::ISO);
			if (getFromBuffer(grid_hash, grid2D))
				return;

			iohandler.readAssimilationData(date, grid2D);
			addToBuffer(grid2D, grid_hash); //the STL containers make a copy
		} else {
			iohandler.readAssimilationData(date, grid2D);
		}
	}
}

void GridsManager::write2DGrid(const Grid2DObject& grid2D, const std::string& name)
{
	iohandler.write2DGrid(grid2D, name);
}

void GridsManager::write2DGrid(const Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write2DGrid(grid2D, parameter, date);
}

void GridsManager::addToBuffer(const Grid2DObject& in_grid2Dobj, const std::string& grid_hash)
{
	if (max_grids==0) return;

	if(IndexBufferedGrids.size() >= max_grids) { //we need to remove the oldest grid
		mapBufferedGrids.erase( mapBufferedGrids.find( IndexBufferedGrids.front() ) );
		IndexBufferedGrids.erase( IndexBufferedGrids.begin() );
	}
	mapBufferedGrids[ grid_hash ] = in_grid2Dobj;
	IndexBufferedGrids.push_back( grid_hash  );
}

bool GridsManager::getFromBuffer(const std::string& grid_hash, Grid2DObject& grid) const
{
	if (IndexBufferedGrids.empty())
		return false;

	const std::map<std::string, Grid2DObject>::const_iterator it = mapBufferedGrids.find( grid_hash );
	if (it != mapBufferedGrids.end()) { //already in map
		grid = (*it).second;
		return true;
	}

	return false;
}

const std::string GridsManager::toString() const {
	ostringstream os;
	os << "<GridsManager>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOHandler& iohandler = " << hex << &iohandler << dec << "\n";
	os << "Processing level = " << processing_level << "\n";
	os << "Max buffered grids = " << max_grids << "\n";

	//cache content
	os << "Current buffer content (" << mapBufferedGrids.size() << " grids):\n";
	std::map<std::string, Grid2DObject>::const_iterator it1;
	for (it1=mapBufferedGrids.begin(); it1 != mapBufferedGrids.end(); ++it1){
		os << setw(10) << "Grid " << it1->first << "\n";
	}

	//dem buffer
	os << "Dem buffer content (" << dem_buffer.size() << " grids):\n";
	for (size_t ii=0; ii<dem_buffer.size(); ++ii){
		std::ostringstream ss;
		ss << dem_buffer[ii].llcorner.printLatLon() << " " << dem_buffer[ii].getNx() << "x" << dem_buffer[ii].getNy() << " @" << dem_buffer[ii].cellsize << "m";
		os << setw(10) << "Dem " << ss.str() << "\n";
	}

	os << "</GridsManager>\n";
	return os.str();
}

} //namespace
