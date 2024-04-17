// SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <meteoio/plugins/GRIBIO_redesign.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/FileUtils.h>


namespace mio {
	using namespace codes;
/**
 * @page gribio GRIBIO
 * @section gribio_format Format
 * *Put here the information about the standard format that is implemented*
 *
 * @section gribio_units Units
 *
 *
 * @section gribio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 * 
 * @note Grid2d and Meteo Files need to be seperated by either pattern or extension
 */

const double GRIBIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const std::string GRIBIO::default_table = "doc/resources/GRIB_param.tbl"; 

#include <algorithm>

static size_t findDate(const std::vector<GRIBFile>& cache, const Date& date) {
    auto it = std::find_if(cache.begin(), cache.end(), [&date](const GRIBFile& file) { return file.isValidDate(date); });

    return (it != cache.end()) ? std::distance(cache.begin(), it) : IOUtils::npos;
}

// ----------------------------- INITIALIZE -----------------------------
GRIBIO::GRIBIO(const std::string& configfile) : cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initialize();
}

GRIBIO::GRIBIO(const Config& cfgreader) : cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initialize();
}

void GRIBIO::initialize() {
	setOptions();
	initTable();
	scanPath(grid2dpath_in, grid2d_ext, grid_2d_pattern, cache_grid2d);
	scanPath(meteopath_in, meteo_ext, meteo_pattern, cache_meteo);
}

void GRIBIO::setOptions()
{
	const std::string tmp = IOUtils::strToUpper(cfg.get("GRID2D", "Input", ""));
	if (tmp == "GRIB") { // keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
		cfg.getValue("GRIB_TABLE", "Input", table_path);

		cfg.getValue("GRIB_DEM_UPDATE", "Input", update_dem, IOUtils::nothrow);
	}

	cfg.getValue("METEOEXT", "Input", meteo_ext, IOUtils::nothrow);
	if (meteo_ext == "none")
		meteo_ext.clear();
	cfg.getValue("METEOPATTERN", "Input", meteo_pattern, IOUtils::nothrow);
	if (meteo_ext == "none")
		meteo_ext.clear();

	cfg.getValue("GRID2DEXT", "Input", grid2d_ext, IOUtils::nothrow);
	if (grid2d_ext == "none")
		grid2d_ext.clear();
	cfg.getValue("GRID2DPATTERN", "Input", grid_2d_pattern, IOUtils::nothrow);
	if (grid2d_ext == "none")
		grid2d_ext.clear();

	cfg.getValue("VERBOSE", "Input", verbose, IOUtils::nothrow);
	cfg.getValue("RECURSIVE", "Input", recursive_search, IOUtils::nothrow);

	if (meteo_ext == grid2d_ext && meteo_pattern == grid_2d_pattern && meteopath_in == grid2dpath_in)
		throw InvalidArgumentException("Meteo and Grid2D files cannot have the same naming and be located in the same place.", AT);
}

void GRIBIO::initTable() {

        if (table_path.empty()) {
#ifdef DEBUG
            std::cerr << "No GRIB table specified, using default table" << std::endl;
#endif
            table_path = default_table;
        }
		parameter_table = GRIBTable(table_path);
}

void GRIBIO::scanPath(const std::string &in_path, const std::string &in_ext, const std::string &in_pattern, std::vector<GRIBFile> &cache) {
	std::list<std::string> dirlist;
	FileUtils::readDirectory(in_path, dirlist, in_ext, recursive_search);

	std::list<std::string>::const_iterator it = dirlist.begin();
	while ((it != dirlist.end())) {
		const std::string filename = *it;
		if (filename.find(grid_2d_pattern) != std::string::npos) {
			const std::string fullpath = grid2dpath_in + "/" + filename;
			cache_grid2d.push_back(GRIBFile(fullpath, parameter_table.getIndexes())); // Files will not be read twice, as we enforce a seperation between meteo and grid
		}
		it++;
	}
}


void GRIBIO::readStations(std::vector<Coords> &vecPoints) {
	cfg.getValue("METEOPATH", "Input", meteopath_in);

	std::vector<std::string> vecStation;
	cfg.getValues("STATION", "INPUT", vecStation);
	for (size_t ii = 0; ii < vecStation.size(); ii++) {
		Coords tmp(coordin, coordinparam, vecStation[ii]);
		if (!tmp.isNodata())
			vecPoints.push_back(tmp);

		std::cerr << "\tRead virtual station " << vecPoints.back().toString(Coords::LATLON) << "\n";
	}
}

Coords GRIBIO::getGeolocalization(double &cellsize_x, double &cellsize_y, const std::map<std::string,double> &grid_params) {
		latitudeOfNorthernPole = grid_params.at("latitudeOfNorthernPole");
        longitudeOfNorthernPole = grid_params.at("longitudeOfNorthernPole");
        double ll_latitude = grid_params.at("ll_latitude");
        double ll_longitude = grid_params.at("ll_longitude");
        double ur_latitude = grid_params.at("ur_latitude");
        double ur_longitude = grid_params.at("ur_longitude");
        double angleOfRotationInDegrees = grid_params.at("angleOfRotationInDegrees");
        double Ni = grid_params.at("Ni");
        double Nj = grid_params.at("Nj");

        if (angleOfRotationInDegrees != 0.) {
            throw InvalidArgumentException("Rotated grids not supported!", AT);
        }

        double ur_lat, ur_lon, ll_lat, ll_lon;
        CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, ur_latitude, ur_longitude, ur_lat, ur_lon);
        double cntr_lat, cntr_lon; // geographic coordinates
        CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, .5 * (ll_latitude + ur_latitude), .5 * (ll_longitude + ur_longitude), cntr_lat, cntr_lon);

        double bearing;
        cellsize_x = CoordsAlgorithms::VincentyDistance(cntr_lat, ll_lon, cntr_lat, ur_lon, bearing) / (double)Ni;
        cellsize_y = CoordsAlgorithms::VincentyDistance(ll_lat, cntr_lon, ur_lat, cntr_lon, bearing) / (double)Nj;

        // determining bearing offset
        double delta_lat, delta_lon; // geographic coordinates
        CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, .5 * (ll_latitude + ur_latitude) + 1., .5 * (ll_longitude + ur_longitude), delta_lat, delta_lon);
        CoordsAlgorithms::VincentyDistance(cntr_lat, cntr_lon, delta_lat, delta_lon, bearing_offset);
        bearing_offset = fmod(bearing_offset + 180., 360.) - 180.; // turn into [-180;180)

        // returning the center point as reference
        Coords cntr(coordin, coordinparam);
        cntr.setLatLon(cntr_lat, cntr_lon, IOUtils::nodata);

        return cntr;
}


// ----------------------------- GRIDDED DATA -----------------------------
    void GRIBIO::read2Dlevel(CodesHandlePtr &h, Grid2DObject &grid_out, const std::map<std::string, double> &grid_params ) {
        std::vector<double> values;
        getGriddedValues(h, values);

        double Ni = grid_params.at("Ni");
        double Nj = grid_params.at("Nj");

        if (!llcorner_initialized) {
            // most important: get cellsize. llcorner will be finalized AFTER aspect ration correction
            double cellsize_x, cellsize_y;
            llcorner = getGeolocalization(cellsize_x, cellsize_y, grid_params); // this is the center cell

            cellsize = (double)Optim::round(std::min(cellsize_x, cellsize_y) * 100.) / 100.; // round to 1cm precision for numerical stability
            if (std::abs(cellsize_x - cellsize_y) / cellsize_x > 1. / 100.) {
                factor_x = cellsize_x / cellsize;
                factor_y = cellsize_y / cellsize;
            }
        }

        grid_out.set(static_cast<size_t>(Ni), static_cast<size_t>(Nj), cellsize, llcorner);
        size_t i = 0;
        for (size_t jj = 0; jj < (unsigned)Nj; jj++) {
            for (size_t ii = 0; ii < (unsigned)Ni; ii++)
                grid_out(ii, jj) = values[i++];
        }

        // cells were not square, we have to resample
        if (factor_x != IOUtils::nodata && factor_y != IOUtils::nodata) {
            grid_out.grid2D = LibResampling2D::Bilinear(grid_out.grid2D, factor_x, factor_y);
        }

        if (!llcorner_initialized) { // take into account aspect ration conversion for computing true llcorner
            llcorner.moveByXY(-.5 * (double)grid_out.getNx() * cellsize, -.5 * (double)grid_out.getNy() * cellsize);
            llcorner_initialized = true;

            grid_out.llcorner = llcorner;
        }
    }

void GRIBIO::read2DGrid(Grid2DObject& grid_out, const std::string& parameter) { // TODO: do we actually need this? And how do i specify where and what?
	throw IOException("Nothing implemented here", AT);
};


std::vector<CodesHandlePtr> GRIBIO::findMessages(const std::string& param_key, const std::string& level_key, const std::string& level_type, const std::string& paramID_string, double paramID_double, long paramID_long, int idx) {
    double npos_double = static_cast<double>(IOUtils::npos);
    long npos_long = static_cast<long>(IOUtils::npos);

    if (!paramID_string.empty() && paramID_double == npos_double && paramID_long == npos_long) {
        return cache_grid2d[idx].listParameterMessages(param_key, paramID_string, level_key, level_type);
    } else if (paramID_string.empty() && paramID_double != npos_double && paramID_long == npos_long) {
        return cache_grid2d[idx].listParameterMessages(param_key, paramID_double, level_key, level_type);
    } else if (paramID_string.empty() && paramID_double == npos_double && paramID_long != npos_long) {
        return cache_grid2d[idx].listParameterMessages(param_key, paramID_long, level_key, level_type);
    } else {
        throw IOException("paramID needs to be either string, double or long", AT);
    }
}


void GRIBIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date) {
	size_t idx = findDate(cache_grid2d, date);
	if (idx == IOUtils::npos) {
		throw IOException("No grid found for the specified date", AT); // TODO: do we need to throw an exception here, or just return?
	}

	std::string param_name = MeteoGrids::getParameterName(parameter);

	// get the paramID	
	std::string paramID_string;
	double paramID_double;
	long paramID_long;
	parameter_table.getParamId(param_name, paramID_string, paramID_double, paramID_long);

	// get additional information from the parameter table
	std::string param_key = parameter_table.getParamKey();
	std::string level_key = parameter_table.getLevelKey();
	long level_no = parameter_table.getLevelNo(param_name);
	std::string level_type = parameter_table.getLevelType(param_name);

	std::vector<CodesHandlePtr> messages = findMessages(param_key, level_key, level_type, paramID_string, paramID_double, paramID_long, idx);

	if (messages.empty()) {
		throw IOException("No messages found for the specified parameter", AT); // TODO: do we need to throw an exception here, or just return?
	}
	
	for (auto &m : messages) {
		long level = 0;
		if (level_no != 0)
			getParameter(m, "level", level);
		if (level == level_no) {
			read2Dlevel(m, grid_out, cache_grid2d[idx].getGridParams());
			return;
		}
	}
	throw IOException("No messages found for the specified level", AT); // TODO: do we need to throw an exception here, or just return?
	return;
};


// ---------------------------- STATION DATA ----------------------------- 
void GRIBIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}



} //namespace
