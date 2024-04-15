// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/GRIBIO.h>

#include <meteoio/FileUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/meteoLaws/Meteoconst.h> //for PI
#include <meteoio/meteoStats/libresampling2D.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <grib_api.h>
#include <iostream>

using namespace std;

namespace mio {

    using namespace codes;
    /**
     * @page gribio GRIBIO
     * @section gribio_format Format and limitations
     *
     * @note
     * - why do we need to have the time information in the filename, and not just in the file itself?
     * - should it be possible to read multiple dates from a single file?
     *
     * @note
     * TODO: temperature vs temperature2m when should be what, with surface and height above ground...
     * 		 add reading arbitrary parameter -->> shouldnt meteogrids have an addParameter as well? -> needs to be done in the handler, when grids are seperate from meteodata
     *
     */

    const double GRIBIO::plugin_nodata = -999.;     // plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
    const double GRIBIO::tz_in = 0.;                // GRIB time zone, always UTC
    const std::string GRIBIO::default_ext = ".grb"; // filename extension

    GRIBIO::GRIBIO(const std::string &configfile)
        : cfg(configfile), grid2dpath_in(), meteopath_in(), vecPts(), cache_meteo_files(), meteo_ext(default_ext), grid2d_ext(default_ext), grid2d_prefix(), coordin(), coordinparam(), VW(), DW(),
          wind_date(), llcorner(), file_index(nullptr), params_to_read(), level_types_to_read(), levels_to_read(), ensemble_to_read(-1), verbose(false), paramIdList(),
          latitudeOfNorthernPole(IOUtils::nodata), longitudeOfNorthernPole(IOUtils::nodata), bearing_offset(IOUtils::nodata), cellsize(IOUtils::nodata), factor_x(IOUtils::nodata),
          factor_y(IOUtils::nodata), meteo_initialized(false), llcorner_initialized(false), update_dem(false), indexed(false), warned_date_in_file(false) {
        setOptions();
    }

    GRIBIO::GRIBIO(const Config &cfgreader)
        : cfg(cfgreader), grid2dpath_in(), meteopath_in(), vecPts(), cache_meteo_files(), meteo_ext(default_ext), grid2d_ext(default_ext), grid2d_prefix(), coordin(), coordinparam(), VW(), DW(),
          wind_date(), llcorner(), file_index(nullptr), params_to_read(), level_types_to_read(), levels_to_read(), ensemble_to_read(-1), verbose(false), paramIdList(),
          latitudeOfNorthernPole(IOUtils::nodata), longitudeOfNorthernPole(IOUtils::nodata), bearing_offset(IOUtils::nodata), cellsize(IOUtils::nodata), factor_x(IOUtils::nodata),
          factor_y(IOUtils::nodata), meteo_initialized(false), llcorner_initialized(false), update_dem(false), indexed(false), warned_date_in_file(false) {
        setOptions();
    }

    void GRIBIO::setOptions() {
        std::string coordout, coordoutparam;
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

        const std::string tmp = IOUtils::strToUpper(cfg.get("GRID2D", "Input", ""));
        if (tmp == "GRIB") { // keep it synchronized with IOHandler.cc for plugin mapping!!
            cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
            cfg.getValue("GRIB_DEM_UPDATE", "Input", update_dem, IOUtils::nothrow);
        }
        cfg.getValue("GRID2DPREFIX", "Input", grid2d_prefix, IOUtils::nothrow);

        cfg.getValue("METEOEXT", "Input", meteo_ext, IOUtils::nothrow);
        if (meteo_ext == "none")
            meteo_ext.clear();

        cfg.getValue("GRID2DEXT", "Input", grid2d_ext, IOUtils::nothrow);
        if (grid2d_ext == "none")
            grid2d_ext.clear();

        cfg.getValue("GRIB_PARAM_IDs", "Input", params_to_read, IOUtils::nothrow);
        cfg.getValue("GRIB_LEVEL_TYPES", "Input", level_types_to_read, IOUtils::nothrow);
        cfg.getValue("GRIB_LEVELS", "Input", levels_to_read, IOUtils::nothrow);
        if (params_to_read.size() != level_types_to_read.size() || params_to_read.size() != levels_to_read.size()) {
            throw IOException("GRIB_PARAM_IDs, GRIB_LEVEL_TYPES and GRIB_LEVELS must have the same number of elements", AT);
        }
        cfg.getValue("GRIB_ENSEMBLE", "Input", ensemble_to_read, IOUtils::nothrow);
        cfg.getValue("VERBOSE", "Input", verbose, IOUtils::nothrow);
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

    Coords GRIBIO::getGeolocalization(double &cell_x, double &cell_y, const std::map<std::string, double> &grid_params) {
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
        cell_x = CoordsAlgorithms::VincentyDistance(cntr_lat, ll_lon, cntr_lat, ur_lon, bearing) / (double)Ni;
        cell_y = CoordsAlgorithms::VincentyDistance(ll_lat, cntr_lon, ur_lat, cntr_lon, bearing) / (double)Nj;

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

    void GRIBIO::read2Dlevel(CodesHandlePtr &h, Grid2DObject &grid_out) {
        std::vector<double> values;
        std::map<std::string, double> grid_params;
        getGriddedValues(h, values, grid_params);

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

    bool GRIBIO::read2DGrid_indexed(const std::string &in_paramId, const std::string &i_levelType, const long &i_level, const Date i_date, Grid2DObject &grid_out) {
        // TODO: implement the ensemble_number, i.e. remove it cause we dont need it, or choose the right one
        std::vector<CodesHandlePtr> handles = getMessages(file_index, in_paramId, ensemble_to_read, i_levelType);
        for (auto &h : handles) {

            double P1, P2;
            Date base_date = getMessageDateGrib(h, P1, P2, tz_in);

            // see WMO code table5 for definitions of timeRangeIndicator. http://dss.ucar.edu/docs/formats/grib/gribdoc/timer.html
            //  0 -> at base_date + P1
            //  1 -> at base_date
            //  2 -> valid between base_date+P1 and base_date+P2
            //  3 -> average within [base_date+P1 , base_date+P2]
            //  4 -> accumulation from base_date+P1 to base_date+P2
            //  5 -> difference (base_date+P2) - (base_date+P1)
            long timeRange = -999;
            getParameter(h, "timeRangeIndicator", timeRange);

            long level = 0;
            if (i_level != 0)
                getParameter(h, "level", level);

            if (level == i_level && timeRange >= 0) {
                if ((i_date.isUndef()) || (timeRange == 0 && i_date == base_date + P1) || (timeRange == 1 && i_date == base_date) ||
                    ((timeRange == 2 || timeRange == 3) && i_date >= base_date + P1 && i_date <= base_date + P2) || ((timeRange == 4 || timeRange == 5) && i_date == base_date + P2)) {
                    read2Dlevel(h, grid_out);
                    if (timeRange == 3)
                        grid_out *= ((P2 - P1) * 24. * 3600.); // convert avg to sum
                    return true;
                }
            }
        }
        return false;
    }

    // try different IDs that give the same Parameter
    bool GRIBIO::read2DGrid_indexed(const std::vector<std::string> &in_paramId, const std::string &i_levelType, const long &i_level, const Date i_date, Grid2DObject &grid_out) {
        for (const auto &paramId : in_paramId) {
            if (read2DGrid_indexed(paramId, i_levelType, i_level, i_date, grid_out))
                return true;
        }
        return false;
    }

    void GRIBIO::read2DGrid(Grid2DObject &grid_out, const std::string &i_name) {
        const std::string filename(grid2dpath_in + "/" + i_name);
        if (!FileUtils::fileExists(filename))
            throw AccessException(filename, AT); // prevent invalid filenames

        std::vector<CodesHandlePtr> handles = getMessages(filename);

        if (handles.empty())
            throw IOException("No grid found in file \"" + filename + "\"", AT);
        if (handles.size() > 1)
            throw IOException("Multiple grids found in file \"" + filename + "\". Please specify which to load.", AT);

        read2Dlevel(handles.front(), grid_out);
    }

    // TODO: rework from here: make it possible to read a file with multiple dates
    void GRIBIO::read2DGrid(Grid2DObject &grid_out, const MeteoGrids::Parameters &parameter, const Date &date) {
        Date UTC_date = date;
        UTC_date.setTimeZone(tz_in);

        const std::string filename(grid2dpath_in + "/" + grid2d_prefix + UTC_date.toString(Date::NUM).substr(0, 10) + grid2d_ext);
        if (!FileUtils::fileExists(filename)) {
            if (FileUtils::fileExists(grid2dpath_in)) {
                read2DGrid(grid2dpath_in, grid_out, parameter, UTC_date);
            }
        }

        read2DGrid(filename, grid_out, parameter, UTC_date);
    }

    void GRIBIO::readWind(const std::string &filename, const Date &date) {
        if (wind_date == date)
            return; // wind fields are already up to date

        if (read2DGrid_indexed(PARAMETER_MAP.at("wind_speed"), "heightAboveGround", 10, date, VW) ||
            read2DGrid_indexed(PARAMETER_GROUPS.at("wind_speed_10m"), "heightAboveGround", 0, date, VW)) { // TODO: check whether the level indexing and key is correct
            if (!read2DGrid_indexed(PARAMETER_MAP.at("wind_direction"), "heightAboveGround", 10, date, DW) &&
                !read2DGrid_indexed(PARAMETER_GROUPS.at("wind_direction_10m"), "heightAboveGround", 0, date, DW)) // TODO: check whether the level indexing and key is correct
                throw NoDataException("Can not read wind direction in file \"" + filename + "\"", AT);
        } else {
            Grid2DObject U, V;
            if (!read2DGrid_indexed(PARAMETER_MAP.at("wind_u_10m"), "heightAboveGround", 0, date, U)); // TODO: check whether the level indexing and key is correct
            // fallback
            read2DGrid_indexed(PARAMETER_MAP.at("wind_u"), "heightAboveGround", 10, date, U);
            if (!read2DGrid_indexed(PARAMETER_MAP.at("wind_v_10m"), "heightAboveGround", 0, date, V)); // TODO: check whether the level indexing and key is correct
            // fallback
			read2DGrid_indexed(PARAMETER_MAP.at("wind_v"), "heightAboveGround", 10, date, V);

            VW.set(U.getNx(), U.getNy(), U.cellsize, U.llcorner);
            DW.set(U.getNx(), U.getNy(), U.cellsize, U.llcorner);
            for (size_t jj = 0; jj < VW.getNy(); jj++) {
                for (size_t ii = 0; ii < VW.getNx(); ii++) {
                    VW(ii, jj) = sqrt(Optim::pow2(U(ii, jj)) + Optim::pow2(V(ii, jj)));
                    DW(ii, jj) = fmod(IOUtils::UV_TO_DW(U(ii, jj), V(ii, jj)) + bearing_offset, 360.); // turn into degrees [0;360)
                }
            }
        }

        wind_date = date;
    }

    void GRIBIO::read2DGrid(const std::string &filename, Grid2DObject &grid_out, const MeteoGrids::Parameters &parameter,
                            const Date &date) { // Parameters should be read in table 2 if available since this table is the one standardized by WMO
        if (!indexed) {
            file_index = indexFile(filename, paramIdList, ensemble_to_read, verbose);
        }

        // Basic meteo parameters
        if (parameter == MeteoGrids::P)
            read2DGrid_indexed(PARAMETER_MAP.at("pressure"), "surface", 0, date, grid_out);
        if (parameter == MeteoGrids::TA) {
            if (!read2DGrid_indexed(PARAMETER_MAP.at("temperature_2m"), "surface", 0, date, grid_out)) {
                read2DGrid_indexed(PARAMETER_MAP.at("temperature"), "heightAboveGround", 2, date, grid_out);
            }; // TODO: checc whether the level indexing and key is correct
        }
        if (parameter == MeteoGrids::RH) {
            if (!read2DGrid_indexed(PARAMETER_MAP.at("relative_humidity_int"), "heightAboveGround", 2, date, grid_out)) { // TODO: is % or fraction?
                // fallback rel hum %
                if (!read2DGrid_indexed(PARAMETER_MAP.at("relative_humidity_per_cent"), "heightAboveGround", 2, date, grid_out)) {
                    // fallback in case no rel hum is available
                    Grid2DObject ta;
                    if (read2DGrid_indexed(PARAMETER_MAP.at("temperature_2m"), "surface", 0, date, ta) ||
                        read2DGrid_indexed(PARAMETER_MAP.at("temperature"), "heightAboveGround", 2, date, ta)) { // TODO: check whether the level indexing and key is correct
                        if (read2DGrid_indexed(PARAMETER_GROUPS.at("dew_point_temperature"), "heightAboveGround", 2, date, grid_out) ||
                            read2DGrid_indexed(PARAMETER_MAP.at("dew_point_temperature_2m"), "surface", 0, date, grid_out)) {
                            for (size_t jj = 0; jj < grid_out.getNy(); jj++) {
                                for (size_t ii = 0; ii < grid_out.getNx(); ii++) {
                                    grid_out(ii, jj) = Atmosphere::DewPointtoRh(grid_out(ii, jj), ta(ii, jj), true);
                                }
                            }
                        }
                    }
                } else {
                    // need to convert to decimal if we found %
                    grid_out /= 100.;
                }
            }
        }

        if (parameter == MeteoGrids::TSS)
            read2DGrid_indexed(PARAMETER_GROUPS.at("snow_surface_temperature"), "depthBelowLand", 0, date, grid_out); // TODO: check whether the level indexing and key is correct
        if (parameter == MeteoGrids::TSG)
            read2DGrid_indexed(PARAMETER_MAP.at("temperature"), "surface", 0, date, grid_out); // T_G

        // hydrological parameters
        if (parameter == MeteoGrids::PSUM)
            read2DGrid_indexed(PARAMETER_MAP.at("total_precipitation"), "surface", 0, date, grid_out); // tp
        if (parameter == MeteoGrids::ROT)
            read2DGrid_indexed(PARAMETER_MAP.at("runoff"), "depthBelowLandLayer", 0, date, grid_out);
        if (parameter == MeteoGrids::SWE)
            if (!read2DGrid_indexed(PARAMETER_MAP.at("SWE_kgm-2"), "surface", 0, date, grid_out)) {
                read2DGrid_indexed(PARAMETER_MAP.at("SWE_m"), "surface", 0, date, grid_out);
                grid_out *= 1000.; // convert to kg/m^2 TODO: is this the correct conversion?
            }                      // W_SNOW
        if (parameter == MeteoGrids::HS) {
            if (!read2DGrid_indexed(PARAMETER_MAP.at("snow_depth"), "surface", 0, date, grid_out)) {
                Grid2DObject snow_density;
                read2DGrid_indexed(PARAMETER_MAP.at("snow_density"), "surface", 0, date, snow_density); // RHO_SNOW
                read2DGrid_indexed(PARAMETER_MAP.at("SWE_kgm-2"), "surface", 0, date, grid_out);        // W_SNOW
                grid_out.grid2D /= snow_density.grid2D;
            }
        }

        // radiation parameters
        if (parameter == MeteoGrids::ALB) {
            if (read2DGrid_indexed(PARAMETER_GROUPS.at("albedo_per_cent"), "surface", 0, date, grid_out)) {
                grid_out /= 100.;
            } else {
                read2DGrid_indexed(PARAMETER_MAP.at("albedo_backup_int"), "surface", 0, date, grid_out); // ALB_S
            };
        }
        if (parameter == MeteoGrids::ILWR) {
            if (read2DGrid_indexed(PARAMETER_MAP.at("long_wave_radiation_flux"), "surface", 0, date, grid_out)) { // long wave
                grid_out.grid2D *= -1.;
            }
        }
        if (parameter == MeteoGrids::TAU_CLD) { // cloudiness
            if (read2DGrid_indexed(PARAMETER_MAP.at("medium_cloud_cover_per_cent"), "surface", 0, date, grid_out))
                grid_out.grid2D /= 100.;
            else
                read2DGrid_indexed(PARAMETER_MAP.at("cloud_cover_per"), "surface", 0, date, grid_out);
        }

        if (parameter == MeteoGrids::ISWR) {
            if (read2DGrid_indexed(PARAMETER_GROUPS.at("short_wave_radiation_flux"), "surface", 0, date, grid_out)) { // short wave <<< why is the -1 there?
                grid_out.grid2D *= -1.;
            } else if (!read2DGrid_indexed(PARAMETER_MAP.at("global_radiation_flux"), "surface", 0, date, grid_out)) { // GLOB
                Grid2DObject diff;
                read2DGrid_indexed(PARAMETER_MAP.at("short_wave_radiation_flux_diffuse"), "surface", 0, date, diff);    // diffuse rad, ASWDIFD_S << 23.201 -- change
                read2DGrid_indexed(PARAMETER_MAP.at("short_wave_radiation_flux_direct"), "surface", 0, date, grid_out); // direct rad, ASWDIR_S<< 22.201 -- change
                grid_out.grid2D += diff.grid2D;
            }
        }

        // DEM parameters
        if (parameter == MeteoGrids::DEM)
            read2DGrid_indexed(PARAMETER_GROUPS.at("DEM"), "surface", 0, date, grid_out); // HSURF
        if (parameter == MeteoGrids::SLOPE) {
            read2DGrid_indexed(PARAMETER_MAP.at("subgrid_slope"), "surface", 0, date, grid_out); // SLO_ANG
            grid_out.grid2D *= Cst::to_deg;
        }
        if (parameter == MeteoGrids::AZI) {
            read2DGrid_indexed(PARAMETER_MAP.at("subgrid_angle_eastward"), "surface", 0, date, grid_out); // TODO: check if it is the right conversion
            for (size_t jj = 0; jj < grid_out.getNy(); jj++) {
                for (size_t ii = 0; ii < grid_out.getNx(); ii++) {
                    grid_out(ii, jj) = fmod(grid_out(ii, jj) * Cst::to_deg + 360. + bearing_offset, 360.); // turn into degrees [0;360)
                }
            }
        }

        // Wind parameters
        if (parameter == MeteoGrids::VW_MAX)
            read2DGrid_indexed(PARAMETER_MAP.at("maximum_wind_velocity"), "heightAboveGround", 10, date, grid_out); // VMAX_10M 10m
        if (parameter == MeteoGrids::W)
            if (!read2DGrid_indexed(PARAMETER_GROUPS.at("vertical_velocity"), "hybrid", 10, date, grid_out)) {
                read2DGrid_indexed(PARAMETER_MAP.at("vertical_velocity"), "heightAboveGround", 10, date, grid_out); // W
            };
        // we need to use VW, DW, correct for re-projection and recompute U,V
        if (parameter == MeteoGrids::U) {
            readWind(filename, date);
            for (size_t jj = 0; jj < grid_out.getNy(); jj++) {
                for (size_t ii = 0; ii < grid_out.getNx(); ii++) {
                    grid_out(ii, jj) = IOUtils::VWDW_TO_U(VW(ii, jj), DW(ii, jj));
                }
            }
        }
        if (parameter == MeteoGrids::V) {
            readWind(filename, date);
            for (size_t jj = 0; jj < grid_out.getNy(); jj++) {
                for (size_t ii = 0; ii < grid_out.getNx(); ii++) {
                    grid_out(ii, jj) = IOUtils::VWDW_TO_V(VW(ii, jj), DW(ii, jj));
                }
            }
        }
        if (parameter == MeteoGrids::DW) {
            readWind(filename, date);
            grid_out = DW;
        }
        if (parameter == MeteoGrids::VW) {
            readWind(filename, date);
            grid_out = VW;
        }

        if (grid_out.empty()) {
            ostringstream ss;
            ss << "No suitable data found for parameter " << MeteoGrids::getParameterName(parameter) << " ";
            ss << "at time step " << date.toString(Date::ISO) << " in file \"" << filename << "\"";
            throw NoDataException(ss.str(), AT);
        }

        // correcting wind speeds
        /*if (parameter==MeteoGrids::U || parameter==MeteoGrids::V || parameter==MeteoGrids::W || parameter==MeteoGrids::VW || parameter==MeteoGrids::VW_MAX) {
            //we need to compute the wind at 7.5m
            Grid2DObject Z0;
            if (read2DGrid_indexed(83.2, 1, 0, date, Z0)) { //Z0
                for (size_t jj=0; jj<grid_out.getNy(); jj++) {
                    for (size_t ii=0; ii<grid_out.getNx(); ii++) {
                        grid_out(ii,jj) = Atmosphere::windLogProfile(grid_out(ii,jj), 10., 7.5, Z0(ii,jj));
                    }
                }
            } else {
                const double wind_factor = Atmosphere::windLogProfile(1., 10., 7.5, 0.03);
                grid_out.grid2D *= wind_factor;
            }
        }*/
    }

    void GRIBIO::readDEM(DEMObject &dem_out) {
        const Date d; // ie: undef. This will be caught when reading the GRIB file
        const std::string filename = cfg.get("DEMFILE", "Input");
        read2DGrid(filename, dem_out, MeteoGrids::DEM, d);
        if (update_dem) {
            dem_out.update();
        } else {
            const int dem_ppt = dem_out.getUpdatePpt();
            if (dem_ppt & DEMObject::SLOPE) {
                Grid2DObject slope;
                read2DGrid(filename, slope, MeteoGrids::SLOPE, d);
                dem_out.slope = slope.grid2D;
                Grid2DObject azi;
                read2DGrid(filename, azi, MeteoGrids::AZI, d);
                dem_out.azi = azi.grid2D;
            }
            if (dem_ppt & DEMObject::NORMAL || dem_ppt & DEMObject::CURVATURE) {
                // we will only update the normals and/or curvatures, then revert update properties
                if (dem_ppt & DEMObject::NORMAL && dem_ppt & DEMObject::CURVATURE)
                    dem_out.setUpdatePpt((DEMObject::update_type)(DEMObject::NORMAL | DEMObject::CURVATURE));
                else if (dem_ppt & DEMObject::NORMAL)
                    dem_out.setUpdatePpt(DEMObject::NORMAL);
                else if (dem_ppt & DEMObject::CURVATURE)
                    dem_out.setUpdatePpt(DEMObject::CURVATURE);

                dem_out.update();
                dem_out.setUpdatePpt((DEMObject::update_type)dem_ppt);
            }

            dem_out.updateAllMinMax();
        }
    }

    void GRIBIO::scanMeteoPath() {
        std::list<std::string> dirlist;
        FileUtils::readDirectory(meteopath_in, dirlist, meteo_ext);
        dirlist.sort();

        // Check date in every filename and cache it
        std::list<std::string>::const_iterator it = dirlist.begin();
        while ((it != dirlist.end())) {
            const std::string &filename = *it;
            const std::string::size_type spos = filename.find_first_of("0123456789");
            if (spos != std::string::npos) {
                Date date;
                IOUtils::convertString(date, filename.substr(spos, 10), tz_in);
                const std::pair<Date, std::string> tmp(date, filename);

                cache_meteo_files.push_back(tmp);
                ++it;
            } else {
                if (!warned_date_in_file) {
                    std::cerr << "Warning: no date found in meteo file names. Assuming multiple dates in file instead.\n";
                    warned_date_in_file = true;
                }
                cache_meteo_files.push_back(std::pair<Date, std::string>(Date(), filename));
            }
        }
    }

    void GRIBIO::readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecMeteo) {
        if (!meteo_initialized) {
            readStations(vecPts);
            scanMeteoPath();
            meteo_initialized = true;
        }

        vecMeteo.clear();

        std::vector<double> lats(vecPts.size());
        std::vector<double> lons(vecPts.size());

        std::vector<StationData> meta; // metadata for meteo time series
        bool meta_ok = false;          // set to true once the metadata have been read

        // find index of first time step
        size_t idx_start = 0;
        bool start_found = false;

        if (!warned_date_in_file) {
            for (idx_start = 0; idx_start < cache_meteo_files.size(); idx_start++) {
                if (dateStart < cache_meteo_files[idx_start].first) {
                    start_found = true;
                    break;
                }
            }

            if (start_found == false) {
                return;
            }
        }

        if (idx_start > 0)
            idx_start--; // start with first element before dateStart (useful for resampling)

        try {
            for (size_t ii = idx_start; ii < cache_meteo_files.size(); ii++) {
                const Date &date = cache_meteo_files[ii].first;
                if (date > dateEnd)
                    break;
                const std::string filename(meteopath_in + "/" + cache_meteo_files[ii].second);

                if (!indexed) {
                    indexFile(filename, paramIdList, ensemble_to_read, verbose); // this will also read geolocalization
                }
                if (!meta_ok) {
                    if (readMeteoMeta(vecPts, meta, lats, lons) == false) {
                        // some points have been removed vecPts has been changed -> re-reading
                        lats.clear();
                        lons.clear();
                        readMeteoMeta(vecPts, meta, lats, lons);
                    }
                    vecMeteo.insert(vecMeteo.begin(), vecPts.size(), std::vector<MeteoData>()); // allocation for the vectors now that we know how many true stations we have
                    meta_ok = true;
                }

                std::vector<MeteoData> Meteo;
                readMeteoStep(meta, lats, lons, date, Meteo);
                for (size_t jj = 0; jj < vecPts.size(); jj++)
                    vecMeteo[jj].push_back(Meteo[jj]);
            }
        } catch (...) {
            throw;
        }
    }

    bool GRIBIO::removeDuplicatePoints(std::vector<Coords> &vecPoints, std::vector<double> &lats, std::vector<double> &lons) { // remove potential duplicates. Returns true if some have been removed
        const size_t npoints = vecPoints.size();
        std::vector<size_t> deletions;
        deletions.reserve(npoints);
        for (size_t ii = 0; ii < npoints; ii++) {
            const double lat = lats[ii];
            const double lon = lons[ii];
            for (size_t jj = ii + 1; jj < npoints; jj++) {
                if (lat == lats[jj] && lon == lons[jj]) {
                    deletions.push_back(jj);
                }
            }
        }

        // we need to erase from the end in order to keep the index unchanged...
        for (size_t ii = deletions.size(); ii-- > 0;) {
            const size_t index = deletions[ii - 1];
            vecPoints.erase(vecPoints.begin() + index);
        }

        if (!deletions.empty())
            return true;
        return false;
    }

    bool GRIBIO::readMeteoMeta(
        std::vector<Coords> &vecPoints, std::vector<StationData> &stations, std::vector<double> &lats,
        std::vector<double> &lons) { // return true if the metadata have been read, false if it needs to be re-read (ie: some points were leading to duplicates -> vecPoints has been changed)
        stations.clear();

        std::vector<CodesHandlePtr> handles = getMessages(file_index, PARAMETER_GROUPS.at("DEM"), ensemble_to_read, "surface");

        if (handles.empty()) {
            throw IOException("No DEM grid found in GRIB file!", AT);
        }

        const size_t npoints = vecPoints.size();

        for (auto &h : handles) {

            if (h == nullptr) {
                throw IOException("Can not find DEM grid in GRIB file!", AT);
            }

            std::map<std::string, double> grid_params = getGridParameters(h);
            latitudeOfNorthernPole = grid_params.at("latitudeOfNorthernPole");
            longitudeOfNorthernPole = grid_params.at("longitudeOfNorthernPole");
            long Ni = static_cast<long>(grid_params.at("Ni"));

            // build GRIB local coordinates for the points
            for (size_t ii = 0; ii < npoints; ii++) {
                CoordsAlgorithms::trueLatLonToRotated(latitudeOfNorthernPole, longitudeOfNorthernPole, vecPoints[ii].getLat(), vecPoints[ii].getLon(), lats[ii], lons[ii]);
            }

            // retrieve nearest points
            std::vector<double> outlats_vec(npoints), outlons_vec(npoints), altitudes_vec(npoints), distances_vec(npoints);
            std::vector<int> indexes_vec(npoints);
            std::tie(outlats_vec, outlons_vec, distances_vec, altitudes_vec, indexes_vec) = getNearestValues_grib(h, lats, lons);

            // remove potential duplicates
            if (removeDuplicatePoints(vecPoints, outlats_vec, outlons_vec) == true)
                return false;

            // fill metadata
            for (size_t ii = 0; ii < npoints; ii++) {
                StationData sd;
                sd.position.setProj(coordin, coordinparam);
                double true_lat, true_lon;
                CoordsAlgorithms::rotatedToTrueLatLon(latitudeOfNorthernPole, longitudeOfNorthernPole, outlats_vec[ii], outlons_vec[ii], true_lat, true_lon);
                sd.position.setLatLon(true_lat, true_lon, altitudes_vec[ii]);
                sd.stationID = "Point_" + IOUtils::toString(indexes_vec[ii]);
                ostringstream ss2;
                ss2 << "GRIB point (" << indexes_vec[ii] % Ni << "," << indexes_vec[ii] / Ni << ")";
                sd.stationName = ss2.str();
                stations.push_back(sd);
            }
        }

        return true;
    }

    bool GRIBIO::readMeteoValues(const std::string &paramId, const std::string &levelType, const long &i_level, const Date &i_date, const size_t &npoints, std::vector<double> &lats,
                                 std::vector<double> &lons, std::vector<double> &values) {
        std::vector<CodesHandlePtr> handles = getMessages(file_index, paramId, ensemble_to_read, levelType);

        if (handles.empty()) {
            throw IOException("No grid found in GRIB file!", AT);
        }

        for (auto &h : handles) {
            if (!h) {
                throw IOException("Unable to create grib handle from index", AT);
            }

            double P1, P2;
            Date base_date = getMessageDateGrib(h, P1, P2, tz_in);

            // see WMO code table5 for definitions of timeRangeIndicator. http://dss.ucar.edu/docs/formats/grib/gribdoc/timer.html
            long timeRange = -999;
            getParameter(h, "timeRangeIndicator", timeRange);

            long level = 0;
            if (i_level != 0)
                getParameter(h, "level", level);
            if (level == i_level && timeRange >= 0) {
                if ((i_date.isUndef()) || (timeRange == 0 && i_date == base_date + P1) || (timeRange == 1 && i_date == base_date) ||
                    ((timeRange == 2 || timeRange == 3) && i_date >= base_date + P1 && i_date <= base_date + P2) || ((timeRange == 4 || timeRange == 5) && i_date == base_date + P2)) {
                    std::vector<double> outlats(npoints);
                    std::vector<double> outlons(npoints);
                    std::vector<double> distances(npoints);
                    std::vector<int> indexes(npoints);
                    std::tie(outlats, outlons, distances, values, indexes) = getNearestValues_grib(h, lats, lons);
                    return true;
                }
            }
        }
        return false;
    }

    // try different IDs that give the same parameter
    bool GRIBIO::readMeteoValues(const std::vector<std::string> &paramId_list, const std::string &levelType, const long &i_level, const Date &i_date, const size_t &npoints, std::vector<double> &lats,
                                 std::vector<double> &lons, std::vector<double> &values) {
        for (const auto &paramId : paramId_list) {
            if (readMeteoValues(paramId, levelType, i_level, i_date, npoints, lats, lons, values))
                return true;
        }
        return false;
    }

    void GRIBIO::fillMeteo(std::vector<double> &values, const MeteoData::Parameters &param, const size_t &npoints, std::vector<MeteoData> &Meteo) {
        for (size_t ii = 0; ii < npoints; ii++) {
            Meteo[ii](param) = values[ii];
        }
    }

    void GRIBIO::readMeteoStep(std::vector<StationData> &stations, std::vector<double> &lats, std::vector<double> &lons, const Date i_date, std::vector<MeteoData> &Meteo) {
        const size_t npoints = stations.size();

        for (size_t ii = 0; ii < npoints; ii++) {
            MeteoData md;
            md.meta = stations[ii];
            md.date = i_date;
            Meteo.push_back(md);
        }

        std::vector<double> values(npoints), values2(npoints); //

        // basic meteorological parameters
        if (readMeteoValues(PARAMETER_MAP.at("pressure"), "surface", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::P, npoints, Meteo); // PS
        if (readMeteoValues(PARAMETER_MAP.at("temperature_2m"), "surface", 0, i_date, npoints, lats, lons, values) ||
            readMeteoValues(PARAMETER_MAP.at("temperature"), "heightAboveGround", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::TA, npoints, Meteo); // TODO: check whether the level indexing and key is correct
        if (readMeteoValues(PARAMETER_GROUPS.at("snow_surface_temperature"), "depthBelowLand", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::TSS, npoints, Meteo);
        if (readMeteoValues(PARAMETER_MAP.at("temperature"), "surface", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::TSG, npoints, Meteo);
        if (readMeteoValues(PARAMETER_MAP.at("relative_humidity_per_cent"), "heightAboveGround", 2, i_date, npoints, lats, lons, values)) {
            for (size_t ii = 0; ii < values.size(); ii++) {
                values[ii] /= 100.; // TODO: check if this is the right conversion
            }
            fillMeteo(values, MeteoData::RH, npoints, Meteo);
        } else if (readMeteoValues(PARAMETER_MAP.at("relative_humidity_int"), "heightAboveGround", 2, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::RH, npoints, Meteo);
        else if (readMeteoValues(PARAMETER_GROUPS.at("dew_point_temperature"), "heightAboveGround", 2, i_date, npoints, lats, lons, values) ||
                 readMeteoValues(PARAMETER_MAP.at("dew_point_temperature_2m"), "surface", 0, i_date, npoints, lats, lons, values)) { // TODO: check indexing of levels
            for (size_t ii = 0; ii < npoints; ii++) {
                if (Meteo[ii](MeteoData::TA) != IOUtils::nodata)
                    Meteo[ii](MeteoData::RH) = Atmosphere::DewPointtoRh(values[ii], Meteo[ii](MeteoData::TA), true);
            }
        }

        // hydrological parameters
        if (readMeteoValues(PARAMETER_MAP.at("total_precipitation"), "surface", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::PSUM, npoints, Meteo); // tp
        if (readMeteoValues(PARAMETER_MAP.at("snow_depth"), "surface", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::HS, npoints, Meteo);
		else if (readMeteoValues(PARAMETER_MAP.at("snow_density"), "surface", 0, i_date, npoints, lats, lons, values)) {
			bool readSWE = readMeteoValues(PARAMETER_MAP.at("SWE_kgm-2"), "surface", 0, i_date, npoints, lats, lons, values2);
			if (!readSWE) { 
				readSWE = readMeteoValues(PARAMETER_MAP.at("SWE_m"), "surface", 0, i_date, npoints, lats, lons, values2);
				if (readSWE) {
					for (size_t ii = 0; ii < values2.size(); ii++) {
						values2[ii] *= 1000.; // convert to kg/m^2 TODO: is this the correct conversion?
					}
				}
			}
			if (readSWE) {
				for (size_t ii = 0; ii < npoints; ii++) {
					Meteo[ii](MeteoData::HS) = values2[ii] / values[ii];
				}
			}
		}

        // radiation parameters
        if (readMeteoValues(PARAMETER_GROUPS.at("short_wave_radiation_flux"), "surface", 0, i_date, npoints, lats, lons, values)) { // long wave
            for (size_t ii = 0; ii < npoints; ii++) {
                Meteo[ii](MeteoData::ISWR) = -values[ii];
            }
        } else if (readMeteoValues(PARAMETER_MAP.at("surface_long_wave_radiation_downwards"), "surface", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::ILWR, npoints, Meteo); // ALWD_S <<< 25.201
        if (readMeteoValues(PARAMETER_MAP.at("global_radiation_flux"), "surface", 0, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::ISWR, npoints, Meteo); // GLOB_H
        else {
            const bool read_dir = (readMeteoValues(PARAMETER_MAP.at("long_wave_radiation_flux"), "surface", 0, i_date, npoints, lats, lons, values)); // O_ASWDIR_S
            const bool read_diff = (readMeteoValues(PARAMETER_MAP.at("global_radiation_flux"), "surface", 0, i_date, npoints, lats, lons, values2)     // O_ASWDIFD_S <<< was 117.2
                                    || readMeteoValues(PARAMETER_MAP.at("surface_solar_radiation_difference"), "surface", 0, i_date, npoints, lats, lons, values2)); // ASWDIFD_S << was 23.201
            if (read_dir && read_diff) {
                for (size_t ii = 0; ii < npoints; ii++) {
                    Meteo[ii](MeteoData::ISWR) = values[ii] + values2[ii];
                }
            }
        }
		bool albedo_int = false;
		bool albedo_per_cent = readMeteoValues(PARAMETER_GROUPS.at("albedo_per_cent"), "surface", 0, i_date, npoints, lats, lons, values);
        if (!albedo_per_cent) { // ALB_RAD
			albedo_int = readMeteoValues(PARAMETER_MAP.at("albedo_backup_int"), "surface", 0, i_date, npoints, lats, lons, values);
        } 
		if (albedo_int || albedo_per_cent) {
			double frac = (albedo_per_cent) ? 100. : 1.;
			for (size_t ii = 0; ii < npoints; ii++) {
                if (Meteo[ii](MeteoData::ISWR) != IOUtils::nodata)
                    Meteo[ii](MeteoData::RSWR) = Meteo[ii](MeteoData::ISWR) * values[ii] / frac;
            }
		}

        // Wind parameters
        if (readMeteoValues(PARAMETER_MAP.at("maximum_wind_velocity"), "heightAboveGround", 10, i_date, npoints, lats, lons, values))
            fillMeteo(values, MeteoData::VW_MAX, npoints, Meteo); // VMAX_10M << was 187.201
        if (readMeteoValues(PARAMETER_GROUPS.at("wind_direction"), "heightAboveGround", 10, i_date, npoints, lats, lons, values) || readMeteoValues(PARAMETER_GROUPS.at("wind_direction_10m"), "surface", 0, i_date, npoints, lats, lons, values)) // TODO: check whether the level indexing and key is correct
            fillMeteo(values, MeteoData::DW, npoints, Meteo); // DD_10M
        else {
            if ((readMeteoValues(PARAMETER_MAP.at("wind_v"), "heightAboveGround", 10, i_date, npoints, lats, lons, values) || readMeteoValues(PARAMETER_MAP.at("wind_v_10m"), "surface", 0, i_date, npoints, lats, lons, values))        // V_10M
                && (readMeteoValues(PARAMETER_MAP.at("wind_u"), "heightAboveGround", 10, i_date, npoints, lats, lons, values2) || readMeteoValues(PARAMETER_MAP.at("wind_u_10m"),"surface", 0, i_date, npoints, lats, lons, values2))) { // U_10M
                for (size_t ii = 0; ii < npoints; ii++) {
                    Meteo[ii](MeteoData::DW) = fmod(IOUtils::UV_TO_DW(values2[ii], values[ii]) + bearing_offset, 360.); // turn into degrees [0;360)
                }
            }
        }
        if (readMeteoValues(PARAMETER_MAP.at("wind_speed"), "heightAboveGround", 10, i_date, npoints, lats, lons, values) || readMeteoValues(PARAMETER_GROUPS.at("wind_speed_10m"), "surface", 0, i_date, npoints, lats, lons, values)) // TODO: check whether the level indexing and key is correct
            fillMeteo(values, MeteoData::VW, npoints, Meteo); // FF_10M
        else {
            if ((readMeteoValues(PARAMETER_MAP.at("wind_v"), "heightAboveGround", 10, i_date, npoints, lats, lons, values) || readMeteoValues(PARAMETER_MAP.at("wind_v_10m"), "surface", 0, i_date, npoints, lats, lons, values))        // V_10M
                && (readMeteoValues(PARAMETER_MAP.at("wind_u"), "heightAboveGround", 10, i_date, npoints, lats, lons, values2) || readMeteoValues(PARAMETER_MAP.at("wind_u_10m"),"surface", 0, i_date, npoints, lats, lons, values2))) { // U_10M
                for (size_t ii = 0; ii < npoints; ii++) {
                    Meteo[ii](MeteoData::VW) = sqrt(Optim::pow2(values[ii]) + Optim::pow2(values2[ii]));
                }
            }
        }
    }

} // namespace
