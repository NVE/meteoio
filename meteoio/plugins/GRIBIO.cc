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
#include <meteoio/FileUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/plugins/GRIBIO.h>
#include <algorithm>


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
     * @todo When add parameter is available for MeteoGrids, support adding unknown parameters via the parameter table for GRIB
     */

    static const double plugin_nodata = IOUtils::nodata; // plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
    static const std::string default_extension = ".grib";
    const std::string GRIBIO::default_table = "doc/resources/GRIB_param.tbl";

    static size_t findDate(const std::vector<GRIBFile> &cache, const Date &date) {
        auto it = std::find_if(cache.begin(), cache.end(), [&date](const GRIBFile &file) { return file.isValidDate(date); });

        return (it != cache.end()) ? std::distance(cache.begin(), it) : IOUtils::npos;
    }

    static void handleConversions(Grid2DObject& grid_out, const double& paramId) { 
        if (paramId == 163) { // slope
            grid_out.grid2D *= 90.;
        } else if (paramId == 162) { // azi
            grid_out.grid2D *= Cst::to_deg;
            grid_out.grid2D -= 90.;
        } else if (paramId == 129) { // geopotential
            grid_out.grid2D /= Cst::gravity;
        } 
    }

    // ----------------------------- INITIALIZE -----------------------------
    GRIBIO::GRIBIO(const std::string &configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(), meteopath_in(), grid2dpath_in(), table_path(), meteo_ext(),
                                                    meteo_pattern(), grid2d_ext(), grid_2d_pattern(), recursive_search(false), verbose(false), update_dem(false), 
                                                    bearing_offset(IOUtils::nodata), latitudeOfNorthernPole(), longitudeOfNorthernPole(), llcorner_initialized(false), llcorner(), 
                                                    cellsize(), factor_x(), factor_y(), grid_initialized(false), meteo_initialized(false), parameter_table(), 
                                                    cache_meteo(), cache_grid2d(), vecPts() {
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
        initialize();
    }

    GRIBIO::GRIBIO(const Config &cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(), meteopath_in(), grid2dpath_in(), table_path(), meteo_ext(),
                                                    meteo_pattern(), grid2d_ext(), grid_2d_pattern(), recursive_search(false), verbose(false), update_dem(false), 
                                                    bearing_offset(IOUtils::nodata), latitudeOfNorthernPole(), longitudeOfNorthernPole(), llcorner_initialized(false), llcorner(), 
                                                    cellsize(), factor_x(), factor_y(), grid_initialized(false), meteo_initialized(false), parameter_table(), 
                                                    cache_meteo(), cache_grid2d(), vecPts() {
        IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
        initialize();
    }

    void GRIBIO::initialize() {
        setOptions();
        initTable();
    }

    void GRIBIO::setOptions() {
        const std::string in_2d = IOUtils::strToUpper(cfg.get("GRID2D", "Input", ""));
        const std::string in_meteo = IOUtils::strToUpper(cfg.get("METEO", "Input", ""));
        const std::string in_dem = IOUtils::strToUpper(cfg.get("DEM", "Input", ""));

        bool meteo_and_grid = in_2d == "GRIB" && in_meteo == "GRIB";

        if (in_2d == "GRIB" || in_meteo == "GRIB" || in_dem == "GRIB") {
            cfg.getValue("GRIB_TABLE", "Input", table_path, IOUtils::nothrow);
            cfg.getValue("VERBOSE", "Input", verbose, IOUtils::nothrow);
            cfg.getValue("RECURSIVE", "Input", recursive_search, IOUtils::nothrow);
        }

        if (in_2d == "GRIB") { // keep it synchronized with IOHandler.cc for plugin mapping!!
            cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
            cfg.getValue("GRID2DEXT", "Input", grid2d_ext, IOUtils::nothrow);
            if (grid2d_ext == "none")
                grid2d_ext.clear();
            cfg.getValue("GRID2DPATTERN", "Input", grid_2d_pattern, IOUtils::nothrow);
            if (grid_2d_pattern == "none")
                grid_2d_pattern.clear();
        }

        if (in_meteo == "GRIB") {
            cfg.getValue("METEOEXT", "Input", meteo_ext, IOUtils::nothrow);
            if (meteo_ext == "none")
                meteo_ext.clear();
            cfg.getValue("METEOPATTERN", "Input", meteo_pattern, IOUtils::nothrow);
            if (meteo_pattern == "none")
                meteo_pattern.clear();
        }

        if (in_dem == "GRIB") {
            cfg.getValue("GRIB_DEM_UPDATE", "Input", update_dem, IOUtils::nothrow);
        }

        if (meteo_and_grid && meteo_ext == grid2d_ext && meteo_pattern == grid_2d_pattern && meteopath_in == grid2dpath_in)
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
#ifdef DEBUG
        parameter_table.printTable();
#endif
    }

    void GRIBIO::scanPath(const std::string &in_path, const std::string &in_ext, const std::string &in_pattern, std::vector<GRIBFile> &cache) {
        std::list<std::string> dirlist;
        FileUtils::readDirectory(in_path, dirlist, in_ext, recursive_search);

        std::list<std::string>::const_iterator it = dirlist.begin();
        while ((it != dirlist.end())) {
            const std::string filename = *it;
            if (filename.find(in_pattern) != std::string::npos) {
                const std::string fullpath = in_path + "/" + filename;
                cache.push_back(GRIBFile(fullpath, parameter_table.getIndexes())); // Files will not be read twice, as we enforce a seperation between meteo and grid
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

            std::cout << "\tRead virtual station " << vecPoints.back().toString(Coords::LATLON) << "\n";
        }
    }

    Coords GRIBIO::getGeolocalization(double &cellsize_x, double &cellsize_y, const std::map<std::string, double> &grid_params) {
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
    void GRIBIO::read2Dlevel(CodesHandlePtr &h, Grid2DObject &grid_out, const std::map<std::string, double> &grid_params) {
        setMissingValue(h, plugin_nodata);
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

        double paramId;
        getParameter(h, "paramId", paramId);
        handleConversions(grid_out, paramId);
        if (verbose) {
            std::cout << "Read " << values.size() << " values from GRIB file" << std::endl;
            std::cout << "Parameter " << paramId << std::endl;
        }
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

        read2Dlevel(handles.front(), grid_out, getGridParameters(handles.front()));
    }

    static std::vector<CodesHandlePtr> findMessages(GRIBFile &file, const std::string &param_key, const std::string &level_key, const std::string &level_type, const std::string &paramID_string,
                                                    double paramID_double, long paramID_long, const bool& verbose, const std::string& param_name) {
        double npos_double = static_cast<double>(IOUtils::npos);
        long npos_long = static_cast<long>(IOUtils::npos);

        if (!paramID_string.empty() && paramID_double == npos_double && paramID_long == npos_long) {
            return file.listParameterMessages(param_key, paramID_string, level_key, level_type);
        } else if (paramID_string.empty() && paramID_double != npos_double && paramID_long == npos_long) {
            return file.listParameterMessages(param_key, paramID_double, level_key, level_type);
        } else if (paramID_string.empty() && paramID_double == npos_double && paramID_long != npos_long) {
            return file.listParameterMessages(param_key, paramID_long, level_key, level_type);
        } else {
            if (verbose) {
                std::cout << "No parameter id provided for "+param_name << std::endl;
            }
            return {};
        }   
    }


    static std::vector<CodesHandlePtr> extractParameterInfoAndFindMessages(GRIBFile &file, const std::string &param_name, const GRIBTable &parameter_table, long &level_no, const bool& verbose) {
        // get the paramID
        std::string paramID_string;
        double paramID_double;
        long paramID_long;
        parameter_table.getParamId(param_name, paramID_string, paramID_double, paramID_long);

        // get additional information from the parameter table
        std::string param_key = parameter_table.getParamKey();
        std::string level_key = parameter_table.getLevelKey();
        level_no = parameter_table.getLevelNo(param_name);
        std::string level_type = parameter_table.getLevelType(param_name);
        // check type of level and where it is lost
        return findMessages(file, param_key, level_key, level_type, paramID_string, paramID_double, paramID_long, verbose, param_name);
    }
    
    static std::vector<CodesHandlePtr> extractParameterInfoAndFindMessages(GRIBFile &file, const MeteoGrids::Parameters &parameter, const GRIBTable &parameter_table, long &level_no, const bool& verbose) {
        std::string param_name = MeteoGrids::getParameterName(parameter);
        return extractParameterInfoAndFindMessages(file, param_name, parameter_table, level_no, verbose);
    };


    void GRIBIO::read2DGrid(Grid2DObject &grid_out, const MeteoGrids::Parameters &parameter, const Date &date) {
        if (!grid_initialized) {
            scanPath(grid2dpath_in, grid2d_ext, grid_2d_pattern, cache_grid2d);
            grid_initialized = true;
        }

        size_t idx = findDate(cache_grid2d, date);
        if (idx == IOUtils::npos) {
            if (verbose) {
                std::cout << "No grid found for the specified date" << std::endl;
            }
            return;
        }
        if (verbose) {
            std::cout << "Reading grid for date " << date.toString() << std::endl;
        }

        long level_no;
        std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(cache_grid2d[idx], parameter, parameter_table, level_no, verbose);

        if (messages.empty()) {
            if (verbose) {
                std::cout << "No messages found for the specified parameter "+MeteoGrids::getParameterName(parameter) << std::endl;
            }
            return;
        }

        for (auto &m : messages) {
            long level = 0;
            getParameter(m, "level", level);
            // if multiple timepoints in a grib file:
            if (getMessageDateGrib(m, 0) != date) {
                if (verbose)
                    std::cout << "No messages found for the specified date "+ getMessageDateGrib(m,0).toString(Date::ISO) << std::endl;
                continue;
            }
            if (level == level_no) {
                read2Dlevel(m, grid_out, cache_grid2d[idx].getGridParams());
                return;
            }
        }
        if (verbose)
            std::cout << "No messages found for the specified level" << std::endl;
        return;
    };

    // ---------------------------- DIGITAL ELEVATION MODEL -----------------------------
    void GRIBIO::processSingleMessage(Grid2DObject &dem_out, GRIBFile &dem_file, const GRIBTable &dem_table, const MeteoGrids::Parameters &parameter) {
        long level_no;
        std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(dem_file, parameter, dem_table, level_no, verbose);

        if (messages.empty()) {
            throw IOException("No messages containing DEM information found." AT);
        }
        if (messages.size() > 1) {
            throw IOException("Multiple messages containing DEM information found in file" + dem_file.getFilename(), AT);
        }

        read2Dlevel(messages.front(), dem_out, dem_file.getGridParams());
    }

    void GRIBIO::readDEM(DEMObject &dem_out) {
        const std::string filename = cfg.get("DEMFILE", "Input");

        GRIBTable dem_table(table_path);
        GRIBFile dem_file(filename, dem_table.getIndexes());

        processSingleMessage(dem_out, dem_file, dem_table, MeteoGrids::DEM);

        if (update_dem) {
            dem_out.update();
        } else {
            const int dem_ppt = dem_out.getUpdatePpt();
            if (dem_ppt & DEMObject::SLOPE) {
                Grid2DObject slope;
                processSingleMessage(slope, dem_file, dem_table, MeteoGrids::SLOPE);
                dem_out.slope = slope.grid2D;
                Grid2DObject azi;
                processSingleMessage(azi, dem_file, dem_table, MeteoGrids::AZI);
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


    // ---------------------------- METEO DATA -----------------------------
    // ---------------------------- STATIC HELPERS
    static bool compareByDate(const GRIBFile &a, const GRIBFile &b) { return a.getStartDate() < b.getStartDate(); }
    static bool compareToDate(const GRIBFile &a, const Date &b) { return a.getStartDate() < b; } // TODO: check if it is < or >

    static std::vector<MeteoData> createMeteoDataVector(const std::vector<StationData>& stations, const Date& date) {
        std::vector<MeteoData> vecMeteo;
        const size_t npoints = stations.size();

        for (size_t ii = 0; ii < npoints; ii++) {
            MeteoData md;
            md.meta = stations[ii];
            md.date = date;
            vecMeteo.push_back(md);
        }

        return vecMeteo;
    }

    static void processMessages(std::vector<CodesHandlePtr>& messages, const long& level_no, std::vector<MeteoData>& vecMeteo, const std::vector<double>& lats, const std::vector<double>& lons, const size_t& npoints, const size_t& par_index) {
        Date curr_date = vecMeteo.front().date;
        for (auto &m : messages) {
            if (getMessageDateGrib(m, 0) != curr_date)
                continue;
            long level = 0;
            getParameter(m, "level", level);
            if (level == level_no) {
                setMissingValue(m, plugin_nodata);
                std::vector<double> outlats_vec(npoints), outlons_vec(npoints), distances_vec(npoints), values(npoints);
                std::vector<int> indexes_vec(npoints);
                getNearestValues_grib(m, lats, lons, outlats_vec, outlons_vec, distances_vec, values, indexes_vec);

                for (size_t ii = 0; ii < npoints; ii++) {
                    vecMeteo[ii](par_index) = values[ii];
                }
            }
        }
    }



    // ---------------------------- METEO READING -----------------------------
    void GRIBIO::readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecvecMeteo) {
        if (!meteo_initialized) {
            readStations(vecPts);
            scanPath(meteopath_in, meteo_ext, meteo_pattern, cache_meteo);
            meteo_initialized = true;
            std::sort(cache_meteo.begin(), cache_meteo.end(), compareByDate);
        }

        vecvecMeteo.clear();

        std::vector<double> lats(vecPts.size());
        std::vector<double> lons(vecPts.size());

        std::vector<StationData> stations;
        bool meta_ok = false;

        auto it = std::lower_bound(cache_meteo.begin(), cache_meteo.end(), dateStart, compareToDate);

        size_t start_idx = 0;
        if (it != cache_meteo.end())
            start_idx = std::distance(cache_meteo.begin(), it);

        if (start_idx > 0)
            start_idx--;

        for (size_t i = start_idx; i < cache_meteo.size(); i++) {
            for (const Date& current_date : cache_meteo[i].getDates()) {
                if (current_date > dateEnd) {
                    break;
                }
                if (!meta_ok) {
                    if (!readMeteoMeta(cache_meteo[i], vecPts, stations, lats, lons)) {
                        // some points have been removed vecPts has been changed -> re-reading
                        lats.clear();
                        lons.clear();
                        readMeteoMeta(cache_meteo[i], vecPts, stations, lats, lons);
                    }
                    vecvecMeteo.insert(vecvecMeteo.begin(), vecPts.size(), std::vector<MeteoData>()); // allocation for the vectors now that we know how many true stations we have
                    meta_ok = true;
                }
                
                std::vector<MeteoData> vecMeteoStations = createMeteoDataVector(stations, current_date);
                const size_t npoints = vecMeteoStations.size();

                for (size_t par_index=0; par_index <= MeteoData::Parameters::lastparam; par_index++) {
                    std::string param_name = MeteoData::getParameterName(par_index);
                    long level_no;
                    std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(cache_meteo[i], param_name, parameter_table, level_no, verbose);

                    if (messages.empty()) {
                        if (verbose) {
                            std::cout << "No messages found for parameter " << param_name << " in file " << cache_meteo[i].getFilename() << std::endl;
                        }
                        continue; 
                    }

                    processMessages(messages, level_no, vecMeteoStations, lats, lons, npoints, par_index);
                }
                for (size_t jj=0; jj<vecPts.size(); jj++)
                    vecvecMeteo[jj].push_back(vecMeteoStations[jj]);
            }
        }
    }


    // ---------------------------- METEO READING HELPERS -----------------------------
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

    bool GRIBIO::readMeteoMeta(GRIBFile& file ,std::vector<Coords>& vecPoints, std::vector<StationData> &stations, std::vector<double> &lats, std::vector<double> &lons) {

        long level_no;
        std::vector<CodesHandlePtr> messages = extractParameterInfoAndFindMessages(file, MeteoGrids::DEM, parameter_table, level_no, verbose);

        if (messages.empty()) {
            throw IOException("No messages containing DEM information found." AT);
        }
        if (messages.size() > 1 && verbose) {
            std::cout << "Multiple messages containing DEM information found in file" << file.getFilename() << std::endl;
        }
        const size_t npoints = vecPoints.size();


        CodesHandlePtr& h = messages.front();
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
        size_t n_lats = lats.size();
        if (n_lats != lons.size())
            throw InvalidArgumentException("lats and lons vectors must have the same size. Something went seriously wront", AT);
        std::vector<double> outlats_vec(n_lats), outlons_vec(n_lats), altitudes_vec(n_lats), distances_vec(n_lats);
        std::vector<int> indexes_vec(n_lats);
        getNearestValues_grib(h, lats, lons, outlats_vec, outlons_vec, distances_vec, altitudes_vec, indexes_vec);

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
            std::ostringstream ss2;
            ss2 << "GRIB point (" << indexes_vec[ii] % Ni << "," << indexes_vec[ii] / Ni << ")";
            sd.stationName = ss2.str();
            stations.push_back(sd);
        }
        return true;
    }  
} // namespace
