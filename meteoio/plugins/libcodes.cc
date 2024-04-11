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

/* https://codes.ecmwf.int/grib/param-db/?ordering=id-- parameter database

*/

/*
NOTES:
- so far this is grib specific, will change to include bufr as well
*/
#include <meteoio/plugins/libcodes.h>
#include <meteoio/FileUtils.h>

#include <string>
#include <cstring>

namespace mio {
    namespace codes {

        // Mapping of parameter names to their GRIB codes
        const std::map<std::string,std::string> PARAMETER_MAP{
            {"dem_geometric_height","300008"},
            {"dem_model_terrain_height","260183"},
            {"dem_geometric_height_above_sea_level","500007"},
            {"dem_orography","228002"},
            {"wind_u_10m","165"}, 
            {"wind_u","131"}, 
            {"wind_v_10m","166"}, 
            {"wind_v","132"}, 
            {"wind_direction","3031"},
            {"wind_direction_10m","260260"},
            {"wind_direction_v2","300031"},
            {"wind_direction_dd_10m","500023"},
            {"wind_direction_dd","500024"},
            {"wind_speed","10"},
            {"wind_speed_10m","207"},
            {"wind_speed_10m_2","160246"},
            {"pressure","54"},
            {"temperature_2m","167"}, // there is two in meteoio, both have 11 in grib 1 table 2 what is the correct way?
            {"temperature","130"},
            {"dew_point_temperature","3017"},
            {"dew_point_temperature_2","300017"},
            {"dew_point_temperature_2m","500017"},
            {"relative_humidity_per_cent","157"},
            {"relative_humidity_int","160157"},
            {"snow_surface_temperature_top","500170"},
            {"snow_surface_temperature","201203"}, 
            {"total_precipitation","228"},
            {"runoff","205"}, //TODO: is this the surf+subsurf runoff?
            {"SWE_m", "141"},
            {"SWE_kgm-2","228141"},
            {"snow_depth", "3066"},
            {"snow_density","33"},
            {"albedo_per_cent","500056"},
            {"albedo_backup_int","243"}, 
            {"albedo_backup_2_per_cent", "160243"},
            {"albedo_backup_3_per_cent","500057"},
            {"albedo_backup_4_per_cent","260509"},
            {"long_wave_radiation_flux","3115"},
            {"net_long_wave_radiation_flux","260099"}, // TODO: probably need a conversion here 
            {"medium_cloud_cover","187"}, // 0-1
            {"medium_cloud_cover_per_cent","3074"},
            {"short_wave_radiation_flux","3116"},
            {"short_wave_radiation_flux_at_surface","174116"},
            {"short_wave_radiation_flux_direct","260262"},
            {"short_wave_radiation_flux_diffuse","260263"},
            {"short_wave_net_radiation_flux","500079"},
            {"short_wave_net_radiation_flux_2","500083"},
            {"global_radiation_flux","3117"}, 
            {"subgrid_slope","163"},
            {"subgrid_angle_eastward","162"},
            {"maximum_wind_velocity","201187"},
            {"vertical_velocity","135"},
            {"vertical_velocity_v1","300040"},
            {"geometric_vertical_velocity","500032"},
            {"geometric_vertival_velocity_2","260238"},
            {"surface_pressure","134"},
            {"logarithm_of_surface_pressure","152"},
            {"surface_short_wave_radiation_downwards","169"}, 
            {"surface_long_wave_radiation_downwards","175"},
            {"surface_solar_radiation_difference","200169"} 
        };

        // clustering parameters that will give exactly the same results
        const std::map<std::string,std::vector<std::string>> PARAMETER_GROUPS {
            {"DEM",{PARAMETER_MAP.at("dem_geometric_height"), PARAMETER_MAP.at("dem_model_terrain_height"), PARAMETER_MAP.at("dem_geometric_height_above_sea_level"), PARAMETER_MAP.at("dem_orography")}},
            {"wind_direction",{PARAMETER_MAP.at("wind_direction"), PARAMETER_MAP.at("wind_direction_v2"), PARAMETER_MAP.at("wind_direction_dd")}},
            {"wind_direction_10m",{PARAMETER_MAP.at("wind_direction_10m"), PARAMETER_MAP.at("wind_direction_dd_10m")}},
            {"wind_speed_10m",{PARAMETER_MAP.at("wind_speed_10m"), PARAMETER_MAP.at("wind_speed_10m_2")}},
            {"dew_point_temperature",{PARAMETER_MAP.at("dew_point_temperature"), PARAMETER_MAP.at("dew_point_temperature_2")}},
            {"snow_surface_temperature",{PARAMETER_MAP.at("snow_surface_temperature_top"), PARAMETER_MAP.at("snow_surface_temperature")}},
            {"albedo_per_cent",{PARAMETER_MAP.at("albedo_per_cent"), PARAMETER_MAP.at("albedo_backup_2_per_cent"), PARAMETER_MAP.at("albedo_backup_3_per_cent"), PARAMETER_MAP.at("albedo_backup_4_per_cent")}},
            {"short_wave_radiation_flux",{PARAMETER_MAP.at("short_wave_radiation_flux"), PARAMETER_MAP.at("short_wave_radiation_flux_at_surface")}},
            {"vertical_velocity",{PARAMETER_MAP.at("vertical_velocity"), PARAMETER_MAP.at("vertical_velocity_v1"), PARAMETER_MAP.at("geometric_vertical_velocity"),PARAMETER_MAP.at("geometric_vertical_velocity_2")}}
        };

        CodesHandlePtr makeUnique(codes_handle *h) {
            CodesHandlePtr ptr(h);
            h = nullptr;
            return ptr;
        }
        CodesIndexPtr makeUnique(codes_index *i) {
            CodesIndexPtr ptr(i);
            i = nullptr;
            return ptr;
        }

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, double& parameterValue) {
            CODES_CHECK(codes_get_double(h.get(), parameterName.c_str(), &parameterValue), 0);
        }
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, long& parameterValue) {
            CODES_CHECK(codes_get_long(h.get(), parameterName.c_str(), &parameterValue), 0);
        }

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, std::string& parameterValue) {
            size_t len = 500;
            char name[500] = {'\0'};
            CODES_CHECK(codes_get_string(h.get(), parameterName.c_str(), name, &len), 0);
            parameterValue = std::string(name);
        }

        // Function to get the list of parameters from a GRIB file, caller needs to free the index
        CodesIndexPtr indexFile(const std::string &filename, std::vector<std::string> &paramIdList, const long& ensembleNumber, bool verbose) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames

            int ret;
            size_t paramIdSize = 0;

            std::string index_str = "paramId,typeOfLevel";
            if (ensembleNumber != -1) index_str += ",number";
            codes_index *index = codes_index_new_from_file(0, filename.c_str(), index_str.c_str(), &ret);

            CODES_CHECK(ret, 0);

            CODES_CHECK(codes_index_get_size(index, "paramId", &paramIdSize), 0);

            std::vector<char *> paramId(paramIdSize);
            CODES_CHECK(codes_index_get_string(index, "paramId", paramId.data(), &paramIdSize), 0);

            paramIdList.reserve(paramId.size());

            for (char* cstr : paramId) {
                if (cstr != nullptr) {
                    paramIdList.push_back(std::string(cstr));
                    delete[] cstr;
                }
            }
            paramId.clear();
            if (verbose) {
                std::cerr << "Found " << paramIdList.size() << " parameters in " << filename << "\n";
                for (const std::string &param : paramIdList) {
                    std::cerr << param << "\n";
                }
            
            }


#ifdef DEBUG
            size_t numberSize, levelsSize;
            std::vector<long> ensembleNumbers;
            std::vector<std::string> levelNumbers;
            CODES_CHECK(codes_index_get_size(index, "number", &numberSize), 0);

            ensembleNumbers.resize(numberSize);
            CODES_CHECK(codes_index_get_long(index, "number", ensembleNumbers.data(), &numberSize), 0);
            
            CODES_CHECK(codes_index_get_size(index, "typeOfLevel", &levelsSize), 0);

            std::vector<char *> levelTypes(levelsSize);
            CODES_CHECK(codes_index_get_string(index, "typeOfLevel", levelTypes.data(), &levelsSize), 0);

            for (char* cstr : levelTypes) {
                if (cstr != nullptr) {
                    levelNumbers.push_back(std::string(cstr));
                    delete[] cstr;
                }
            }
            levelTypes.clear();
            std::cerr << "Found " << paramIdList.size() << " parameters in " << filename << "\n";
            for (const std::string &param : paramIdList) {
                std::cerr << param << "\n";
            }
            std::cerr << "Found " << ensembleNumbers.size() << " ensemble numbers in " << filename << "\n";
            for (const long &number : ensembleNumbers) {
                std::cerr << number << "\n";
            }
            std::cerr << "Found " << levelNumbers.size() << " level numbers in " << filename << "\n";
            for (const std::string &number : levelNumbers) {
                std::cerr << number << "\n";
            }
#endif

            return makeUnique(index);
        }

        std::vector<CodesHandlePtr> getMessages(CodesIndexPtr &index, const std::string &paramID, const long &ensembleNumber, const std::string &levelType) {
            codes_index *raw = index.get();
            if (codes_index_select_string(raw, "paramId", paramID.c_str()) != 0) {
                return {};
            };
            if (ensembleNumber != -1)
                CODES_CHECK(codes_index_select_long(raw, "number", ensembleNumber), 0);
            CODES_CHECK(codes_index_select_string(raw, "typeOfLevel", levelType.c_str()), 0);

            codes_handle *h = nullptr;
            int ret;
            std::vector<CodesHandlePtr> handles;
            while ((h = codes_handle_new_from_index(raw, &ret)) != nullptr) {
                if (!h)
                    throw IOException("Unable to create grib handle from index", AT);
                if (ret != 0 && ret != CODES_END_OF_INDEX) {
                    throw IOException("Error reading message: Errno " + std::to_string(ret), AT);
                }
#ifdef DEBUG
                size_t len = 500;
                char name[len] = {'\0'};
                std::cerr << "Found message for " << paramID << " " << ensembleNumber << " " << levelType << "\n";
                CODES_CHECK(codes_get_string(h, "name", name, &len), 0);
                std::cerr << "With name " << name << "\n";
                CODES_CHECK(codes_get_string(h, "shortName", name, &len), 0);
                std::cerr << "With shortName " << name << "\n";
                if (levelType != 1) {
                    long level;
                    GRIB_CHECK(codes_get_long(h, "level", &level), 0);
                    std::cerr << "With level " << level << "\n";
                }
#endif
                handles.push_back(makeUnique(h));
            }
            return handles;
        }

        std::vector<CodesHandlePtr> getMessages(CodesIndexPtr &index, const std::vector<std::string> &paramID_list, const long &ensembleNumber, const std::string &levelType) {
            for (const std::string &paramID : paramID_list) {
                std::vector<CodesHandlePtr> handles = getMessages(index, paramID, ensembleNumber, levelType);
                if (!handles.empty())
                    return handles;
            }
        }

        std::vector<CodesHandlePtr> getMessages(const std::string &filename, ProductKind product) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames
            errno = 0;
            FILE *fp = fopen(filename.c_str(),"r");
            if (fp==nullptr) {
                std::ostringstream ss;
                ss << "Error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
                throw AccessException(ss.str(), AT);
            }
            std::vector<CodesHandlePtr> handles = getMessages(fp, product);
           
            fclose(fp);
            return handles;
        }

        std::vector<CodesHandlePtr> getMessages(FILE* in_file, ProductKind product) {
            codes_handle *h = nullptr;
            int err = 0;
            std::vector<CodesHandlePtr> handles;
            while ((h = codes_handle_new_from_file(0, in_file, product, &err)) != nullptr) {
                if (!h)
                    throw IOException("Unable to create handle from file", AT);
                if (err != 0) {
                    codes_handle_delete(h);
                    throw IOException("Error reading message: Errno " + std::to_string(err), AT);
                }
                handles.push_back(makeUnique(h));
            }
            return handles;
        }

        Date getMessageDateGrib(CodesHandlePtr &h, double &d1, double &d2, const double &tz_in) {
            Date base;
            long validityDate, validityTime;
            getParameter(h, "validityDate", validityDate);
            getParameter(h, "validityTime", validityTime);

            const int year = static_cast<int>(validityDate / 10000), month = static_cast<int>(validityDate / 100 - year * 100), day = static_cast<int>(validityDate - month * 100 - year * 10000);
            const int hour = static_cast<int>(validityTime / 100), minutes = static_cast<int>(validityTime - hour * 100); // HACK: handle seconds!
            base.setDate(year, month, day, hour, minutes, tz_in);

            // reading offset to base date/time, as used for forecast, computed at time t for t+offset
            long startStep, endStep;
            std::string stepUnits;
            getParameter(h, "stepUnits", stepUnits);
            getParameter(h, "startStep", startStep);
            getParameter(h, "endStep", endStep);

            double step_units; // in julian, ie. in days

            std::string stepUnitsStr(stepUnits);
            if (stepUnitsStr == "s") {
                step_units = 1. / (24. * 60. * 60.);
            } else if (stepUnitsStr == "m") {
                step_units = 1. / (24. * 60.);
            } else if (stepUnitsStr == "h") {
                step_units = 1. / 24.;
            } else if (stepUnitsStr == "D") {
                step_units = 1.;
            } else if (stepUnitsStr == "3H") {
                step_units = 3. / 24.;
            } else if (stepUnitsStr == "6H") {
                step_units = 6. / 24.;
            } else if (stepUnitsStr == "12H") {
                step_units = 12. / 24.;
            } else if (stepUnitsStr == "M") {
                step_units = 30.44;
            } else if (stepUnitsStr == "Y") {
                step_units = 365.25;
            } else if (stepUnitsStr == "10Y") {
                step_units = 10 * 365.25;
            } else if (stepUnitsStr == "C") {
                step_units = 10 * 100 * 365.25;
            } else {
                std::ostringstream ss;
                ss << "GRIB file using stepUnits=" << stepUnits << ", which is not supported";
                throw InvalidFormatException(ss.str(), AT);
            }

            d1 = static_cast<double>(startStep) * step_units;
            d2 = static_cast<double>(endStep) * step_units;
            return base;
        }

        Date getMessageDateBUFR(CodesHandlePtr &h){

        };

        StationData getStationData(CodesHandlePtr &h) {

        };


        std::map<std::string, double> getGridParameters(CodesHandlePtr &h_unique ) {
            // getting transformation parameters
            long Ni, Nj;
            getParameter(h_unique, "Ni", Ni);
            getParameter(h_unique, "Nj", Nj);

            double angleOfRotationInDegrees, latitudeOfSouthernPole, longitudeOfSouthernPole, latitudeOfNorthernPole, longitudeOfNorthernPole;
            getParameter(h_unique, "angleOfRotationInDegrees", angleOfRotationInDegrees);

            getParameter(h_unique, "latitudeOfSouthernPoleInDegrees", latitudeOfSouthernPole);
            getParameter(h_unique, "longitudeOfSouthernPoleInDegrees", longitudeOfSouthernPole);
            latitudeOfNorthernPole = -latitudeOfSouthernPole;
            longitudeOfNorthernPole = longitudeOfSouthernPole + 180.;

            // determining llcorner, urcorner and center coordinates
            double ll_latitude, ll_longitude, ur_latitude, ur_longitude;
            getParameter(h_unique, "latitudeOfFirstGridPointInDegrees", ll_latitude);
            getParameter(h_unique, "longitudeOfFirstGridPointInDegrees", ll_longitude);

            getParameter(h_unique, "latitudeOfLastGridPointInDegrees", ur_latitude);
            getParameter(h_unique, "longitudeOfLastGridPointInDegrees", ur_longitude);

            std::map<std::string, double> gridParams = {
                {"Ni", static_cast<double>(Ni)},
                {"Nj", static_cast<double>(Nj)},
                {"Nj", Nj},
                {"ll_latitude", ll_latitude},
                {"ll_longitude", ll_longitude},
                {"ur_latitude", ur_latitude},
                {"ur_longitude", ur_longitude},
                {"angleOfRotationInDegrees", angleOfRotationInDegrees},
                {"latitudeOfSouthernPole", latitudeOfSouthernPole},
                {"longitudeOfSouthernPole", longitudeOfSouthernPole},
                {"latitudeOfNorthernPole", latitudeOfNorthernPole},
                {"longitudeOfNorthernPole", longitudeOfNorthernPole}
            };
            return gridParams;
        }

        std::map<std::string,double> getGriddedValues(CodesHandlePtr &h, std::vector<double> &values) {
            std::map<std::string, double> gridParams = getGridParameters(h);
            getGriddedValues(h, values, gridParams);
            return gridParams;
        }

        void getGriddedValues(CodesHandlePtr &h, std::vector<double> &values, std::map<std::string,double> &gridParams) {
            codes_handle *raw = h.get();

            if (gridParams.empty()) {
                gridParams = getGridParameters(h);
            }

            size_t values_len;
            CODES_CHECK(codes_get_size(raw, "values", &values_len), 0);
            double Ni = gridParams.at("Ni"), Nj = gridParams.at("Nj");
            if (values_len != (unsigned)(Ni * Nj)) {
                std::ostringstream ss;
                ss << "Declaring grid of size " << Ni << "x" << Nj << "=" << Ni * Nj << " ";
                ss << "but containing " << values_len << " values. This is inconsistent!";
                throw InvalidArgumentException(ss.str(), AT);
            }
            
            values.resize(values_len);
            GRIB_CHECK(codes_get_double_array(raw, "values", values.data(), &values_len), 0);

        }

        // Retrieve Meteo Specific Data
        /**
         * Retrieves the nearest values from a GRIB message based on the given latitude and longitude coordinates.
         *
         * @param h The handle to the GRIB file.
         * @param in_lats The input latitude coordinates.
         * @param in_lons The input longitude coordinates.
         * @return A tuple containing vectors with.
         *         - the output latitude coordinates, 
         *         - output longitude coordinates,
         *         - distances between the points,
         *         - values at the points,
         *         - indexes of the nearest values.
         */
        std::tuple<std::vector<double>&&, std::vector<double>&&, std::vector<double>&&, std::vector<double>&&, std::vector<int>&&> getNearestValues_grib(CodesHandlePtr &h, const std::vector<double> &in_lats, const std::vector<double> &in_lons) {
            codes_handle *raw = h.get();

            size_t npoints = in_lats.size();
            std::vector<double> out_lats(npoints), out_lons(npoints), distances(npoints), values(npoints);
            std::vector<int> indexes(npoints);


            CODES_CHECK(codes_grib_nearest_find_multiple(raw, 0, in_lats.data(), in_lons.data(), static_cast<long>(npoints), out_lats.data(), out_lons.data(), values.data(), distances.data(), indexes.data()), 0);

            return std::make_tuple(std::move(out_lats), std::move(out_lons), std::move(distances), std::move(values), std::move(indexes));
        }

    } // namespace codes

} // namespace mio