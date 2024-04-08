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
#include "libcodes.h"
#include <meteoio/FileUtils.h>

#include <string>

namespace mio {
    namespace codes {

        const std::map<std::string,std::string> PARAMETER_MAP{
            {"geometric_height","3008"},
            {"10m_wind_u","165"}, // is it necessary to specify 10 levl here then?
            {"10m_wind_v","166"}, // is it necessary to specify 10 levl here then?
            {"wind_direction","3031"},
            {"wind_speed","10"},
            {"pressure","54"},
            {"temperature_2m","167"}, // there is two in meteoio, both have 11 in grib 1 table 2 what is the correct way?
            {"dew_point_temperature","3017"},
            {"relative_humidity","157"},
            {"snow_surface_temperature","201203"}, //TODO: or 500170
            {"total_precipitation","228"},
            {"runoff","205"}, //TODO: is this the surf+subsurf runoff?
            {"SWE_m", "141"},
            {"SWE_kgm-2","228141"},
            {"snow_depth", "3066"},
            {"snow_density","33"},
            {"albedo","32,174..."}, //TODO: is it snow albedo, or uv, or...?
            {"long_wave_radiation_flux","3115"},
            {"medium_cloud_cover","187"}, // 0-1
            {"short_wave_radiation_flux","3116"},
            {"diffuse_radiation_albedo","228245"},
            {"direct_radiation_albedo","228244"},
            {"global_radiation_flux","3117"}, // is this glob->111.250 and whats the difference to 117.2?
            {"subgrid_slope","163"},
            {"subgrid_angle_eastward","162"},
            {"maximum_wind_velocity","201187"},
            {"geometric_vertival_velocity","260238"},
            {"surface_pressure","134"},
            {"water_equivalent_of_accumulated_snow_depth(deprecated)","260056"}, // what to use instead??
            {"surface_short_wave_radiation_downwards","169"}, // is this the correct one, or mean...
            {"surface_solar_radiation_difference","200176"} // or is it 200047?
        };

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, double& parameterValue) {
            CODES_CHECK(codes_get_double(h.get(), parameterName.c_str(), &parameterValue), 0);
        }
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, long& parameterValue) {
            CODES_CHECK(codes_get_long(h.get(), parameterName.c_str(), &parameterValue), 0);
        }

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, std::string& parameterValue) {
            size_t len = 500;
            char name[len] = {'\0'};
            CODES_CHECK(codes_get_string(h.get(), parameterName.c_str(), name, &len), 0);
            parameterValue = std::string(name);
        }

        // Function to get the list of parameters from a GRIB file, caller needs to free the index
        CodesIndexPtr indexFile(const std::string &filename, std::vector<std::string> &paramIdList, const long& ensembleNumber, bool verbose) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames

            int ret;
            size_t paramIdSize, values_len = 0;

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

//         CodesIndexPtr indexFile(const std::string &filename, std::vector<std::string> &paramIdList, std::vector<long> &ensembleNumbers, std::vector<long> &levelNumbers, std::vector<double> &datesList) {
//             if (!FileUtils::fileExists(filename))
//                 throw AccessException(filename, AT); // prevent invalid filenames

//             int ret;
//             size_t paramIdSize, numberSize, levelsSize, values_len = 0;

//             codes_index *index = codes_index_new_from_file(0, filename.c_str(), "paramId,number,typeOfLevel,dataDate", &ret);

//             CODES_CHECK(ret, 0);

//             CODES_CHECK(codes_index_get_size(index, "paramId", &paramIdSize), 0);

//             std::vector<char *> paramId(paramIdSize);
//             CODES_CHECK(codes_index_get_string(index, "paramId", paramId.data(), &paramIdSize), 0);

//             paramIdList.reserve(paramId.size());

//             for (char* cstr : paramId) {
//                 if (cstr != nullptr) {
//                     paramIdList.push_back(std::string(cstr));
//                     delete[] cstr;
//                 }
//             }
//             paramId.clear();

//             CODES_CHECK(codes_index_get_size(index, "number", &numberSize), 0);

//             ensembleNumbers.resize(numberSize);
//             CODES_CHECK(codes_index_get_long(index, "number", ensembleNumbers.data(), &numberSize), 0);
            
//             CODES_CHECK(codes_index_get_size(index, "typeOfLevel", &levelsSize), 0);

//             levelNumbers.resize(levelsSize);
//             CODES_CHECK(codes_index_get_long(index, "typeOfLevel", levelNumbers.data(), &levelsSize), 0);
// #ifdef DEBUG
//             std::cerr << "Found " << paramIdList.size() << " parameters in " << filename << "\n";
//             for (const std::string &param : paramIdList) {
//                 std::cerr << param << "\n";
//             }
//             std::cerr << "Found " << ensembleNumbers.size() << " ensemble numbers in " << filename << "\n";
//             for (const long &number : ensembleNumbers) {
//                 std::cerr << number << "\n";
//             }
//             std::cerr << "Found " << levelNumbers.size() << " level numbers in " << filename << "\n";
//             for (const long &number : levelNumbers) {
//                 std::cerr << number << "\n";
//             }
// #endif

//             return makeUnique(index);
//         }

        std::vector<CodesHandlePtr> getMessages(CodesIndexPtr &index, const std::string &paramID, const long &ensembleNumber, const std::string &levelType) {
            codes_index *raw = index.get();
            CODES_CHECK(codes_index_select_string(raw, "paramId", paramID.c_str()), 0);
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

        std::vector<CodesHandlePtr> getMessages(const std::string &filename, ProductKind product = PRODUCT_GRIB) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames
            errno = 0;
            FILE *fp = fopen(filename.c_str(),"r");
            if (fp==nullptr) {
                std::ostringstream ss;
                ss << "Error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
                throw AccessException(ss.str(), AT);
            }

            codes_handle *h = nullptr;
            int err = 0;
            std::vector<CodesHandlePtr> handles;
            while ((h = codes_handle_new_from_file(0, fp, product, &err)) != nullptr) {
                if (!h)
                    throw IOException("Unable to create grib handle from file", AT);
                if (err != 0) {
                    codes_handle_delete(h);
                    throw IOException("Error reading message: Errno " + std::to_string(err), AT);
                }
                handles.push_back(makeUnique(h));
            }
            fclose(fp);
            return handles;
        }


        Date getMessageDate(CodesHandlePtr &h, double &d1, double &d2, const double &tz_in) {
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


        std::map<std::string, double> getGridParameters(CodesHandlePtr &h_unique ) {
            // getting transformation parameters
            long Ni, Nj;
            getParameter(h_unique, "Ni", Ni);
            getParameter(h_unique, "Nj", Nj);

            double angleOfRotationInDegrees, latitudeOfSouthernPole, longitudeOfSouthernPole, latitudeOfNorthernPole, longitudeOfNorthernPole;
            getParameter(h_unique, "angleOfRotationInDegrees", angleOfRotationInDegrees);
            if (angleOfRotationInDegrees != 0.) {
                throw InvalidArgumentException("Rotated grids not supported!", AT);
            }

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
            long Ni = gridParams.at("Ni"), Nj = gridParams.at("Nj");
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