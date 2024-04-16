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
- needs documentatino and refactoring
*/
#include <meteoio/plugins/libcodes.h>
#include <meteoio/FileUtils.h>

#include <string>
#include <cstring>

namespace mio {
    namespace codes {

        // ------------------- GRIB PARAMETER MAPPING -------------------
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

        // ------------------- BUFR -------------------
        const std::map<std::string, std::string> BUFR_PARAMETER {
                {"P","pressure"}, 
                {"TA","airTemperature"}, // TODO: or is it airTemperatureAt2M	
                {"RH","relativeHumidity"}, 
                {"TSG","groundTemperature"}, 
                {"TSS","snowTemperature"}, 
                {"HS","totalSnowDepth"}, 
                {"VW","windSpeed"}, 
                {"DW","windDirection"},
                {"VW_MAX","maximumWindSpeedMeanWind"}, //TODO: or is it maximumWindGustSpeed
                {"RSWR",""}, // TODO: which parameters are tey
                {"ISWR",""}, // TODO: which parameters are tey
                {"ILWR",""}, // TODO: which parameters are tey
                {"TAU_CLD","cloudCoverTotal	"}, // TODO: is in per_cent
                {"PSUM","totalPrecipitationOrTotalWaterEquivalent	"},  
                {"PSUM_PH","precipitationType"} // should the type be mapped to the phase?
        };

        static const std::vector<int> FLAG_TO_EPSG = {4326, 4258, 4269, 4314};

        // ------------------- POINTER HANDLING -------------------
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

        // ------------------- FILE HANDLING -------------------
        // Function to get the list of parameters from a GRIB file, caller needs to free the index
        CodesIndexPtr indexFile(const std::string &filename, const std::vector<std::string>& index_keys, bool verbose) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames

            std::string indexing_string;
            for (const std::string &key : index_keys) {
                if (!indexing_string.empty()) {
                    indexing_string += ",";
                }
                indexing_string += key;
            }

            int ret;
            codes_index *index = codes_index_new_from_file(0, filename.c_str(), indexing_string.c_str(), &ret);

            CODES_CHECK(ret, 0);

            if (verbose) {
                std::cerr << "Indexing " << filename << " with keys " << indexing_string << "\n";
                for (const std::string &key : index_keys) {
                    size_t size;
                    CODES_CHECK(codes_index_get_size(index, key.c_str(), &size), 0);
                    std::cerr << "Found " << size << " " << key << " in " << filename << "\n";

                    std::vector<char *> values(size);
                    CODES_CHECK(codes_index_get_string(index, key.c_str(), values.data(), &size), 0);

                    std::cerr << "Values:\n";
                    for (char* cstr : values) {
                        if (cstr != nullptr) {
                            std::cerr << cstr << "\n";
                            delete[] cstr;
                        }
                    }   


                }
            }
            return makeUnique(index);
        }

        // ------------------- GETTERS -------------------
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, double& parameterValue) {
            CODES_CHECK(codes_get_double(h.get(), parameterName.c_str(), &parameterValue), 0);
        }
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, long& parameterValue) {
            CODES_CHECK(codes_get_long(h.get(), parameterName.c_str(), &parameterValue), 0);
        }

        // casts long to int
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, int& parameterValue) {
            long tmp;
            CODES_CHECK(codes_get_long(h.get(), parameterName.c_str(), &tmp), 0);
            parameterValue = static_cast<int>(tmp);
        }

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, std::string& parameterValue) {
            size_t len = 500;
            char name[500] = {'\0'};
            CODES_CHECK(codes_get_string(h.get(), parameterName.c_str(), name, &len), 0);
            parameterValue = std::string(name);
        }

        // ------------------- MESSAGE HANDLING -------------------
        std::vector<CodesHandlePtr> getMessages(const std::string &filename, ProductKind product) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames
            errno = 0;
            FILE *fp = fopen(filename.c_str(),"r");
            
            return getMessages(fp, product);
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

        // TODO: check if this gives the correct date, and what is d1, d2??
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

        // assumes UTC 0, need to convert to local time after
        Date getMessageDateBUFR(CodesHandlePtr &h, const double &tz_in){
            std::vector<std::string> parameters = {"year", "month", "day", "hour", "minute", "second"};
            std::vector<int> values(parameters.size(), -1);

            for (size_t i = 0; i < parameters.size(); ++i) {
                getParameter(h, parameters[i], values[i]);
            }

            if (values[0] == -1 || values[1] == -1 || values[2] == -1 || values[3] == -1) return Date();
            if (values[5] == -1) return Date(values[0], values[1], values[2], values[3], values[4], tz_in);
            
            Date base(values[0], values[1], values[2], values[3], values[4], values[5], tz_in);
            return base;
        };

        StationData getStationDataBUFR(CodesHandlePtr &h, const std::string &ref_coords, std::string& error) {
            StationData stationData;

            double latitude, longitude, altitude;
            getParameter(h, "latitude", latitude);
            getParameter(h, "longitude", longitude);

            std::vector<std::string> height_keys = {"heightOfStation","height","elevation"};
            getParameter(h, height_keys, altitude);

            std::vector<std::string> id_keys = {"stationNumber", "nationalStationNumber", "stationID","stationId"};
            std::vector<std::string> name_keys = {"shortStationName","stationOrSiteName","longStationName"};
            std::string stationID, stationName;
            getParameter(h, id_keys, stationID);
            getParameter(h, name_keys, stationName);

            long ref_flag = 999;
            getParameter(h, "coordinateReferenceSystem", ref_flag);
            Coords position;

            if (ref_flag == 999 && ref_coords.empty()) {
                error = "No reference coordinates found in BUFR file and none provided in the configuration";
                return stationData;
            }

            if(ref_flag == 65535) {
                if (ref_coords.empty()) {
                    error = "Missing reference coordinates in BUFR file, and none provided in configuration";
                    return stationData;
                } else {
                    Coords ref_coords_obj(ref_coords);
                    ref_coords_obj.setLatLon(latitude, longitude, altitude);
                    position = ref_coords_obj;
                }
            } else {
                if (ref_flag > 3) {
                    std::ostringstream ss;
                    ss << "Unsuppported reference flag " << ref_flag << " in BUFR file";
                    throw InvalidFormatException(ss.str(), AT);
                } else {
                    position.setEPSG(FLAG_TO_EPSG[ref_flag]);
                    position.setPoint(latitude, longitude, altitude);
                }
            }
            
            stationData.setStationData(position, stationID, stationName);            
            return stationData;
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