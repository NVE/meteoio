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
- needs documentatino
*/
#include <meteoio/plugins/libcodes.h>
#include <meteoio/FileUtils.h>

#include <string>
#include <cstring>

namespace mio {
    namespace codes {

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
            int err = codes_get_double(h.get(), parameterName.c_str(), &parameterValue);
            if (err != 0) {
                throw IOException("Error reading parameter " + parameterName + ": Errno " + std::to_string(err), AT);
            }
        }
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, long& parameterValue) {
            int err = codes_get_long(h.get(), parameterName.c_str(), &parameterValue);
            if (err != 0) {
                throw IOException("Error reading parameter " + parameterName + ": Errno " + std::to_string(err), AT);
            }
        }

        // casts long to int
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, int& parameterValue) {
            long tmp;
            int err = codes_get_long(h.get(), parameterName.c_str(), &tmp);
            if (err != 0) {
                throw IOException("Error reading parameter " + parameterName + ": Errno " + std::to_string(err), AT);
            }
            parameterValue = static_cast<int>(tmp);
        }

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, std::string& parameterValue) {
            size_t len = 500;
            char name[500] = {'\0'};
            int err = codes_get_string(h.get(), parameterName.c_str(), name, &len);
            if (err != 0) {
                throw IOException("Error reading parameter " + parameterName + ": Errno " + std::to_string(err), AT);
            }
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

        void unpackMessage(CodesHandlePtr& m) {
            /* We need to instruct ecCodes to expand the descriptors
            i.e. unpack the data values */
            CODES_CHECK(codes_set_long(m.get(),"unpack",1),0);  
        }

        // TODO: check if this gives the correct date, and what is d1, d2??
        Date getMessageDateGrib(CodesHandlePtr &h, const double &tz_in) {
            Date base;
            long validityDate, validityTime;
            getParameter(h, "validityDate", validityDate);
            getParameter(h, "validityTime", validityTime);

            const int year = static_cast<int>(validityDate / 10000), month = static_cast<int>(validityDate / 100 - year * 100), day = static_cast<int>(validityDate - month * 100 - year * 10000);
            const int hour = static_cast<int>(validityTime / 100), minutes = static_cast<int>(validityTime - hour * 100); // HACK: handle seconds!
            base.setDate(year, month, day, hour, minutes, tz_in);

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
            try {
                getParameter(h_unique, "angleOfRotationInDegrees", angleOfRotationInDegrees);
            } catch (...) {
                angleOfRotationInDegrees = 0.; // angle of rotation is not there
            }

            try {
                getParameter(h_unique, "latitudeOfSouthernPoleInDegrees", latitudeOfSouthernPole);
                getParameter(h_unique, "longitudeOfSouthernPoleInDegrees", longitudeOfSouthernPole);
            } catch (...) {
                latitudeOfSouthernPole = -90.; // default values for when it is not rotated
                longitudeOfSouthernPole = 0.;
            }
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
            if (gridParams.empty()) {
                gridParams = getGridParameters(h);
            }

            size_t values_len;
            CODES_CHECK(codes_get_size(h.get(), "values", &values_len), 0);
            double Ni = gridParams.at("Ni"), Nj = gridParams.at("Nj");
            if (values_len != (unsigned)(Ni * Nj)) {
                std::ostringstream ss;
                ss << "Declaring grid of size " << Ni << "x" << Nj << "=" << Ni * Nj << " ";
                ss << "but containing " << values_len << " values. This is inconsistent!";
                throw InvalidArgumentException(ss.str(), AT);
            }
            
            values.resize(values_len);
            GRIB_CHECK(codes_get_double_array(h.get(), "values", values.data(), &values_len), 0);

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
        void getNearestValues_grib(CodesHandlePtr &h, const std::vector<double> &in_lats, const std::vector<double> &in_lons, std::vector<double> &out_lats, std::vector<double> &out_lons, std::vector<double> &distances, std::vector<double> &values, std::vector<int> &indexes) {
            size_t npoints = in_lats.size();
            std::cout << "Finding nearest values\n";
            CODES_CHECK(codes_grib_nearest_find_multiple(h.get(), 0, in_lats.data(), in_lons.data(), static_cast<long>(npoints), out_lats.data(), out_lons.data(), values.data(), distances.data(), indexes.data()), 0);
            std::cout << "Found nearest values\n";
        }

        // ------------------- SETTERS -------------------
        void setMissingValue(CodesHandlePtr &message, double missingValue) {
            CODES_CHECK(codes_set_double(message.get(), "missingValue", missingValue), 0);
        }
        bool selectParameter(codes_index* raw, const std::string& param_key, const std::string& paramId) {
            return codes_index_select_string(raw, param_key.c_str(), paramId.c_str()) == 0;
        };
        
        bool selectParameter(codes_index* raw, const std::string& param_key, const double& paramId) {
            return codes_index_select_double(raw, param_key.c_str(), paramId) == 0;
        };
        bool selectParameter(codes_index* raw, const std::string& param_key, const long& paramId) {
            return codes_index_select_long(raw, param_key.c_str(), paramId) == 0;
        };

    } // namespace codes

} // namespace mio