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

/* https://codes.ecmwf.int/grib/param-db/?ordering=id&encoding=grib2 -- parameter database

*/

#include "libcodes.h"
#include <meteoio/FileUtils.h>

#include <string>

namespace mio {
    namespace codes {

        void getParameter(codes_handle* h, const std::string& parameterName, double& parameterValue) {
            CODES_CHECK(codes_get_double(h, parameterName.c_str(), &parameterValue), 0);
        }
        void getParameter(codes_handle* h, const std::string& parameterName, long& parameterValue) {
            CODES_CHECK(codes_get_long(h, parameterName.c_str(), &parameterValue), 0);
        }
        void getParameter(codes_handle* h, const std::string& parameterName, std::string& parameterValue) {
            size_t len = 500;
            char name[len] = {'\0'};
            CODES_CHECK(codes_get_string(h, parameterName.c_str(), name, &len), 0);
            parameterValue = std::string(name);
        }

        // Function to get the list of parameters from a GRIB file, caller needs to free the index
        CodesIndexPtr getListOfParamaters(const std::string &filename, std::vector<std::string> &paramIdList, std::vector<long> &ensembleNumbers) {
            if (!FileUtils::fileExists(filename))
                throw AccessException(filename, AT); // prevent invalid filenames

            int ret;
            size_t paramIdSize, numberSize, values_len = 0;

            codes_index *index = codes_index_new_from_file(0, filename.c_str(), "paramId,number,indicatorOfTypeOfLevel", &ret);

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

            CODES_CHECK(codes_index_get_size(index, "number", &numberSize), 0);

            std::vector<long> ensembleNumbers(numberSize);
            CODES_CHECK(codes_index_get_long(index, "number", ensembleNumbers.data(), &numberSize), 0);
#ifdef DEBUG
            std::cerr << "Found " << paramIdList.size() << " parameters in " << filename << "\n";
            for (const std::string &param : paramIdList) {
                std::cerr << param << "\n";
            }
#endif

            return makeUnique(index);
        }

        Date getDate(CodesHandlePtr h, double &d1, double &d2, const double &tz_in) {
            codes_handle *raw = h.get();
            Date base;
            long dataDate, dataTime;
            getParameter(raw, "dataDate", dataDate);
            getParameter(raw, "dataTime", dataTime);

            const int year = static_cast<int>(dataDate / 10000), month = static_cast<int>(dataDate / 100 - year * 100), day = static_cast<int>(dataDate - month * 100 - year * 10000);
            const int hour = static_cast<int>(dataTime / 100), minutes = static_cast<int>(dataTime - hour * 100); // HACK: handle seconds!
            base.setDate(year, month, day, hour, minutes, tz_in);

            // reading offset to base date/time, as used for forecast, computed at time t for t+offset
            long startStep, endStep;
            std::string stepUnits;
            getParameter(raw, "stepUnits", stepUnits);
            getParameter(raw, "startStep", startStep);
            getParameter(raw, "endStep", endStep);

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
        }

        std::vector<CodesHandlePtr> getMessages(CodesIndexPtr index, const std::string &paramID, const long &ensembleNumber, const long &levelType) {
            codes_index *raw = index.get();
            CODES_CHECK(codes_index_select_string(raw, "paramId", paramID.c_str()), 0);
            CODES_CHECK(codes_index_select_long(raw, "number", ensembleNumber), 0);
            CODES_CHECK(codes_index_select_long(raw, "indicatorOfTypeOfLevel", levelType), 0);

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

        std::map<std::string, double> getGridParameters(codes_handle* h ) {
            // getting transformation parameters
            double angleOfRotationInDegrees, latitudeOfSouthernPole, longitudeOfSouthernPole, latitudeOfNorthernPole, longitudeOfNorthernPole;
            getParameter(h, "angleOfRotationInDegrees", angleOfRotationInDegrees);
            if (angleOfRotationInDegrees != 0.) {
                throw InvalidArgumentException("Rotated grids not supported!", AT);
            }

            getParameter(h, "latitudeOfSouthernPoleInDegrees", latitudeOfSouthernPole);
            getParameter(h, "longitudeOfSouthernPoleInDegrees", longitudeOfSouthernPole);
            latitudeOfNorthernPole = -latitudeOfSouthernPole;
            longitudeOfNorthernPole = longitudeOfSouthernPole + 180.;

            // determining llcorner, urcorner and center coordinates
            double ll_latitude, ll_longitude, ur_latitude, ur_longitude;
            getParameter(h, "latitudeOfFirstGridPointInDegrees", ll_latitude);
            getParameter(h, "longitudeOfFirstGridPointInDegrees", ll_longitude);

            getParameter(h, "latitudeOfLastGridPointInDegrees", ur_latitude);
            getParameter(h, "longitudeOfLastGridPointInDegrees", ur_longitude);

            std::map<std::string, double> gridParams = {
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

        void getGriddedValues(CodesHandlePtr h, std::vector<double> &values, long &Ni, long &Nj, std::map<std::string,double> &gridParams) {
            codes_handle *raw = h.get();
            size_t values_len;
            CODES_CHECK(codes_get_size(raw, "values", &values_len), 0);
            getParameter(raw, "Ni", Ni);
            getParameter(raw, "Nj", Nj);
            if (values_len != (unsigned)(Ni * Nj)) {
                std::ostringstream ss;
                ss << "Declaring grid of size " << Ni << "x" << Nj << "=" << Ni * Nj << " ";
                ss << "but containing " << values_len << " values. This is inconsistent!";
                throw InvalidArgumentException(ss.str(), AT);
            }
            
            gridParams = getGridParameters(raw);

            values.resize(values_len);
            GRIB_CHECK(codes_get_double_array(raw, "values", values.data(), &values_len), 0);

        }

        // Retrieve Meteo Specific Data

    } // namespace codes

} // namespace mio