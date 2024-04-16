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
#include <meteoio/plugins/GRIBFile.h>
#include <meteoio/dataClasses/MeteoData.h>

#include <fstream>


namespace mio {
    using namespace codes;

    // ------------------------- GRIBTable -------------------------
    const std::string GRIBTable::default_table = "doc/resources/GRIB_param.tbl"; 

    GRIBTable::GRIBTable(const std::string &in_filename)
        : filename(in_filename), param_indexing(), level_indexing(), param_table(), param_table_double(), param_table_long(), level_type_table(), level_no_table(),
          parameter_id_type(PARAM_TYPE::STRING), known_params() {
        // Open the file

        if (filename.empty()) {
#ifdef DEBUG
            std::cerr << "No GRIB table specified, using default table" << std::endl;
#endif
            filename = default_table;
        }

        init_known_params();
        readTable();
    }

    void GRIBTable::init_known_params() {
        for (size_t i=0; i<MeteoData::nrOfParameters; i++) {
            std::string param_name = MeteoData::getParameterName(i);
            known_params.insert(param_name);
            param_table[param_name] = "";
            param_table_double[param_name] = 0.0;
            param_table_long[param_name] = 0;
            level_type_table[param_name] = "";
            level_no_table[param_name] = 0.0;
        }
    };

    // ------------------------- reading the GRIB table -------------------------
    void GRIBTable::readTable() {
        std::ifstream file(filename);

        if (!file.is_open()) {
            throw AccessException("Could not open GRIB table: " + filename);
        }

        // Read the file line by line and parse it
        std::string line;
        std::vector<std::string> unknown_parameters;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#' || line[0] == ' ') {
                continue;
            }

            std::vector<std::string> values = IOUtils::split(line, ':');

            if (parseIndexing(values))
                continue;
            if (parseParamType(values))
                continue;

            // only comments or parameters left so fill the parameter table
            fillParameterTables(values, unknown_parameters);
        }
#ifdef DEBUG
        if (!unknown_parameters.empty()) {
            std::cerr << "Unknown parameters in GRIB table: ";
            for (const auto &p : unknown_parameters) {
                std::cerr << p << " ";
            }
            std::cerr << std::endl;
        }
#endif
    }

    bool GRIBTable::parseIndexing(const std::vector<std::string> &line_vals) {
        if (line_vals[0] == "index") {
            if (line_vals.size() != 3) {
                std::stringstream ss;
                ss << "Invalid number of indexing values to retrieve data, got: ";
                for (const auto &v : line_vals) {
                    ss << v << " ";
                }
                ss << std::endl;
                ss << "Expected: index:param:levelType" << std::endl;
                throw InvalidFormatException(ss.str(), AT);
            }
            param_indexing = line_vals[1];
            level_indexing = line_vals[2];
            return true;
        }
        return false;
    }

    bool GRIBTable::parseParamType(const std::vector<std::string> &line_vals) {
        if (line_vals[0] == "paramIdType") {
            if (line_vals.size() != 2) {
                std::stringstream ss;
                ss << "Invalid parameter type specification, got: ";
                for (const auto &v : line_vals) {
                    ss << v << " ";
                }
                ss << std::endl;
                ss << "Expected: paramIdType:[doubele|str|long]" << std::endl;
                throw InvalidFormatException(ss.str(), AT);
            }
            if (line_vals[1] == "double") {
                parameter_id_type = PARAM_TYPE::DOUBLE;
            } else if (line_vals[1] == "string") {
                parameter_id_type = PARAM_TYPE::STRING;
            } else if (line_vals[1] == "long") {
                parameter_id_type = PARAM_TYPE::LONG;
            } else {
                std::stringstream ss;
                ss << "Invalid parameter type specification, got: " << line_vals[1] << std::endl;
                ss << "Expected: paramIdType:[doubele|str]" << std::endl;
                throw InvalidFormatException(ss.str(), AT);
            }
            return true;
        }
        return false;
    }

    void GRIBTable::fillParameterTables(const std::vector<std::string> &line_vals, std::vector<std::string> &unknown_params) {
        std::string key, paramID, levelType;
        double levelNo;
        if (line_vals.size() != 4) {
            std::stringstream ss;
            ss << "Invalid number of values to retrieve data, got: ";
            for (const auto &v : line_vals) {
                ss << v << " ";
            }
            ss << std::endl;
            ss << "Expected: paramName:paramID:levelType:levelNo" << std::endl;
            throw InvalidFormatException(ss.str(), AT);
        }

        key = line_vals[0];
        paramID = line_vals[1];
        levelType = line_vals[2];
        if (!IOUtils::convertString(levelNo, line_vals[3])) {
            throw InvalidArgumentException("Invalid level number: " + line_vals[3], AT);
        };

        if (known_params.find(key) == known_params.end()) {
            unknown_params.push_back(key);
        }

        if (parameter_id_type == PARAM_TYPE::DOUBLE) {
            double paramID_double;
            if (!IOUtils::convertString(paramID_double, paramID)) {
                throw InvalidArgumentException("Could not convert parameter ID to number (as specified in paramType): " + paramID, AT);
            }
            param_table_double[key] = paramID_double;
        } else if (parameter_id_type == PARAM_TYPE::STRING) {
            param_table[key] = paramID;
        } else if (parameter_id_type == PARAM_TYPE::LONG) {
            long paramID_long;
            if (!IOUtils::convertString(paramID_long, paramID)) {
                throw InvalidArgumentException("Could not convert parameter ID to number (as specified in paramType): " + paramID, AT);
            }
            param_table_long[key] = paramID_long;
        } else {
            throw InvalidArgumentException("Invalid parameter type, if this appears, please contact the developers", AT);
        }
        level_type_table[key] = levelType;
        level_no_table[key] = levelNo;
    }

    // ------------------------- GETTERS -------------------------
    std::vector<std::string> GRIBTable::getIndexes() const {
        return {param_indexing, level_indexing};
    }

    void GRIBTable::getParamId(const std::string &param_name, std::string &paramId, double &paramId_num, long &paramId_long) const {
        if (parameter_id_type == PARAM_TYPE::DOUBLE) {
            paramId_num = param_table_double.at(param_name);
            paramId = "";
            paramId_long = 0;
        } else if (parameter_id_type == PARAM_TYPE::STRING) {
            paramId = param_table.at(param_name);
            paramId_num = 0.0;
            paramId_long = 0;
        } else if (parameter_id_type == PARAM_TYPE::LONG) {
            paramId_long = param_table_long.at(param_name);
            paramId = "";
            paramId_num = 0.0;
        } else {
            throw InvalidArgumentException("Invalid parameter type, if this appears, please contact the developers", AT);
        }   
    }

    std::string GRIBTable::getLevelType(const std::string &param_name) const { return level_type_table.at(level_name); }

    double GRIBTable::getLevelNo(const std::string &param_name) const { return level_no_table.at(level); }



    // ------------------------- GRIBFile -------------------------
    GRIBFile::GRIBFile(const std::string &in_filename, const std::vector<std::string> &indexes)
        : filename(in_filename), file(), grid_params() {

            file = indexFile(filename, indexes, false); // true = verbose output
            checkValidity();
        }
    
    bool GRIBFile::checkValidity() {
        std::vector<CodesHandlePtr> messages = getMessages(filename, PRODUCT_GRIB);
        if (messages.empty()) {
            throw AccessException("No messages found in GRIB file: " + filename);
        }

        Date timepoint;
        for (auto &m : messages) {
            std::map<std::string, double> new_grid_params = getGridParameters(m);
            if (grid_params.empty()) {
                grid_params = new_grid_params;
            } else if (!new_grid_params.empty() && grid_params != new_grid_params) {
                    throw InvalidFormatException("Grid parameters do not match in GRIB file: " + filename);
            }

            double d1, d2;
            if (timepoint == Date()) {
                timepoint = getMessageDateGrib(m, d1, d2, 0);
            } else {
                Date new_timepoint = getMessageDateGrib(m, d1, d2, 0);
                if (timepoint != new_timepoint && new_timepoint != Date()) {
                    throw InvalidFormatException("There is more than 1 timepoint in GRIB file: " + filename);
                }
            
            }
        }
        if (grid_params.empty()) {
            throw InvalidFormatException("No grid parameters found in GRIB file: " + filename);
        }
    }  
} // namespace mio