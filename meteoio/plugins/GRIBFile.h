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

#ifndef GRIBFILE_H
#define GRIBFILE_H

#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/plugins/libcodes.h>

#include <string>
#include <vector>
#include <unordered_set>

namespace mio
{

using namespace codes;

class GRIBTable {

    enum class PARAM_TYPE {
        DOUBLE,
        LONG,
        STRING,
    };

    public:
        GRIBTable(const std::string &in_filename);

        // GETTERS
        std::vector<std::string> getIndexes() const;

        void GRIBTable::getParamId(const std::string &param_name, std::string &paramId, double &paramId_num, long &paramId_long) const;
        std::string getLevelType(const std::string &level_name) const;
        double getLevelNo(const std::string &level) const;
    
    private:
        std::string filename;

        std::string param_indexing;
        std::string level_indexing;

        std::map<std::string, std::string> param_table;
        std::map<std::string, double> param_table_double;
        std::map<std::string, long> param_table_long;
        std::map<std::string, std::string> level_type_table;
        std::map<std::string, double> level_no_table;

        PARAM_TYPE parameter_id_type;

        std::unordered_set<std::string> known_params;

        // initialization
        void init_known_params();

        void readTable();
        bool parseIndexing(const std::vector<std::string> &line_vals);
        bool parseParamType(const std::vector<std::string> &line_vals);
        void fillParameterTables(const std::vector<std::string> &line_vals, std::vector<std::string>& unknown_params);

        static const std::string default_table;


};

// Needs to contain all variables for 1 timepoint, and only 1 timepoint
class GRIBFile {
    public:
        GRIBFile(const std::string &in_filename, const std::vector<std::string> &indexes);

        template <typename T>
        std::vector<CodesHandlePtr> listParameterMessages(const std::string& param_key, const T& paramID, const std::string& level_key,const std::string &levelType) {
            return getMessages(file, param_key, paramID, level_key, levelType);
        };

    private:
        std::string filename;
        CodesIndexPtr file;

        std::map<std::string, double> grid_params;

        bool checkValidity();

};

} // namespace mio

#endif // GRIBFILE_H
