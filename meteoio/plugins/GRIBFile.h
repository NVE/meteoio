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

namespace mio
{

using namespace codes;

class GRIBTable {
    enum class KEY_TYPE {
        PARAM,
        LEVEL,
        ENSEMBLE
    };

    public:
        GRIBTable(const std::string &filename);

        std::string getIndexString() const;
        std::string getIndexingKey(KEY_TYPE &key) const;

        std::string getParamId(const std::string &param_name) const;
        std::string getLevelType(const std::string &level_name) const;
        double getLevelNo(const std::string &level) const;
    
    private:
        std::string filename;

        std::string param_indexing;
        std::string level_indexing;
        std::string ensemble_indexing;

        std::map<std::string, std::string> param_table;
        std::map<std::string, double> param_table_double;
        std::map<std::string, std::string> level_type_table;
        std::map<std::string, double> level_no_table;

        bool ensemble_average;
        size_t ensemble_no;

};

class GRIBFile {
    public:
        GRIBFile(const std::string &filename, const GRIBTable &table);

        std::vector<CodesHandlePtr> listParameterMessages(const std::string& param_name);
    
    private:
        std::string filename;
        GRIBTable table;
        CodesIndexPtr file;

        std::map<std::string, double> grid_params;

        bool checkValidity();

};

} // namespace mio

#endif // GRIBFILE_H
