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

#include <meteoio/IOInterface.h>

#include <eccodes.h>
#include <vector>

#include <iostream>

namespace mio
{   
    namespace codes {
        
        struct HandleDeleter {
            void operator()(codes_handle* h) const {
                codes_handle_delete(h);
            }
        };

        struct IndexDeleter {
            void operator()(codes_index* i) const {
                codes_index_delete(i);
            }
        };
        using CodesHandlePtr = std::unique_ptr<codes_handle, HandleDeleter>;
        using CodesIndexPtr = std::unique_ptr<codes_index, IndexDeleter>;
        CodesHandlePtr makeUnique(codes_handle* h) {
            CodesHandlePtr ptr(h);
            h = nullptr;
            return ptr;
        }
        CodesIndexPtr makeUnique(codes_index* i) {
            CodesIndexPtr ptr(i);
            i = nullptr;
            return ptr;
        }


        CodesIndexPtr indexFile(const std::string &filename, std::vector<std::string> &paramIdList, std::vector<long> &ensembleNumbers, std::vector<long> &levelTypes);
        std::vector<CodesHandlePtr> getMessages(CodesIndexPtr &index, const std::string &paramID, const long& ensembleNumber, const long& levelType);
        std::vector<CodesHandlePtr> getMessages(const std::string &filename);

        Date getMessageDate(CodesHandlePtr &h, double &d1, double &d2, const double &tz_in);

        std::map<std::string, double> getGridParameters(CodesHandlePtr &h_unique );
        void getGriddedValues(CodesHandlePtr &h, std::vector<double> &values, std::map<std::string,double> &gridParams);
        std::map<std::string,double> getGriddedValues(CodesHandlePtr &h, std::vector<double> &values);

        std::tuple<std::vector<double>&&, std::vector<double>&&, std::vector<double>&&, std::vector<double>&&, std::vector<int>&&> getNearestValues_grib(CodesHandlePtr &h, const std::vector<double> &in_lats, const std::vector<double> &in_lons);

        void getParameter(CodesHandlePtr &h, const std::string& parameterName, double &param_value);
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, long &param_value);
        void getParameter(CodesHandlePtr &h, const std::string& parameterName, std::string& param_value);

    }
} // namespace mio
