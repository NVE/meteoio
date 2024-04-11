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

#ifndef BUFRFILE_H
#define BUFRFILE_H

#include <meteoio/IOInterface.h>
#include <meteoio/plugins/libcodes.h>

#include <string>

namespace mio {

using namespace codes;

class BUFRFile {
    public:
        BUFRFile(const std::string &filename);
        
        void readMetaData();
        void readData(std::vector<MeteoData> &vecMeteo, const std::vector<std::string> &additional_params = {});

    private:
        std::string filename;
        StationData meta_data;
        Date start_date;
        Date end_date;

        std::unique_ptr<FILE> in_file;
        std::vector<CodesHandlePtr> messages;

};

} //namespace


#endif // BUFRFILE_H