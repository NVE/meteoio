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

#include "BUFRFile.h"

namespace mio {


// initialize start_date very low, and end_date very high
BUFRFile::BUFRFile(const std::string &filename) : filename(filename) , meta_data(), start_date(), end_date(), in_file(fopen(filename.c_str(), "r")), messages() {
    if (!in_file) {
        throw AccessException("Could not open file " + filename);
    };

    readMetaData();
};

void BUFRFile::readMetaData() {
    messages = getMessages(in_file.get(), PRODUCT_BUFR);
    StationData meta_data;
    for (CodesHandlePtr &message : messages) {
        if (meta_data.getStationID().empty()) {
            // fill metadata
        } else {
            // check if metadata is consistent
        };
        
        Date date = getMessageDateBUFR(message);
        if (date < start_date) {
            start_date = date;
        };
        if (date > end_date) {
            end_date = date;
        };

    };
};

void BUFRFile::readData(std::vector<MeteoData> &vecMeteo, const std::vector<std::string> &additional_params) {
    // Read the data
    // try for all meteo standard parameters and additional_params
    // ...
};

} //namespace mio