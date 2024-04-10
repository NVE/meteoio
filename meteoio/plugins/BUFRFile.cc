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

BUFRFile::BUFRFile(const std::string &filename) : filename(filename) , meta_data(), start_date(), end_date() {
    in_file = std::unique_ptr<FILE>(fopen(filename.c_str(), "r"));
    if (!in_file) {
        throw AccessException("Could not open file " + filename);
    };

    readMetaData();
};

void BUFRFile::readMetaData() {
    // Read the metadata
    // ...
};

void BUFRFile::readData(std::vector<MeteoData> &vecMeteo, const std::vector<std::string> &additional_params) {
    // Read the data
    // try for all meteo standard parameters and additional_params
    // ...
};

} //namespace mio