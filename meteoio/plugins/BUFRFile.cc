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

/**
 * @todo Maybe we cannot save all messages in memory, but need to read them one by one
 * 
 * @note For reading profiles, it will be possible to read all iterations in a subset (only way to allow subsets)
**/

#include <meteoio/plugins/BUFRFile.h>

namespace mio {


// initialize start_date very low, and end_date very high
BUFRFile::BUFRFile(const std::string &filename, const std::string& ref_coords) : filename(filename) , meta_data(), start_date(), end_date(), dates(), tz(), in_file(fopen(filename.c_str(), "r")), messages() {
    if (!in_file) {
        throw AccessException("Could not open file " + filename);
    };

    readMetaData(ref_coords);
};

void BUFRFile::readMetaData(const std::string& ref_coords) {
    messages = getMessages(in_file.get(), PRODUCT_BUFR);
    StationData meta_data;
    std::string error_msg;
    for (CodesHandlePtr &message : messages) {    
        long num_subsets = -1;
        getParameter(message, "numberOfSubsets", num_subsets);
        if (num_subsets > 1) {
            throw IOException("Subsets in BUFR file " + filename + " are not supported", AT);
        };
        unpackMessage(message);
        Date date = getMessageDateBUFR(message);
        if (date < start_date) {
            start_date = date;
        };
        if (date > end_date) {
            end_date = date;
        };
        auto success = dates.insert(date);
        if (success.second == false) {
            throw IOException("Duplicate date in file " + filename);
        };

        StationData new_meta = getStationDataBUFR(message, ref_coords, error_msg);
        if (!meta_data.isValid()) {
            meta_data = new_meta;
            continue;
        };

        if (!new_meta.isEmpty() && meta_data != new_meta) {
            throw IOException("Inconsistent station data in file " + filename);
        }
    };

    if (!meta_data.isValid()) {
        std::ostringstream msg;
        msg << "No valid station data in file " << filename;
        if (!error_msg.empty()) {
            msg << "\n(" << error_msg << ")";
        }
        throw IOException(msg.str(),AT);
    }
    tz = IOUtils::nodata; // TODO: find timezone from location
    return;
};


static void fillFromMessage(MeteoData &md, CodesHandlePtr &message, const std::map<std::string,std::string> &additional_params) {
    // fill the MeteoData object with data from the message and additional_params
    for (size_t id = 0; id < md.getNrOfParameters(); id++) {
        std::string param_name = md.getParameterName(id);
        getParameter(message, BUFR_PARAMETER.at(param_name), md(id));
        if (param_name == "TAU_CLD") {
            md(id) = md(id) / 100.0;
        }
    }
    for (const auto &params : additional_params) {
        std::string param_key = params.first;
        std::string param_name = params.second;
        double value = IOUtils::nodata;
        getParameter(message, param_key, value);
        size_t param_id = md.addParameter(param_name);
        md(param_id) = value;
    }
};

void BUFRFile::readData(std::vector<MeteoData> &vecMeteo, const std::map<std::string,std::string> &additional_params) {
    // Read the data
    // try for all meteo standard parameters and additional_params
    vecMeteo.clear();   // TODO: there is no possibility of vecMeteo being non-empty, correct?
    vecMeteo.resize(messages.size());
    for (CodesHandlePtr &message : messages) {
        Date date = getMessageDateBUFR(message, tz); // TODO: can use Dates set, but if there is subsets there might be more dates than subsets
        MeteoData md;
        md.setDate(date);
        md.meta = meta_data;
        fillFromMessage(md, message, additional_params);
        vecMeteo.push_back(md);
    };
};

} //namespace mio