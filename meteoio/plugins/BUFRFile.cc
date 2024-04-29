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

    using namespace codes;
    // initialize start_date very low, and end_date very high
    BUFRFile::BUFRFile(const std::string &in_filename, const std::string &ref_coords)
        : filename(in_filename), meta_data(), start_date(), end_date(), tz(), isMeteoTimeseries(false), isProfile(false), subset_numbers(), in_file(fopen(filename.c_str(), "r")), station_ids_in_file()  {
        if (!in_file) {
            throw AccessException("Could not open file " + filename);
        };

        readMetaData(ref_coords);
    };


    // ------------------ STATIC HELPERS ------------------
    static void getStationIdfromMessage(CodesHandlePtr &h, std::string &station_id, std::string &station_name, const size_t& subsetNumber) {
        std::vector<std::string> id_keys = {"stationNumber", "wigosIdentifierSeries", "shortStationName"};
        std::vector<std::string> name_keys = {"stationOrSiteName", "longStationName"};
        std::string stationID, stationName;
        getParameter(h, id_keys, stationID, subsetNumber);
        getParameter(h, name_keys, stationName, subsetNumber);

        if (stationID.empty()) {
            if (stationName.empty()) {
                throw AccessException("No station ID or name found in BUFR file, need to identify stations", AT);
            } else {
                stationID = stationName;
            }
        }
        station_id = stationID;
        station_name = stationName;
        return;
    };

    // Read all the necessary metadata from a BUFR file
    static StationData getStationDataBUFR(CodesHandlePtr &h, const std::string &ref_coords, std::string &error, const size_t &subsetNumber) {
        std::string subset_prefix = getSubsetPrefix(subsetNumber);
        StationData stationData;

        double latitude, longitude, altitude;
        getParameter(h, subset_prefix + "latitude", latitude);
        getParameter(h, subset_prefix + "longitude", longitude);

        std::vector<std::string> height_keys = {"heightOfStation", "height", "elevation"};
        getParameter(h, height_keys, altitude, subsetNumber);

        std::string stationID, stationName;
        getStationIdfromMessage(h, stationID, stationName, subsetNumber);

        long ref_flag = 999;
        getParameter(h, subset_prefix + "coordinateReferenceSystem", ref_flag);
        Coords position;

        if (ref_flag == 999 && ref_coords.empty()) {
            error = "No reference coordinates found in BUFR file and none provided in the configuration";
            return stationData;
        }

        if (ref_flag == 65535) {
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

    static void fillFromMessage(MeteoData &md, CodesHandlePtr &message, const std::vector<std::string> &additional_params, const size_t &subsetNumber) {
        std::string subset_prefix = getSubsetPrefix(subsetNumber);
        // fill the MeteoData object with data from the message and additional_params
        for (size_t id = 0; id < md.getNrOfParameters(); id++) {
            std::string param_name = md.getParameterName(id);
            getParameter(message, subset_prefix + BUFR_PARAMETER.at(param_name), md(id));
            if (param_name == "TAU_CLD") {
                md(id) = md(id) / 100.0;
            }
        }
        for (const auto &params : additional_params) {
            double value = IOUtils::nodata;
            getParameter(message, subset_prefix + params, value);
            size_t param_id = md.addParameter(params);
            md(param_id) = value;
        }
    };

    static bool areMultipleSubsetsInMessage(CodesHandlePtr &h, double &num_subsets) {
        getParameter(h, "numberOfSubsets", num_subsets);
        return num_subsets > 1;
    };

    // ------------------ GET ALL STATION DATA ------------------
    void BUFRFile::readMetaData(const std::string &ref_coords) {
        std::vector<CodesHandlePtr> messages = getMessages(in_file.get(), PRODUCT_BUFR);
        const double NO_SUBSETS_FOUND = -1;
        size_t message_counter = 1;
        std::map<std::string, std::set<Date>> station_dates;

        for (auto &message : messages) {
            double num_subsets = getNumSubsets(message, NO_SUBSETS_FOUND);
            subset_numbers.push_back(static_cast<size_t>(num_subsets));

            unpackMessage(message);

            processSubsets(message, ref_coords, subset_numbers.back(), station_dates);
            message_counter++;
        }
        if (meta_data.empty()) {
            throw IOException("No valid stations found in file " + filename);
        }
        tz = IOUtils::nodata; // TODO: find timezone from location
        isMeteoTimeseries = true;
    }

    // ------------------ STATION HELPERS ------------------
    double BUFRFile::getNumSubsets(CodesHandlePtr &message, double default_value) {
        bool hasMultipleSubsets = areMultipleSubsetsInMessage(message, default_value);
        return hasMultipleSubsets ? default_value : 1;
    }

    void BUFRFile::processSubsets(CodesHandlePtr &message, const std::string &ref_coords, const size_t& num_subsets, std::map<std::string, std::set<Date>> &station_dates) {
        std::string error_msg;
        for (size_t sb_id = 0; sb_id < num_subsets; sb_id++) {
            StationData new_meta = getStationDataBUFR(message, ref_coords, error_msg, sb_id);
            Date new_date = getMessageDateBUFR(message, sb_id, tz);
            std::string station_id = new_meta.getStationID();

            if (isNewStation(station_id)) {
                processNewStation(station_id, new_meta, new_date, station_dates);
            } else {
                processExistingStation(station_id, new_meta, new_date, station_dates);
            }

            updateDateRange(new_date);
        }
    }

    void BUFRFile::processNewStation(const std::string &station_id, const StationData &new_meta, const Date &new_date, std::map<std::string, std::set<Date>> &station_dates) {
        if (!new_meta.isValid())
            throw IOException("No valid station metadata for station " + station_id + " file " + filename);
        station_ids_in_file.push_back(station_id);
        meta_data[station_id] = new_meta;

        // create a set to be sure we done have duplicate timepoints
        std::set<Date> dates;
        dates.insert(new_date);
        station_dates[station_id] = dates;
    }

    void BUFRFile::processExistingStation(const std::string &station_id, const StationData &new_meta, const Date &new_date, std::map<std::string, std::set<Date>> &station_dates) {
        // check if the date is already in the set
        auto success = station_dates[station_id].insert(new_date);
        if (!success.second) {
            throw IOException("Duplicate date for station " + station_id + " in file " + filename);
        }

        // check consistency between station data
        if (!new_meta.isEmpty() && meta_data[station_id] != new_meta) {
            throw IOException("Inconsistent station data for station " + station_id + " in file " + filename);
        }
    }

    void BUFRFile::updateDateRange(const Date &new_date) {
        if (new_date < start_date) {
            start_date = new_date;
        };
        if (new_date > end_date) {
            end_date = new_date;
        };
    }

    // ------------------ READ DATA ------------------
    void BUFRFile::readData(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping, const std::vector<std::string> &additional_params) {
        // a bufr file can contain multiple messages, each message can contain multiple subsets
        // the order of where stations are stored at what date is quite arbitrary, so we need to map the station ids to the correct index in vecMeteo
        // as well as keeping track between calls, what station ids we already have
        auto messages = getMessages(in_file.get(), PRODUCT_BUFR);
        updateStationIdsMapping(vecMeteo, station_ids_mapping);

        for (size_t msg_id = 0; msg_id < messages.size(); msg_id++) {
            processMessage(vecMeteo, station_ids_mapping, additional_params, messages[msg_id], msg_id);
        }
    }

    // ------------------ HELPERS ------------------
    void BUFRFile::updateStationIdsMapping(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping) {
        for (const auto& station : station_ids_in_file) {
            if (station_ids_mapping.find(station) == station_ids_mapping.end()) {
                METEO_SET empty_set;
                station_ids_mapping[station] = vecMeteo.size();
                vecMeteo.push_back(empty_set);
            }
        }
    }

    void BUFRFile::processMessage(std::vector<METEO_SET> &vecMeteo, std::map<std::string, size_t>& station_ids_mapping, const std::vector<std::string> &additional_params, CodesHandlePtr &message, size_t msg_id) {
        // TODO: if we need speedup cache the date and meta data (station id) for each message and subset
        for (size_t sub_id = 0; sub_id < subset_numbers[msg_id]; sub_id++) {
            std::string station_id, stationName;
            getStationIdfromMessage(message, station_id, stationName, sub_id);
            size_t index_in_vecMeteo = station_ids_mapping[station_id];
            
            Date date = getMessageDateBUFR(message, sub_id,tz);
            MeteoData md;
            md.setDate(date);
            md.meta = meta_data[station_id];
            fillFromMessage(md, message, additional_params, sub_id);
            vecMeteo[index_in_vecMeteo].push_back(md);
        }
    }
} // namespace mio