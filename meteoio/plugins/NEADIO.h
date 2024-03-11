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
#ifndef NEAD_H
#define NEAD_H

#include <meteoio/FileUtils.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/plugins/NEADHelper.h>
#include <meteoio/plugins/libacdd.h>

#include <vector>

namespace mio {

    using namespace NEAD;
    /**
     * @class NEADIO
     * @brief A class to read and write NEAD files
     *
     * @section TODO (when available)
     * - file extension in inishell and plugin is .nead at the moment, is probably something else
     * - if meteoparam metadata is available use it to pass along long_name ... and also write it out
     *     
     * @section notes Notes on Nead
     * - geometry cant handle UTM
     * - timestamp should be iso and required for stations
     * - in pyNEAD the file format is forced to be UTF-8
     * - geometry = column name (is something like "lat,lon,alt" allowed?) --> fix TODOs if so
     * - should the file extension actually be .csv?
     * - xy + latlon in geometry
     * - the example files do not follow the format with add_offset and scale_factor...
     *
     * @ingroup plugins
     * @author Patrick Leibersperger
     * @date   2024-02-015
     */
    class NEADIO : public IOInterface {
    public:
        NEADIO(const std::string &configfile);
        NEADIO(const NEADIO &);
        NEADIO(const Config &cfgreader);

        virtual void readStationData(const Date &date, std::vector<StationData> &vecStation);
        virtual void readMeteoData(const Date &dateStart, const Date &dateEnd, std::vector<std::vector<MeteoData>> &vecMeteo);

        virtual void writeMeteoData(const std::vector<std::vector<MeteoData>> &vecMeteo, const std::string &name = "");

    private:
        // input section
        const Config cfg;
        std::string coordin, coordinparam, coordout, coordoutparam; // projection parameters
        bool snowpack_slopes;
        bool read_sequential;

        // file information
        std::vector<NEADFile> stations_files;

        // output section
        ACDD acdd_metadata;
        double TZ_out;
        std::string outpath;
        bool allow_overwrite;
        bool allow_append;
        char out_delimiter;

        // constants
        const std::string NEAD_version = "1.0";
        const std::string NEAD_firstline = "# NEAD " + NEAD_version + " UTF-8";
        static const double snVirtualSlopeAngle;

        // constructor helpers
        void parseInputSection();
        void parseOutputSection();

        // read helpers
        void identify_fields(const std::vector<std::string> &fields, std::vector<size_t> &indexes, MeteoData &md,
                                     const std::string &geometry_field);
        double getSnowpackSlope(const std::string &id);
        void read_meta_data(const NEADFile &current_file, StationData &meta);
        void setMetaDataPosition(const NEADFile &current_file, StationData &meta, const double &nodata_value);
        void setMetaDataSlope(const NEADFile &current_file, StationData &meta, const double &nodata_value);

        void readDataSequential(NEADFile &current_file);
        std::vector<MeteoData> createMeteoDataVector(NEADFile &current_file, std::vector<Date> &date_vec,
                                                             std::vector<geoLocation> &location_vec);
        void setMeteoDataLocation(MeteoData &tmp_md, geoLocation &loc, NEADFile &current_file, double nodata);
        void setMeteoDataFields(MeteoData &tmp_md, NEADFile &current_file, Date &date, std::vector<size_t> &indexes, double nodata);

        // write helpers
        void prepareOutfile(NEADFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists);
        void handleNewFile(NEADFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists);
        void handleFileAppend(NEADFile &outfile, const std::vector<MeteoData> &vecMeteo);

        std::string getGeometry(const geoLocation &loc);
        bool checkLocationConsistency(const std::vector<MeteoData> &vecMeteo);
        bool createFilename(NEADFile &outfile, const StationData &station, size_t station_num);
        void createMetaDataSection(NEADFile &current_file, const std::vector<MeteoData> &vecMeteo);
        void createFieldsSection(NEADFile &current_file, const std::vector<MeteoData> &vecMeteo);
        void writeToFile(const NEADFile &outfile);
    };

} // namespace
#endif
