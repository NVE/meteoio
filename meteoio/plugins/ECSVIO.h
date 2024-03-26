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
#ifndef ECSV_H
#define ECSV_H

#include <meteoio/FileUtils.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/plugins/ECSVHelper.h>
#include <meteoio/plugins/libacdd.h>

#include <vector>

namespace mio {

    using namespace ECSV;
    /**
     * @class ECSVIO
     * @brief A class to read and write ECSV files
     *
     * @section TODO (when available)
     * - if meteoparam metadata is available use it to pass along long_name ... and also write it out
     *     
     *
     * @ingroup plugins
     * @author Patrick Leibersperger
     * @date   2024-02-015
     */
    class ECSVIO : public IOInterface {
    public:
        ECSVIO(const std::string &configfile);
        ECSVIO(const ECSVIO &);
        ECSVIO(const Config &cfgreader);

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
        std::vector<ECSVFile> stations_files;

        // output section
        ACDD acdd_metadata;
        double TZ_out;
        std::string outpath;
        bool allow_overwrite;
        bool allow_append;
        char out_delimiter;

        // constants
        const std::string ECSV_version = "1.0";
        const std::string ECSV_firstline = "# ECSV " + ECSV_version + " UTF-8";
        static const double snVirtualSlopeAngle;

        // constructor helpers
        void parseInputSection();
        void parseOutputSection();

        // read helpers
        void identify_fields(const std::vector<std::string> &fields, std::vector<size_t> &indexes, MeteoData &md,
                                     const std::string &geometry_field);
        double getSnowpackSlope(const std::string &id);
        void read_meta_data(const ECSVFile &current_file, StationData &meta);
        void setMetaDataPosition(const ECSVFile &current_file, StationData &meta, const double &nodata_value);
        void setMetaDataSlope(const ECSVFile &current_file, StationData &meta, const double &nodata_value);

        void readDataSequential(ECSVFile &current_file);
        std::vector<MeteoData> createMeteoDataVector(ECSVFile &current_file, std::vector<Date> &date_vec,
                                                             std::vector<geoLocation> &location_vec);
        void setMeteoDataLocation(MeteoData &tmp_md, geoLocation &loc, ECSVFile &current_file, double nodata);
        void setMeteoDataFields(MeteoData &tmp_md, ECSVFile &current_file, Date &date, std::vector<size_t> &indexes, double nodata);

        // write helpers
        void prepareOutfile(ECSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists);
        void handleNewFile(ECSVFile &outfile, const std::vector<MeteoData> &vecMeteo, bool file_exists);
        void handleFileAppend(ECSVFile &outfile, const std::vector<MeteoData> &vecMeteo);

        std::string getGeometry(const geoLocation &loc);
        bool checkLocationConsistency(const std::vector<MeteoData> &vecMeteo);
        bool createFilename(ECSVFile &outfile, const StationData &station, size_t station_num);
        void createMetaDataSection(ECSVFile &current_file, const std::vector<MeteoData> &vecMeteo);
        void createFieldsSection(ECSVFile &current_file, const std::vector<MeteoData> &vecMeteo);
        void writeToFile(const ECSVFile &outfile);
    };

} // namespace
#endif
