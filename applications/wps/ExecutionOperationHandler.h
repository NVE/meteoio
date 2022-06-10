// SPDX-License-Identifier: LGPL-3.0-or-later
/*
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
 */
/*  This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EXECUTION_OPERATION_HANDLER_H
#define EXECUTION_OPERATION_HANDLER_H

#include <fstream>
#include <sys/stat.h>
#include "oatpp/core/base/Environment.hpp"
#include "WpsOperationHandler.h"
#include "util/UUID.h"
#include "Timeseries.h"
#include "exceptions/BadRequestException.h"
#include "exceptions/InternalServerError.h"

using namespace std;

class ExecutionOperationHandler : public WpsOperationHandler
{
public:
    ExecutionOperationHandler(string job_directory);

    string handleOperation(rapidxml_ns::xml_node<> *root_node) override;

private:
    const char *NS_WPS = "http://www.opengis.net/wps/2.0";
    const char *OPERATION_EXECUTION = "Execute";
    const char *EXECUTION_MODE = "mode";
    const char *EXECUTION_INPUT = "Input";
    const char *EXECUTION_INPUT_ID = "id";
    const char *EXECUTION_INPUT_DATA = "Data";
    const char *EXECUTION_INPUT_DATA_ENCODING = "encoding";
    const char *EXECUTION_INPUT_DATA_MIME_TYPE = "mimeType";

    string _job_directory;

    struct ExecutionInput
    {
        string id = "";
        string rawData = "";
        string encoding = "UTF-8";
        string mimeType = "text/plain";
    };

    struct ExecuteRequest
    {
        string jobId = UUID::generate();
        string expirationDate = ""; // put_time(&std::localtime(&time(nullptr)), "%d-%m-%Y %H-%M-%S");
        string mode = "";
        vector<ExecutionInput> inputs = {};
        string output = "";
    };

    struct DateRange
    {
        Date dateBegin;
        Date dateEnd;
    };

    ExecuteRequest readExecuteRequest(rapidxml_ns::xml_node<> *root_node);

    ExecutionInput readExecutionInput(rapidxml_ns::xml_node<> *input);

    shared_ptr<Timeseries> timeseriesFromExecuteRequest(ExecuteRequest &executeRequest, string &workingDir, string &resultDir);

    inline string getResponseBody(ExecuteRequest &executeRequest)
    {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
               "<wps:Result"
               "	xmlns:wps=\"http://www.opengis.net/wps/2.0\""
               "	xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
               "	<wps:JobID>" + executeRequest.jobId + "</wps:JobID>" +
               (executeRequest.expirationDate.empty() ? "" : 
               "	<wps:ExpirationDate>" + executeRequest.expirationDate + "</wps:ExpirationDate>"
               ) +
               "	<wps:Output id=\"RESULT\">"
               "       <wps:Reference xlink:href=\"/results/" + executeRequest.jobId + "/result.zip\"/>"
               "	</wps:Output>"
               "</wps:Result>";
    }

    inline void checkEncoding(ExecutionInput &input)
    {
        if (input.encoding != "UTF-8")
        {
            throw BadRequestException("Currently, only input with encoding 'UTF-8' is supported");
        }
    }

    inline string createInputFile(string &workingDir, ExecutionInput &input)
    {
        string filepath = workingDir + "/" + input.id;
        ofstream out(filepath.c_str());
        out << input.rawData;
        out.close();
        return filepath;
    }

    inline Date getDate(const std::string &date_str, const double &TZ)
    {
        Date parsedDate;
        if (date_str == "NOW")
        { // interpret user provided start date
            parsedDate.setFromSys();
            parsedDate.setTimeZone(TZ);
            parsedDate.rnd(10, mio::Date::DOWN); // rounding 10' down
        }
        else
        {
            mio::IOUtils::convertString(parsedDate, date_str, TZ);
        }

        return parsedDate;
    }

    inline DateRange getDateRange(string &begin_date_str, string &end_date_str, double &duration, Config &cfg)
    {
        bool setStart = !begin_date_str.empty();
        bool setEnd = !end_date_str.empty();
        bool setDuration = duration != IOUtils::nodata;
        const bool validDateRange = (setStart && setEnd && !setDuration) || (setStart && !setEnd && setDuration) || (!setStart && setEnd && setDuration);
        if (!validDateRange)
            throw BadRequestException("You must specify either {'begin_date' and 'end_date'}, or {'begin_date' and 'duration'} or {'end_date' and 'duration'}!");

        DateRange dateRange;
        const double TZ = cfg.get("TIME_ZONE", "Input"); // get user provided input time_zone
        // the date range specification has been validated above
        if (!begin_date_str.empty())
            dateRange.dateBegin = getDate(begin_date_str, TZ);
        if (!end_date_str.empty())
            dateRange.dateEnd = getDate(end_date_str, TZ);
        else
            dateRange.dateEnd = dateRange.dateBegin + duration;
        if (dateRange.dateBegin.isUndef())
            dateRange.dateBegin = dateRange.dateEnd - duration;
        return dateRange;
    }

    inline void createDir(string &dir)
    {
        OATPP_LOGI("ExecutionOperationHandler", "Creating dir: %s", dir.c_str());
        const int dir_err = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
        {
            throw InternalServerError("Error creating directory!");
        }
    }

    inline string getWorkingDirectory(ExecuteRequest &executeRequest)
    {
        return _job_directory + "/" + executeRequest.jobId;
    }

    inline string getResultDirectory(ExecuteRequest &executeRequest)
    {
        return _job_directory + "/" + executeRequest.jobId + "/result";
    }

    inline Config createConfig(string &workingDir, string &resultDir, string &cfgfile)
    {
        Config cfg;
        // add working dir so users can use it as a reference -> ${INPUT::CWD}
        cfg.addKey("CWD", "INPUT", workingDir);

        cfg.addFile(cfgfile);

        if (cfg.keyExists("METEOPATH", "INPUT"))
            cfg.addKey("METEOPATH", "INPUT", workingDir);
        if (cfg.keyExists("GRID2DPATH", "INPUT"))
            cfg.addKey("GRID2DPATH", "INPUT", workingDir);

        if (cfg.keyExists("METEOPATH", "OUTPUT"))
            cfg.addKey("METEOPATH", "OUTPUT", resultDir);
        if (cfg.keyExists("GRID2DPATH", "OUTPUT"))
            cfg.addKey("GRID2DPATH", "OUTPUT", resultDir);
        return cfg;
    }
};

#endif // EXECUTION_OPERATION_HANDLER_H