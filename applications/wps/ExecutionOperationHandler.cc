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

#include "ExecutionOperationHandler.h"
#include "util/Zip.h"

ExecutionOperationHandler::ExecutionOperationHandler(string job_directory) : _job_directory(job_directory)
{
}

string ExecutionOperationHandler::handleOperation(rapidxml_ns::xml_node<> *root_node)
{
    ExecuteRequest executeRequest = readExecuteRequest(root_node);
    OATPP_LOGI("ExecutionOperationHandler", "Processing execution operation for jobId='%s'", executeRequest.jobId.c_str());

    string workingDir = getWorkingDirectory(executeRequest);
    string resultDir = getResultDirectory(executeRequest);
    createDir(workingDir);
    createDir(resultDir);
    shared_ptr<Timeseries> timeseries = timeseriesFromExecuteRequest(executeRequest, workingDir, resultDir);
    timeseries->run();

    Zip::zipDirectory(resultDir);

    // Return output
    return getResponseBody(executeRequest);
}

ExecutionOperationHandler::ExecuteRequest ExecutionOperationHandler::readExecuteRequest(rapidxml_ns::xml_node<> *root_node)
{
    ExecuteRequest executeRequest;
    executeRequest.mode = root_node->first_attribute(EXECUTION_MODE)->value();
    if (executeRequest.mode != "sync" && executeRequest.mode != "auto")
        throw BadRequestException("Only execution mode 'sync' is currently supported");

    for (rapidxml_ns::xml_node<> *input = root_node->first_node_ns(NS_WPS, EXECUTION_INPUT); input; input = input->next_sibling_ns(NS_WPS, EXECUTION_INPUT))
    {
        executeRequest.inputs.push_back(readExecutionInput(input));
    }
    return executeRequest;
}

ExecutionOperationHandler::ExecutionInput ExecutionOperationHandler::readExecutionInput(rapidxml_ns::xml_node<> *input)
{
    rapidxml_ns::xml_node<> *data = input->first_node_ns(NS_WPS, EXECUTION_INPUT_DATA);
    if (!data)
        throw BadRequestException("Only input of type '" + string(EXECUTION_INPUT_DATA) + "' is supported");

    ExecutionInput executionInput;
    executionInput.id = input->first_attribute(EXECUTION_INPUT_ID)->value();

    executionInput.rawData = data->value();
    if (data->first_attribute(EXECUTION_INPUT_DATA_ENCODING))
        executionInput.encoding = data->first_attribute(EXECUTION_INPUT_DATA_ENCODING)->value();
    if (data->first_attribute(EXECUTION_INPUT_DATA_MIME_TYPE))
        executionInput.encoding = data->first_attribute(EXECUTION_INPUT_DATA_MIME_TYPE)->value();
    return executionInput;
}

shared_ptr<Timeseries> ExecutionOperationHandler::timeseriesFromExecuteRequest(ExecuteRequest &executeRequest, string &workingDir, string &resultDir)
{
    double samplingRate = IOUtils::nodata;
    size_t outputBufferSize = 0;

    string cfgfile = "";
    string begin_date_str, end_date_str;
    double duration = IOUtils::nodata;

    for (ExecutionInput &input : executeRequest.inputs)
    {
        checkEncoding(input);

        if (input.id == "cfg.ini")
        {
            cfgfile = createInputFile(workingDir, input);
        }
        else if (input.id == "begin_date")
        {
            begin_date_str = input.rawData;
        }
        else if (input.id == "end_date")
        {
            end_date_str = input.rawData;
        }
        else if (input.id == "duration")
        {
            mio::IOUtils::convertString(duration, input.rawData);
        }
        else if (input.id == "sampling_rate")
        {
            mio::IOUtils::convertString(samplingRate, input.rawData);
        }
        else if (input.id == "output_buffer_size")
        {
            mio::IOUtils::convertString(outputBufferSize, input.rawData);
        }
        else
        {
            createInputFile(workingDir, input);
        }
    }

    if (cfgfile == "")
        throw BadRequestException("You must provide a 'cfg.ini' file!");

    Config cfg = createConfig(workingDir, resultDir, cfgfile);

    DateRange dateRange = getDateRange(begin_date_str, end_date_str, duration, cfg);

    // we don't overwrite request options if set
    if (samplingRate == IOUtils::nodata)
        samplingRate = cfg.get("SAMPLING_RATE_MIN", "Output", 60.);
    samplingRate /= 24. * 60; // convert to sampling rate in days

    shared_ptr<Timeseries> timeseries = std::make_shared<Timeseries>(cfg, dateRange.dateBegin, dateRange.dateEnd);
    timeseries->setSamplingRate(samplingRate);
    timeseries->setOutputBufferSize(outputBufferSize);
    return timeseries;
}