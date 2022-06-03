#ifndef REQUESTHANDLER_H
#define REQUESTHANDLER_H

#include "oatpp/web/server/HttpRequestHandler.hpp"
#include "rapidxml_ns-1.13.2/rapidxml_ns.hpp"
#include <iostream>
#include <fstream>
#include "BadRequestException.h"
#include "InternalServerError.h"
#include "UUID.h"
#include <iomanip>
#include <ctime>
#include "Timeseries.h"
#include <sys/stat.h>

using namespace std;

// Custom request handler
class RequestHandler : public oatpp::web::server::HttpRequestHandler
{
public:

    RequestHandler(unsigned int default_timeout_secs, string job_directory)
        : default_timeout_secs(default_timeout_secs), job_directory(job_directory) 
    {}

    // Process incoming requests and return responses
    shared_ptr<OutgoingResponse> handle(const shared_ptr<IncomingRequest> &request) override
    {
        string body = request->readBodyToString();
        rapidxml_ns::xml_document<> doc;

        vector<char> body_copy(body.begin(), body.end());
        body_copy.push_back('\0');

        doc.parse<0>(&body_copy[0]);

        rapidxml_ns::xml_node<> *root_node = doc.first_node();
        string operationName = root_node->local_name();
        string operationNameNs = root_node->namespace_uri();

        if (operationName == OPERATION_EXECUTION && operationNameNs == NS_WPS)
        {
            try
            {
                ExecuteRequest executeRequest = readExecuteRequest(root_node);                
                Timeseries timeseries = timeseriesFromExecuteRequest(executeRequest);
                timeseries.run();

                // TODO
                // Return output

                string responseBody =
                    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                    "<wps:Result"
                    "	xmlns:wps=\"http://www.opengis.net/wps/2.0\""
                    "	xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
                    "	"
                    "	<wps:JobID>" + executeRequest.jobId + "</wps:JobID>"
                    // "	<wps:ExpirationDate>executeRequest.expirationDate</wps:ExpirationDate>"
                    "	<wps:Output id=\"RESULT\">"
                    "       <wps:Data>10.0</wps:Data>"
                    "	</wps:Output>"
                    "	"
                    "</wps:Result>";
                shared_ptr<OutgoingResponse> response = ResponseFactory::createResponse(Status::CODE_200, responseBody);
                response->putHeader("Content-Type", "text/xml");
                return response;
            }
            catch (BadRequestException &e)
            {
                return ResponseFactory::createResponse(Status::CODE_400, e.what());
            }
            catch (exception &e)
            {
                return ResponseFactory::createResponse(Status::CODE_500, e.what());
            }
        }

        return ResponseFactory::createResponse(Status::CODE_400, "Operation " + operationName + " is not supported");
    }

private:
    const char *NS_WPS = "http://www.opengis.net/wps/2.0";
    const char *OPERATION_EXECUTION = "Execute";
    const char *EXECUTION_MODE = "mode";
    const char *EXECUTION_INPUT = "Input";
    const char *EXECUTION_INPUT_ID = "id";
    const char *EXECUTION_INPUT_DATA = "Data";
    const char *EXECUTION_INPUT_DATA_ENCODING = "encoding";
    const char *EXECUTION_INPUT_DATA_MIME_TYPE = "mimeType";

    unsigned int default_timeout_secs = 60;
    string job_directory = "/tmp/jobs";

    struct ExecutionInput
    {
        string id;
        string rawData;
        string encoding = "UTF-8";
        string mimeType = "text/plain";
    };

    struct ExecuteRequest
    {
        string jobId = UUID::generate();
        // string expirationDate = put_time(&std::localtime(&time(nullptr)), "%d-%m-%Y %H-%M-%S");
        string mode;
        vector<ExecutionInput> inputs;
        string output;
    };

    inline ExecuteRequest readExecuteRequest(rapidxml_ns::xml_node<> *root_node)
    {
        ExecuteRequest executeRequest;
        executeRequest.mode = root_node->first_attribute(EXECUTION_MODE)->value();
        if(executeRequest.mode != "sync" && executeRequest.mode != "auto")
            throw BadRequestException("Only execution mode 'sync' is currently supported");

        for (rapidxml_ns::xml_node<> *input = root_node->first_node_ns(NS_WPS, EXECUTION_INPUT); input; input = input->next_sibling_ns(NS_WPS, EXECUTION_INPUT))
        {
            rapidxml_ns::xml_node<> *data = input->first_node_ns(NS_WPS, EXECUTION_INPUT_DATA);
            if (data)
            {
                ExecutionInput executionInput;
                executionInput.id = input->first_attribute(EXECUTION_INPUT_ID)->value();

                executionInput.rawData = data->value();
                if (data->first_attribute(EXECUTION_INPUT_DATA_ENCODING))
                    executionInput.encoding = data->first_attribute(EXECUTION_INPUT_DATA_ENCODING)->value();
                if (data->first_attribute(EXECUTION_INPUT_DATA_MIME_TYPE))
                    executionInput.encoding = data->first_attribute(EXECUTION_INPUT_DATA_MIME_TYPE)->value();
                executeRequest.inputs.push_back(executionInput);
            }
            else
            {
                throw BadRequestException("Only input of type '" + string(EXECUTION_INPUT_DATA) + "' is supported");
            }
        }
        return executeRequest;
    }

    inline void checkEncoding(ExecutionInput &input) {
        if(input.encoding != "UTF-8") {
            throw BadRequestException("Currently, only input with encoding 'UTF-8' is supported");
        }
    }

    inline string createInputFile(string &workingDir, ExecutionInput &input) {
        string filepath = workingDir + "/" + input.id;
        ofstream out(filepath.c_str());
        out << input.rawData;
        out.close();
        return filepath;
    }

    inline Date getDate(const std::string& date_str, const double& TZ)
    {
        Date parsedDate;
        if (date_str == "NOW") { //interpret user provided start date
            parsedDate.setFromSys();
            parsedDate.setTimeZone(TZ);
            parsedDate.rnd(10, mio::Date::DOWN); //rounding 10' down
        } else {
            mio::IOUtils::convertString(parsedDate, date_str, TZ);
        }
        
        return parsedDate;
    }

    inline string createWorkingDir(string &jobId) {
        string workingDir = job_directory + "/" + jobId;
        const int dir_err = mkdir(workingDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
        {
            throw InternalServerError("Error creating directory!");
        }
        return workingDir;
    }

    Timeseries timeseriesFromExecuteRequest(ExecuteRequest &executeRequest)
    {
        Config cfg;
        Date dateBegin, dateEnd;
        double samplingRate = IOUtils::nodata;
	    size_t outputBufferSize = 0;
	    unsigned int timeoutSecs = default_timeout_secs;
        string workingDir = createWorkingDir(executeRequest.jobId);

        std::string cfgfile = "";
        std::string begin_date_str, end_date_str;
        double duration = IOUtils::nodata;
        int longindex=0, opt=-1;
        bool setStart = false, setEnd = false, setDuration = false;
        
        for (ExecutionInput &input: executeRequest.inputs) {
            checkEncoding(input);

            if(input.id == "cfg.ini") {
                cfgfile = createInputFile(workingDir, input);
            }
            else if(input.id == "begin_date") {
                begin_date_str = input.rawData;
                setStart = true;
            }
            else if(input.id == "end_date") {
                end_date_str = input.rawData;
                setEnd = true;
            }
            else if(input.id == "duration") {
                mio::IOUtils::convertString(duration, input.rawData);
                setDuration = true;
            }
            else if(input.id == "sampling_rate") {
                mio::IOUtils::convertString(samplingRate, input.rawData);
            }
            else if(input.id == "output_buffer_size") {
                mio::IOUtils::convertString(outputBufferSize, input.rawData);
            }
            else if(input.id == "timeout_secs") {
                mio::IOUtils::convertString(timeoutSecs, input.rawData);
            }
            else {
                createInputFile(workingDir, input);
            }            
        }

        if(cfgfile == "") {
            throw BadRequestException("You must provide a 'cfg' ini file!");
        }

        const bool validDateRange = (setStart && setEnd && !setDuration) || (setStart && !setEnd && setDuration) || (!setStart && setEnd && setDuration);
        if (!validDateRange) {
            throw BadRequestException("You must specify either {'begin_date' and 'end_date'}, or {'begin_date' and 'duration'} or {'end_date' and 'duration'}!");
        }
        
        cfg.addFile(cfgfile);
        const double TZ = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
        
        //the date range specification has been validated above
        if (!begin_date_str.empty()) dateBegin = getDate( begin_date_str, TZ );
        if (!end_date_str.empty()) 
            dateEnd = getDate( end_date_str, TZ );
        else
            dateEnd = dateBegin + duration;
        if (dateBegin.isUndef())
            dateBegin = dateEnd - duration;
        
        //we don't overwrite command line options if set
        if (samplingRate==IOUtils::nodata)
            samplingRate = cfg.get("SAMPLING_RATE_MIN", "Output", 60.);
        samplingRate /= 24.*60; //convert to sampling rate in days

        Timeseries timeseries(cfg, dateBegin, dateEnd);
        timeseries.setSamplingRate(samplingRate);
        timeseries.setTimeoutSecs(timeoutSecs);
        timeseries.setOutputBufferSize(outputBufferSize);
        return timeseries;
    }

};

#endif // REQUESTHANDLER_H