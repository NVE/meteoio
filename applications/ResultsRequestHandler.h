#ifndef RESULTS_REQUESTHANDLER_H
#define RESULTS_REQUESTHANDLER_H

#include "oatpp/web/server/HttpRequestHandler.hpp"
#include "oatpp/core/data/stream/FileStream.hpp"
#include "oatpp/web/protocol/http/outgoing/StreamingBody.hpp"
#include <exception>

using namespace std;

// Custom request handler
class ResultsRequestHandler : public oatpp::web::server::HttpRequestHandler
{
public:
    ResultsRequestHandler(string job_directory) : _job_directory(job_directory) 
    {}

    // Process incoming requests and return responses
    shared_ptr<OutgoingResponse> handle(const shared_ptr<IncomingRequest> &request) override
    {
        string jobId = request->getPathVariable("jobId");
        string filepath = _job_directory + "/" + jobId + "/result.zip";
        try {
            auto body = std::make_shared<oatpp::web::protocol::http::outgoing::StreamingBody>(
                std::make_shared<oatpp::data::stream::FileInputStream>(filepath.c_str())
            );
            shared_ptr<OutgoingResponse> response =  OutgoingResponse::createShared(Status::CODE_200, body);
            response->putHeader("Content-Type", "application/zip");
            return response;
        }
        catch (exception &e)
        {
            cout << "Requested file '" << filepath << "' could not be served." << endl;
            cout << e.what() << endl;
            return ResponseFactory::createResponse(Status::CODE_404, "Not found");
        }
    }

private:
    string _job_directory = "/tmp/jobs";
};

#endif // RESULTS_REQUESTHANDLER_H