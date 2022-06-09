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

#ifndef RESULTS_REQUESTHANDLER_H
#define RESULTS_REQUESTHANDLER_H

#include <exception>
#include "oatpp/web/server/HttpRequestHandler.hpp"
#include "oatpp/core/data/stream/FileStream.hpp"
#include "oatpp/web/protocol/http/outgoing/StreamingBody.hpp"

using namespace std;

// Custom request handler
class ResultsRequestHandler : public oatpp::web::server::HttpRequestHandler
{
public:
    ResultsRequestHandler(string job_directory) : _job_directory(job_directory)
    {
    }

    // Process incoming requests and return responses
    shared_ptr<OutgoingResponse> handle(const shared_ptr<IncomingRequest> &request) override
    {
        string jobId = request->getPathVariable("jobId");
        string filepath = _job_directory + "/" + jobId + "/result.zip";
        try
        {
            auto body = std::make_shared<oatpp::web::protocol::http::outgoing::StreamingBody>(
                std::make_shared<oatpp::data::stream::FileInputStream>(filepath.c_str()));
            shared_ptr<OutgoingResponse> response = OutgoingResponse::createShared(Status::CODE_200, body);
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
    string _job_directory;
};

#endif // RESULTS_REQUESTHANDLER_H