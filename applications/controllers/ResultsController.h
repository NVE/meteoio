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

#ifndef RESULTS_CONTROLLER_H
#define RESULTS_CONTROLLER_H

#include <exception>
#include "oatpp/web/server/api/ApiController.hpp"
#include "oatpp/core/macro/codegen.hpp"
#include "oatpp/core/macro/component.hpp"
#include "oatpp/core/data/stream/FileStream.hpp"
#include "oatpp/web/protocol/http/outgoing/StreamingBody.hpp"

#include OATPP_CODEGEN_BEGIN(ApiController) //<-- Begin Codegen

using namespace std;

// Custom request handler
class ResultsController : public oatpp::web::server::api::ApiController
{
public:
    ResultsController(string job_directory, OATPP_COMPONENT(std::shared_ptr<ObjectMapper>, objectMapper))
        : oatpp::web::server::api::ApiController(objectMapper), _job_directory(job_directory)
    {
    }

public:
    ENDPOINT("GET", "/results/{jobId}/result.zip", results, PATH(String, jobId))
    {
        string filepath = "";
        try
        {
            filepath = _job_directory + "/" + jobId + "/result.zip";
            auto body = std::make_shared<oatpp::web::protocol::http::outgoing::StreamingBody>(
                std::make_shared<oatpp::data::stream::FileInputStream>(filepath.c_str()));
            shared_ptr<OutgoingResponse> response = OutgoingResponse::createShared(Status::CODE_200, body);
            response->putHeader("Content-Type", "application/zip");
            return response;
        }
        catch (exception &e)
        {
            OATPP_LOGW("ResultsRequestHandler", "Requested file '%s' could not be served.", filepath.c_str());
            OATPP_LOGW("ResultsRequestHandler", e.what());
            return ResponseFactory::createResponse(Status::CODE_404, "Not found");
        }
    }

private:
    string _job_directory;
};

#include OATPP_CODEGEN_END(ApiController) //<-- End Codegen

#endif // RESULTS_CONTROLLER_H