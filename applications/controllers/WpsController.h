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

#ifndef WPS_CONTROLLER_H
#define WPS_CONTROLLER_H

#include <iostream>
#include <fstream>
#include "oatpp/web/server/api/ApiController.hpp"
#include "oatpp/core/macro/codegen.hpp"
#include "oatpp/core/macro/component.hpp"
#include "rapidxml_ns-1.13.2/rapidxml_ns.hpp"
#include "exceptions/BadRequestException.h"
#include "wps/ExecutionOperationHandler.h"

#include OATPP_CODEGEN_BEGIN(ApiController) //<-- Begin Codegen

using namespace std;

// Custom request handler
class WpsController : public oatpp::web::server::api::ApiController
{
public:
    WpsController(string job_directory, OATPP_COMPONENT(std::shared_ptr<ObjectMapper>, objectMapper))
        : oatpp::web::server::api::ApiController(objectMapper), _job_directory(job_directory)
    {
    }

public:
    ENDPOINT("POST", "/wps", wps, BODY_STRING(String, body))
    {
        try
        {
            rapidxml_ns::xml_document<> doc;

            vector<char> body_copy(body->begin(), body->end());
            body_copy.push_back('\0');

            doc.parse<0>(&body_copy[0]);

            rapidxml_ns::xml_node<> *root_node = doc.first_node();
            string operationName = root_node->local_name();
            string operationNameNs = root_node->namespace_uri();

            if (operationName == OPERATION_EXECUTION && operationNameNs == NS_WPS)
            {
                ExecutionOperationHandler executionOperationHandler(_job_directory);
                string responseBody = executionOperationHandler.handleOperation(root_node);
                auto response = createResponse(Status::CODE_200, responseBody);
                response->putHeader("Content-Type", "text/xml");
                return response;
            }
            OATPP_LOGW("WpsRequestHandler", "Operation '%s' is not supported", operationName.c_str());
            return createResponse(Status::CODE_400, "Operation '" + operationName + "' is not supported");
        }
        catch (BadRequestException &e)
        {
            OATPP_LOGW("WpsRequestHandler", e.what());
            return createResponse(Status::CODE_400, e.what());
        }
        catch (exception &e)
        {
            OATPP_LOGE("WpsRequestHandler", e.what());
            return createResponse(Status::CODE_500, e.what());
        }
    }

private:
    const char *NS_WPS = "http://www.opengis.net/wps/2.0";
    const char *OPERATION_EXECUTION = "Execute";

    string _job_directory;
};

#include OATPP_CODEGEN_END(ApiController) //<-- End Codegen

#endif // WPS_CONTROLLER_H