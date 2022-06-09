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

#ifndef WPS_REQUESTHANDLER_H
#define WPS_REQUESTHANDLER_H

#include <iostream>
#include <fstream>
#include "oatpp/web/server/HttpRequestHandler.hpp"
#include "rapidxml_ns-1.13.2/rapidxml_ns.hpp"
#include "exceptions/BadRequestException.h"
#include "wps/ExecutionOperationHandler.h"

using namespace std;

// Custom request handler
class WpsRequestHandler : public oatpp::web::server::HttpRequestHandler
{
public:
    WpsRequestHandler(string job_directory) : _job_directory(job_directory)
    {
    }

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

        try
        {
            if (operationName == OPERATION_EXECUTION && operationNameNs == NS_WPS)
            {
                ExecutionOperationHandler executionOperationHandler(_job_directory);
                return executionOperationHandler.handleOperation(root_node);
            }
        }
        catch (BadRequestException &e)
        {
            cout << e.what() << endl;
            return ResponseFactory::createResponse(Status::CODE_400, e.what());
        }
        catch (exception &e)
        {
            cerr << e.what() << endl;
            return ResponseFactory::createResponse(Status::CODE_500, e.what());
        }

        return ResponseFactory::createResponse(Status::CODE_400, "Operation " + operationName + " is not supported");
    }

private:
    const char *NS_WPS = "http://www.opengis.net/wps/2.0";
    const char *OPERATION_EXECUTION = "Execute";

    string _job_directory;
};

#endif // WPS_REQUESTHANDLER_H