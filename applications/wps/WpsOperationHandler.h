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

#ifndef WPS_OPERATION_HANDLER_H
#define WPS_OPERATION_HANDLER_H

#include "oatpp/web/protocol/http/outgoing/ResponseFactory.hpp"
#include "oatpp/web/protocol/http/outgoing/Response.hpp"
#include "rapidxml_ns-1.13.2/rapidxml_ns.hpp"

using namespace std;

class WpsOperationHandler
{
public:
    virtual ~WpsOperationHandler(){}

    // Process incoming requests and return responses
    virtual shared_ptr<oatpp::web::protocol::http::outgoing::Response> handleOperation(rapidxml_ns::xml_node<> *root_node) = 0;
};

#endif // WPS_OPERATION_HANDLER_H