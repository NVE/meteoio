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

#ifndef BADREQUESTEXCEPTION_H
#define BADREQUESTEXCEPTION_H

#include <exception>
#include <string>

using namespace std;

class BadRequestException : public exception
{
public:
    BadRequestException(const string &msg) : message(msg) {}

    virtual const char *what() const throw()
    {
        return message.c_str();
    }

private:
    const string message;
};

#endif // BADREQUESTEXCEPTION_H