// SPDX-License-Identifier: LGPL-3.0-or-later
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FSTREAM_H
#define FSTREAM_H

#include <sstream>
#include <fstream>


namespace mio {
std::string cutPathToCWD(const std::string &path);
std::string limitAccess(const char* filename);

void createTree(const char* filename, bool verbose = false);

class ofstream : public std::ofstream
{
	private:
		std::string initialize(const char* filename);
    public:
        ofstream(){};
        ofstream(const char* filename, std::ios_base::openmode mode = std::ios_base::out);
        void open(const char* filename, std::ios_base::openmode mode = std::ios_base::out);
};
}

#endif