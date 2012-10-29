/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/IOPlugin.h>

namespace mio {

const std::string IOPlugin::header="<IOPlugin>              libname, classname,    &class,      &lib</IOPlugin>";

IOPlugin& IOPlugin::operator=(const IOPlugin& source)
{
	//since this represents a plugin on a given machine/node, since the pointers point to entry points
	//in the compiled code, they should remain valid and therefore can be copied
	if(this != &source) {
		libname = source.libname;
		classname = source.classname;
		io = source.io;
		dynLibrary = source.dynLibrary;
	}
	return *this;
}

std::ostream& operator<<(std::ostream& os, const IOPlugin& data) {
	const unsigned int pt_w=8;
	os << "<IOPlugin>" << std::setw(21) << data.libname << "," << std::setw(10) << data.classname;
	os << "," << std::showbase << std::setw(pt_w+2) << data.io;
	os << "," << std::showbase << std::setw(pt_w+2) << data.dynLibrary << "</IOPlugin>\n";
	return os;
}

} //end namespace
