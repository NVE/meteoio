/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __IOPLUGIN_H__
#define __IOPLUGIN_H__

#include "DynamicLibrary.h"
#include "IOInterface.h"

class IOPlugin {
	public:
		std::string libname;
		std::string classname;
		IOInterface *io;
		DynamicLibrary *dynLibrary;
		
		IOPlugin(std::string _s1, std::string _s2, IOInterface *p1, DynamicLibrary *p2) : libname(_s1), classname(_s2), io(p1), dynLibrary(p2){}
		IOPlugin() : libname(""), classname(""), io(NULL), dynLibrary(NULL){}
};

#endif
