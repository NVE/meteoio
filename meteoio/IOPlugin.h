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

#include <meteoio/DynamicLibrary.h>
#include <meteoio/IOInterface.h>

namespace mio {

/**
 * @class IOPlugin
 * @brief A helper class representing a MeteoIO plugin in the context of the dynamic loading as implemented in class
 * 	IOHandler. Each IOPlugin object represents an implementation of IOInterface or in more tangible terms it represents
 *	an instance of IOInterface in use. The two pointer member variables *io and *dynLibrary represent the opened
 *	dynamic library and the loaded object (from within that library). These pointers are essential for knowing
 *	how to deallocate the loaded object and library.
 *
 * @author Thomas Egger
 * @date   2009-08-11
 */
class IOPlugin {
	public:
		std::string libname; ///< A string representing the file to be loaded, e.g. "libgeotopio.so", can be empty
		std::string classname; ///< Classname of the object to be loaded from that dynamic library (e.g. "A3DIO")
		IOInterface *io; ///< The pointer to the actual dynamically loaded instance of IOInterface
		DynamicLibrary *dynLibrary; ///< The pointer to the opened dynamic library

		/**
		 * @brief The main constructor for the IOPlugin class
		 *
		 * @param i_s1 A std::string representing the file to be opened (or "" if plugin is statically linked)
		 * @param i_s2 A std::string that is the classname of the object to be loaded (e.g. "A3DIO", "GSNIO")
		 * @param p1  A pointer to the loaded object of type IOInterface (or NULL)
		 * @param p2  A pointer to the loaded dynamic library (or NULL)
		 */
		IOPlugin(std::string i_s1, std::string i_s2, IOInterface *p1, DynamicLibrary *p2) : libname(i_s1), classname(i_s2), io(p1), dynLibrary(p2){}
		IOPlugin() : libname(""), classname(""), io(NULL), dynLibrary(NULL){}
		IOPlugin(const IOPlugin& c) : libname(c.libname), classname(c.classname), io(c.io), dynLibrary(c.dynLibrary){}

		IOPlugin& operator=(const IOPlugin& source);

		friend std::ostream& operator<<(std::ostream& os, const IOPlugin& data);
		static const std::string header; //to contain a helpful header for understanding the output of <<

};

} //end namespace

#endif
