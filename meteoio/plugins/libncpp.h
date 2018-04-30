/***********************************************************************************/
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
#ifndef LIBNCPP_H
#define LIBNCPP_H

#include <meteoio/dataClasses/Grid2DObject.h>

#include <netcdf.h>
#include <string>
#include <vector>

namespace ncpp {
	void open_file(const std::string& filename, const int& omode, int& ncid);
	void create_file(const std::string& filename, const int& cmode, int& ncid);
	void end_definitions(const std::string& filename, const int& ncid);
	void close_file(const std::string& filename, const int& ncid);
	
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const int& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
	bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name);
	
	void read_data(const int& ncid, const std::string& varname, const int& varid, const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data);
	void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data);
	void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols, const size_t& pos_start, const double * const data);
} // end namespace

#endif
