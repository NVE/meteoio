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
#include "NetCDFIO.h"

using namespace std;

namespace mio {
/**
 * @page netcdf NetCDF
 * @section netcdf_format Format
 * *Put here the informations about the standard format that is implemented*
 *
 * @section netcdf_units Units
 *
 *
 * @section netcdf_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

const double NetCDFIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), coordin(""), coordinparam(""), coordout(""), coordoutparam("")
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), coordin(""), coordinparam(""), coordout(""), coordoutparam("")
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

NetCDFIO::~NetCDFIO() throw()
{

}

void NetCDFIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*name_in*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::read2DGrid(Grid2DObject& /*grid_out*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);

}

void NetCDFIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const std::string& varname)
{
	int ncid, varid;
	vector<int> dimid, dim_varid; 
	vector<string> dimname;
	vector<size_t> dimlen;

	open_file(filename, NC_NOWRITE, ncid);
	get_variable(ncid, varname, varid);
	get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

	if (dimid.size()!=2 || dimlen[0]==0 || dimlen[1]==0)
		throw IOException("Variable '" + varname + "' may only have two dimensions and both have to have length >0", AT);

	cout << "Dimensions: " << dimlen[0] << " x " << dimlen[1] << endl;

	double *lat = new double[dimlen[0]];
	double *lon = new double[dimlen[1]];
	double *grid = new double[dimlen[0]*dimlen[1]];
	
	read_data(ncid, varname, varid, grid);
	read_data(ncid, dimname[0], dim_varid[0], lat);
	read_data(ncid, dimname[1], dim_varid[1], lon);

	cout << grid[0] << endl;
	cout << lat[0] << endl;
	cout << lon[0] << endl;

	cout << "Starting to copy to Grid2DObject...";

	Coords location(coordin, coordinparam);
	location.setLatLon(lat[0], lon[0], grid[0]);
	double cellsize = calculate_cellsize(dimlen[0], dimlen[1], lat, lon);
	grid_out.set(dimlen[1], dimlen[0], cellsize, location);

	for (size_t kk=0; kk < dimlen[0]; kk++) {
		for (size_t ll=0; ll < dimlen[1]; ll++) {
			grid_out(ll, kk) = IOUtils::standardizeNodata(grid[kk*dimlen[1] + ll], plugin_nodata);
		}
	}

	close_file(filename, ncid);
	cout << "Finished" << endl;
	delete[] lat;
	delete[] lon;
	delete[] grid;
}

double NetCDFIO::calculate_cellsize(const size_t& latlen, const size_t& lonlen, double* const& lat, double* const& lon)
{
	return 30.0;
}


void NetCDFIO::open_file(const std::string& filename, const int& omode, int& ncid)
{
	int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_variable(const int& ncid, const std::string& varname, int& varid)
{
	int status = nc_inq_varid(ncid, varname.c_str(), &varid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve varid for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_dimension(const int& ncid, const std::string& varname, const int& varid, 
                             std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen)
{
	int dimids[NC_MAX_VAR_DIMS], ndimsp;

	int status = nc_inq_var(ncid, varid, NULL, NULL, &ndimsp, dimids, NULL);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimensions for variable '" + varname + "': " + nc_strerror(status), AT);

	for (int ii=0; ii<ndimsp; ii++) {
		int dimvarid;
		size_t length=0;
		char name[NC_MAX_NAME+1];

		status = nc_inq_dimname(ncid, dimids[ii], name);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);

		status = nc_inq_dimlen(ncid, dimids[ii], &length);
		if (status != NC_NOERR) throw IOException("Could not read dimension length for '" + string(name)  + "':" + nc_strerror(status), AT);

		status = nc_inq_varid(ncid, name, &dimvarid);
		if (status != NC_NOERR)
			throw IOException("Could not retrieve varid for variable '" + string(name) + "': " + nc_strerror(status), AT);

		dimid.push_back(dimids[ii]);
		dim_varid.push_back(dimvarid);
		dimname.push_back(string(name));
		dimlen.push_back(length);
	}
}

void NetCDFIO::read_data(const int& ncid, const std::string& varname, const int& varid, double*& data)
{
	int status = nc_get_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::close_file(const std::string& filename, const int& ncid)
{
	int status = nc_close(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close netcdf file  '" + filename + "': " + nc_strerror(status), AT);

}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	string filename(""), varname("");
	
	cfg.getValue("DEMFILE", "Input", filename);
	cfg.getValue("DEMVAR", "Input", varname);

	read2DGrid_internal(dem_out, filename, varname);
}

void NetCDFIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                             const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                              const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::write2DGrid(const Grid2DObject& /*grid_in*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::cleanup() throw()
{

}

} //namespace
