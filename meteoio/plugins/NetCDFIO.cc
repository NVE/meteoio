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
#include <meteoio/Timer.h>

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
 * - DEMFILE: The filename of the file containing the DEM; [Input] section 
 * - DEMVAR: The variable name of the DEM within the DEMFILE; [Input] section 
 * - GRID2DFILE: the NetCDF file which shall be used for gridded input/output; [Input] and [Output] section
 */

const double NetCDFIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const std::string NetCDFIO::time_str = "time";
const std::string NetCDFIO::lat_str = "lat";
const std::string NetCDFIO::lon_str = "lon";
const std::string NetCDFIO::z_str = "z";
const std::string NetCDFIO::ta_str = "temperature";
const std::string NetCDFIO::rh_str = "humidity";

	NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), coordin(""), coordinparam(""), coordout(""), coordoutparam(""), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), coordin(""), coordinparam(""), coordout(""), coordoutparam(""), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

NetCDFIO::~NetCDFIO() throw()
{

}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const std::string& arguments)
{
	vector<string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() == 2) {
		read2DGrid_internal(grid_out, vec_argument[0], vec_argument[1]);
	} else {
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::read2DGrid is filename:varname", AT);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	string filename("");
	cfg.getValue("GRID2DFILE", "Input", filename);

	string varname = get_varname(parameter);

	int ncid, varid;
	vector<int> dimid, dim_varid; 
	vector<string> dimname;
	vector<size_t> dimlen;

	open_file(filename, NC_NOWRITE, ncid);
	get_variable(ncid, varname, varid);
	get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

	if (dimid.size()!=3 || dimlen[0]<1 || dimlen[1]<2 || dimlen[2]<2)
		throw IOException("Variable '" + varname + "' may only have three dimensions, all have to have length >1", AT);

	cout << "Dimensions: " << dimlen[1] << " x " << dimlen[2] << " (" << dimlen[0]  << "timesteps)" << endl;

	size_t pos = find_record(ncid, NetCDFIO::time_str, dimid[0], date.getModifiedJulianDate());
	if (pos == IOUtils::npos) 
		throw IOException("No record for date " + date.toString(Date::ISO), AT);

	double *lat = new double[dimlen[1]];
	double *lon = new double[dimlen[2]];
	double *grid = new double[dimlen[1]*dimlen[2]];

	read_data(ncid, varname, varid, pos, dimlen[1], dimlen[2], grid);
	read_data(ncid, dimname[1], dim_varid[1], lat);
	read_data(ncid, dimname[2], dim_varid[2], lon);

	copy_grid(dimlen[1], dimlen[2], lat, lon, grid, grid_out);

	close_file(filename, ncid);

	delete[] lat;
	delete[] lon;
	delete[] grid;

	cout << "Grid2DObject: " << grid_out.ncols << " x " << grid_out.nrows << "  Cellsize: " << grid_out.cellsize << endl;
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

	if (dimid.size()!=2 || dimlen[0]<2 || dimlen[1]<2)
		throw IOException("Variable '" + varname + "' may only have two dimensions and both have to have length >1", AT);

	cout << "Dimensions: " << dimlen[0] << " x " << dimlen[1] << endl;

	double *lat = new double[dimlen[0]];
	double *lon = new double[dimlen[1]];
	double *grid = new double[dimlen[0]*dimlen[1]];
	
	read_data(ncid, varname, varid, grid);
	read_data(ncid, dimname[0], dim_varid[0], lat);
	read_data(ncid, dimname[1], dim_varid[1], lon);

	copy_grid(dimlen[0], dimlen[1], lat, lon, grid, grid_out);

	close_file(filename, ncid);

	delete[] lat;
	delete[] lon;
	delete[] grid;

	cout << "Grid2DObject: " << grid_out.ncols << " x " << grid_out.nrows << "  Cellsize: " << grid_out.cellsize << endl;
}

void NetCDFIO::copy_grid(const size_t& latlen, const size_t& lonlen, double*& lat, double*& lon, double*& grid, Grid2DObject& grid_out)
{
	
	cout << "Lat: " << latlen << endl;
	cout << "Lon: " << lonlen << endl;
	cout << "Starting to copy to Grid2DObject...";

	Coords location(coordin, coordinparam);
	location.setLatLon(lat[0], lon[0], grid[0]);

	double resampling_factor = IOUtils::nodata;
	double cellsize = calculate_cellsize(latlen, lonlen, lat, lon, resampling_factor);

	cout << "Detected a cellsize of: " << cellsize << endl;

	grid_out.set(lonlen, latlen, cellsize, location);

	for (size_t kk=0; kk < latlen; kk++) {
		for (size_t ll=0; ll < lonlen; ll++) {
			grid_out(ll, kk) = IOUtils::standardizeNodata(grid[kk*lonlen + ll], plugin_nodata);
		}
	}

	if (resampling_factor != IOUtils::nodata) {
		/*
		cout << "(0,0): " << grid_out.grid2D(0,0) << endl;
		cout << "(0,1): " << grid_out.grid2D(0,1) << endl;
		cout << "(1,0): " << grid_out.grid2D(1,0) << endl;
		*/
		cout << "Resampling required... " << endl;
		grid_out.grid2D = ResamplingAlgorithms2D::BilinearResampling(grid_out.grid2D, resampling_factor, 1.0);
		/*
		cout << "(0,0): " << grid_out.grid2D(0,0) << endl;
		cout << "(0,1): " << grid_out.grid2D(0,1) << endl;
		cout << "(1,0): " << grid_out.grid2D(1,0) << endl;
		*/
		grid_out.ncols = grid_out.grid2D.getNx();
		grid_out.nrows = grid_out.grid2D.getNy();
	}
	//cout << "Finished" << endl;
}

double NetCDFIO::calculate_cellsize(const size_t& latlen, const size_t& lonlen, 
                                    double* const& lat, double* const& lon, double& factor)
{
	cout << setprecision(9) << setw(20) << endl;
	cout << "Lat[0]: " << lat[0] << "   Lat[end]: " << lat[latlen-1] << endl;
	cout << "Lon[0]: " << lon[0] << "   Lon[end]: " << lon[lonlen-1] << endl;

	double alpha = 0.;
	double distanceX = Coords::cosineDistance(lat[0], lon[0], lat[0], lon[lonlen-1], alpha);
	//cout << "AlphaX: " << alpha << endl;
	double distanceY = Coords::cosineDistance(lat[0], lon[0], lat[latlen-1], lon[0], alpha);
	//cout << "AlphaY: " << alpha << endl;
	double checkX = Coords::cosineDistance(lat[0], lon[0], lat[0], lon[1], alpha);
	double checkY = Coords::cosineDistance(lat[0], lon[0], lat[1], lon[0], alpha);

	bool equal = IOUtils::checkEpsilonEquality(distanceX, distanceY, 1.0);
	bool check1 = IOUtils::checkEpsilonEquality(checkX, distanceX, 1.0);
	bool check2 = IOUtils::checkEpsilonEquality(checkY, distanceY, 1.0);

	cout << endl;

	cout << "Latlen: " << Coords::lat_degree_lenght(lat[latlen/2]) << endl;
	cout << "Lonlen: " << Coords::lon_degree_lenght(lat[latlen/2]) << endl;

	cout << "DistanceX: " <<  distanceX << endl;;
	cout << "DistanceY: " <<  distanceY << endl;;

	cout << "CellsizeX: " << (distanceX/lonlen) << endl;
	cout << "CellsizeY: " << (distanceY/latlen) << endl;

	cout << "CheckX: " << checkX << endl;
	cout << "CheckY: " << checkY << endl;


	//HACK: Check which cellsize is larger, use that one
	if (equal && check1 && check2) {
		return distanceY/latlen;
	} else {
		factor =  (distanceX/lonlen) / (distanceY/latlen);
		//cout << "Factor: " << factor << endl;
		return distanceY/latlen;
	}
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

void NetCDFIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	if (!vecMetaData.empty()) {
		vecStation = vecMetaData;
		return;
	}

	string filename("");
	cfg.getValue("METEOFILE", "Input", filename);

	int ncid;

	open_file(filename, NC_NOWRITE, ncid);
	readMetaData(ncid, vecMetaData);
	close_file(filename, ncid);

	vecStation = vecMetaData;
}

void NetCDFIO::readMetaData(const int& ncid, std::vector<StationData>& vecStation)
{
	int vid_alt, vid_lat, vid_lon, dimid;
	size_t dimlen;

	get_dimension(ncid, "Number_of_points", dimid, dimlen);
	get_variable(ncid, "ZS", vid_alt);
	get_variable(ncid, IOUtils::strToUpper(NetCDFIO::lat_str), vid_lat);
	get_variable(ncid, IOUtils::strToUpper(NetCDFIO::lon_str), vid_lon);

	double *alt = new double[dimlen];
	double *lat = new double[dimlen];
	double *lon = new double[dimlen];

	read_data(ncid, "ZS", vid_alt, alt);
	read_data(ncid, NetCDFIO::lat_str, vid_lat, lat);
	read_data(ncid, NetCDFIO::lon_str, vid_lon, lon);

	//Parse to StationData objects
	Coords location(coordin, coordinparam);
	ostringstream ss;
	string id, name;
	for (size_t ii=0; ii<dimlen; ii++) {
		location.setLatLon(lat[ii], lon[ii], alt[ii]);
		
		ss << (ii+1);
		id = ss.str();
		ss.str("");

		ss << "Station " << (ii +1);
		name = ss.str();
		ss.str("");
		
		StationData tmp(location, id, name);
		vecStation.push_back(tmp);
	}

	delete[] alt; delete[] lat; delete[] lon;
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

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const std::string& arguments)
{
	vector<string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() != 2)
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::write2DGrid is filename:varname", AT);

	const string& filename = vec_argument[0];
	const string& varname = vec_argument[1];

	bool exists = IOUtils::fileExists(filename);
	
	double *lat_array = new double[grid_in.nrows];
	double *lon_array = new double[grid_in.ncols];
	double *data = new double[grid_in.nrows * grid_in.ncols];

	calculate_dimensions(grid_in, lat_array, lon_array);
	fill_data(grid_in, data);

	int ncid, did_lat, did_lon, vid_lat, vid_lon, vid_var;
	bool create_dimensions(false), create_variable(false);

	if (exists) {
		open_file(filename, NC_WRITE, ncid);

		//check of lat/lon are defined and consistent
		if (check_dim_var(ncid, NetCDFIO::lat_str) && check_dim_var(ncid, NetCDFIO::lon_str)) {
			check_consistency(ncid, grid_in, lat_array, lon_array, did_lat, did_lon, vid_lat, vid_lon);
		} else {
			create_dimensions = true;
		}

		if (check_variable(ncid, varname)) { // variable exists
			get_variable(ncid, varname, vid_var);

			vector<int> dimid, dim_varid;
			vector<string> dimname;
			vector<size_t> dimlen;

			get_dimension(ncid, vec_argument[1], vid_var, dimid, dim_varid, dimname, dimlen);

			if ((dimname[0] != NetCDFIO::lat_str) || (dimname[1] != NetCDFIO::lon_str) || (dimlen[0]!=grid_in.nrows) || (dimlen[1]!=grid_in.ncols))
				throw IOException("Variable '" + vec_argument[1]  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
		} else {
			create_variable = true;
		}

		start_definitions(filename, ncid);
	} else {
		create_file(filename, NC_CLASSIC_MODEL, ncid);
		add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");

		create_variable = create_dimensions = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);

	if (create_variable) {
		add_2D_variable(ncid, vec_argument[1], NC_DOUBLE, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, vec_argument[1]);
	}

	end_definitions(filename, ncid);

	if (create_dimensions) {
		write_data(ncid, NetCDFIO::lat_str, vid_lat, lat_array);
		write_data(ncid, NetCDFIO::lon_str, vid_lon, lon_array);
	}

	write_data(ncid, vec_argument[1], vid_var, data);

	close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	string filename("");
	cfg.getValue("GRID2DFILE", "Output", filename);

	string varname = get_varname(parameter);
	bool exists = IOUtils::fileExists(filename);
	
	double *lat_array = new double[grid_in.nrows];
	double *lon_array = new double[grid_in.ncols];
	double *data = new double[grid_in.nrows * grid_in.ncols];

	calculate_dimensions(grid_in, lat_array, lon_array);
	fill_data(grid_in, data);

	int ncid, did_lat, did_lon, did_time, vid_lat, vid_lon, vid_var, vid_time;
	bool create_dimensions(false), create_variable(false), create_time(false);

	if (exists) {
		open_file(filename, NC_WRITE, ncid);

		//check of lat/lon are defined and consistent
		if (check_dim_var(ncid, NetCDFIO::lat_str) && check_dim_var(ncid, NetCDFIO::lon_str)) {
			check_consistency(ncid, grid_in, lat_array, lon_array, did_lat, did_lon, vid_lat, vid_lon);
		} else {
			create_dimensions = true;
		}

		//check if a time dimension/variable already exists
		if (check_dim_var(ncid, NetCDFIO::time_str)) {
			get_dimension(ncid, NetCDFIO::time_str, did_time);
			get_variable(ncid, NetCDFIO::time_str, vid_time);			
		} else {
			create_time = true;
		}

		if (check_variable(ncid, varname)) { // variable exists
			get_variable(ncid, varname, vid_var); //HACK check dimensionality
		} else {
			create_variable = true;
		}

		start_definitions(filename, ncid);
	} else {
		create_file(filename, NC_CLASSIC_MODEL, ncid);
		add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");

		create_variable = create_dimensions = create_time = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);
	if (create_time) create_time_dimension(ncid, did_time, vid_time);

	if (create_variable) {
		add_2D_record(ncid, varname, NC_DOUBLE, did_time, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, varname);
	}

	end_definitions(filename, ncid);

	if (create_dimensions) {
		write_data(ncid, NetCDFIO::lat_str, vid_lat, lat_array);
		write_data(ncid, NetCDFIO::lon_str, vid_lon, lon_array);
	}

	size_t pos_start = append_record(ncid, NetCDFIO::time_str, vid_time, date.getModifiedJulianDate());
	write_data(ncid, varname, vid_var, grid_in, pos_start, data);

	close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

std::string NetCDFIO::get_varname(const MeteoGrids::Parameters& parameter)
{
	string varname("varname");

	if (parameter == MeteoGrids::TA) varname = NetCDFIO::ta_str;
	else if (parameter == MeteoGrids::RH) varname = NetCDFIO::rh_str;
	else if (parameter == MeteoGrids::DEM) varname = NetCDFIO::z_str;

	return varname;
}

void NetCDFIO::create_latlon_dimensions(const int& ncid, const Grid2DObject& grid_in, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	define_dimension(ncid, NetCDFIO::lat_str, grid_in.nrows, did_lat);
	add_1D_variable(ncid, NetCDFIO::lat_str, NC_DOUBLE, did_lat, vid_lat);
	add_attributes_for_variable(ncid, vid_lat, NetCDFIO::lat_str);

	define_dimension(ncid, NetCDFIO::lon_str, grid_in.ncols, did_lon);
	add_1D_variable(ncid, NetCDFIO::lon_str, NC_DOUBLE, did_lon, vid_lon);
	add_attributes_for_variable(ncid, vid_lon, NetCDFIO::lon_str);
}

void NetCDFIO::create_time_dimension(const int& ncid, int& did_time, int& vid_time)
{
	define_dimension(ncid, NetCDFIO::time_str, NC_UNLIMITED, did_time);
	add_1D_variable(ncid, NetCDFIO::time_str, NC_DOUBLE, did_time, vid_time); // julian day
	add_attributes_for_variable(ncid, vid_time, NetCDFIO::time_str);
}

void NetCDFIO::fill_data(const Grid2DObject& grid, double*& data)
{
	for (size_t kk=0; kk<grid.nrows; kk++) {
		for (size_t ll=0; ll<grid.ncols; ll++) {
			data[kk*grid.ncols + ll] = grid.grid2D(ll,kk);
		}
	}
}

void NetCDFIO::add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname)
{
	if (varname == NetCDFIO::lat_str) {
		add_attribute(ncid, varid, "standard_name", "latitude");
		add_attribute(ncid, varid, "long_name", "latitude");
		add_attribute(ncid, varid, "units", "degreesnorth");
	} else if (varname == NetCDFIO::lon_str) {
		add_attribute(ncid, varid, "standard_name", "longitude");
		add_attribute(ncid, varid, "long_name", "longitude");
		add_attribute(ncid, varid, "units", "degreeseast");
	} else if (varname == NetCDFIO::z_str) {
		add_attribute(ncid, varid, "standard_name", "altitude");
		add_attribute(ncid, varid, "long_name", "height above mean sea level");
		add_attribute(ncid, varid, "units", "m");
		add_attribute(ncid, varid, "positive", "up");
		add_attribute(ncid, varid, "axis", "Z");
	} else if (varname == NetCDFIO::time_str) {
		add_attribute(ncid, varid, "standard_name", NetCDFIO::time_str);
		add_attribute(ncid, varid, "long_name", NetCDFIO::time_str);
		add_attribute(ncid, varid, "units", "days since 1858-11-17 00:00:00");
	}
}

void NetCDFIO::calculate_dimensions(const Grid2DObject& grid, double*& lat_array, double*& lon_array)
{
	lat_array[0] = grid.llcorner.getLat();
	lon_array[0] = grid.llcorner.getLon();

	Coords tmp_coord(grid.llcorner);

	for (size_t ii=1; ii<grid.nrows; ii++) {
		tmp_coord.setGridIndex(0, ii, IOUtils::nodata, true); // one step North
	     grid.gridify(tmp_coord);
		lat_array[ii] = tmp_coord.getLat();
	}
	
	for (size_t ii=1; ii<grid.ncols; ii++) {
		tmp_coord.setGridIndex(ii, 0, IOUtils::nodata, true); // one step East
		grid.gridify(tmp_coord);
		lon_array[ii] = tmp_coord.getLon();
	}
}

void NetCDFIO::check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
                                 int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	size_t latlen, lonlen;

	get_dimension(ncid, NetCDFIO::lat_str, did_lat, latlen);
	get_dimension(ncid, NetCDFIO::lon_str, did_lon, lonlen);

	get_variable(ncid, NetCDFIO::lat_str, vid_lat);
	get_variable(ncid, NetCDFIO::lon_str, vid_lon);

	if ((latlen != grid.nrows) || (lonlen != grid.ncols))
		throw IOException("Error while writing grid - grid size and lat/lon coordinates are inconsistent", AT);

	double *lat = new double[grid.nrows];
	double *lon = new double[grid.ncols];

	read_data(ncid, NetCDFIO::lat_str, vid_lat, lat);
	read_data(ncid, NetCDFIO::lon_str, vid_lon, lon);

	for (size_t ii=0; ii<latlen; ii++) {
		if (lat_array[ii] != lat[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	for (size_t ii=0; ii<lonlen; ii++) {
		if (lon_array[ii] != lon[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}
}

//
// NetCDF C Library wrappers
//
void NetCDFIO::open_file(const std::string& filename, const int& omode, int& ncid)
{
	int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void NetCDFIO::create_file(const std::string& filename, const int& cmode, int& ncid)
{
	int status = nc_create(filename.c_str(), cmode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not create netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_variable(const int& ncid, const std::string& varname, int& varid)
{
	int status = nc_inq_varid(ncid, varname.c_str(), &varid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve varid for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_dimension(const int& ncid, const std::string& dimname, int& dimid)
{
	int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimid for dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_dimension(const int& ncid, const std::string& dimname, int& dimid, size_t& dimlen)
{
	int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimid for dimension '" + dimname + "': " + nc_strerror(status), AT);

	status = nc_inq_dimlen(ncid, dimid, &dimlen);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve length for dimension '" + dimname + "': " + nc_strerror(status), AT);
}

bool NetCDFIO::check_variable(const int& ncid, const std::string& varname)
{
	int varid;
	int status = nc_inq_varid(ncid, varname.c_str(), &varid);

	if (status != NC_NOERR) return false;

	return true;
}

bool NetCDFIO::check_dim_var(const int& ncid, const std::string& dimname)
{
	int dimid;
	int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR) return false;

	return check_variable(ncid, dimname);
}

size_t NetCDFIO::get_1D_var_len(const int& ncid, const std::string& varname)
{
	int dimidp;
	size_t length = 0;

	int status = nc_inq_dimid(ncid, varname.c_str(), &dimidp);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimid for dimension '" + varname + "': " + nc_strerror(status), AT);

	status = nc_inq_dimlen(ncid, dimidp, &length);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dim length for dimension '" + varname + "': " + nc_strerror(status), AT);

	return length;
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

void NetCDFIO::read_data(const int& ncid, const std::string& varname, const int& varid,
                         const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data)
{
	size_t start[] = {pos, 0, 0};
     size_t count[] = {1, latlen, lonlen};

	int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::read_data(const int& ncid, const std::string& varname, const int& varid, double*& data)
{
	int status = nc_get_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::write_data(const int& ncid, const std::string& varname, const int& varid, double*& data)
{
	int status = nc_put_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::write_data(const int& ncid, const std::string& varname, const int& varid, const Grid2DObject& grid, const size_t& pos_start, double*& data)
{
	size_t start[] = {pos_start, 0, 0};
	size_t count[] = {1, grid.nrows, grid.ncols};

	int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR) {
		throw IOException("Could not write variable '" + varname + "': " + string(nc_strerror(status)), AT);
	}
}

size_t NetCDFIO::find_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);

	//check if record already exists
	if (dimlen > 0) {
		double *timesteps = new double[dimlen];
		read_data(ncid, varname, varid, timesteps);

		for (size_t ii=0; ii<dimlen; ii++) {
			if (timesteps[ii] == data) {
				delete[] timesteps;
				return ii;
			}
		}
		
		delete[] timesteps;
	}

	return IOUtils::npos; // data not found
}

size_t NetCDFIO::append_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid, status;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);

	//check if record already exists
	if (dimlen > 0) {
		double last_value = IOUtils::nodata;
		const size_t index_read[] = {dimlen-1};

		status = nc_get_var1_double(ncid, varid, index_read, &last_value);
		if (status != NC_NOERR)
			throw IOException("Could not retrieve last value for record '" + varname + "': " + nc_strerror(status), AT);
		cout << "Last time value: " << last_value << endl;

		if (last_value == data) return (dimlen - 1); //The timestamp already exists
	}

	const size_t index[] = {dimlen}; //append at the end

	status = nc_put_var1_double(ncid, varid, index, &data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for record '" + varname + "': " + nc_strerror(status), AT);

	return dimlen;
}

void NetCDFIO::define_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid)
{
	int status = nc_def_dim(ncid, dimname.c_str(), length, &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not define dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value)
{
	int status = nc_put_att_text(ncid, varid, attr_name.c_str(), attr_value.size(), attr_value.c_str());
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid)
{
	int status = nc_def_var(ncid, varname.c_str(), xtype, 1, &dimid, &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid)
{
	vector<int> dimids;
	dimids.push_back(dimid1);
	dimids.push_back(dimid2);

	int status = nc_def_var(ncid, varname.c_str(), xtype, 2, &dimids[0], &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_2D_record(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record, const int& dimid1, const int& dimid2, int& varid)
{
	vector<int> dimids;
	dimids.push_back(dimid_record); // has to be the first one, the slowest changing index
	dimids.push_back(dimid1);
	dimids.push_back(dimid2);


	int status = nc_def_var(ncid, varname.c_str(), xtype, 3, &dimids[0], &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::start_definitions(const std::string& filename, const int& ncid)
{
	int status = nc_redef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void NetCDFIO::end_definitions(const std::string& filename, const int& ncid)
{
	int status = nc_enddef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void NetCDFIO::close_file(const std::string& filename, const int& ncid)
{
	int status = nc_close(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close netcdf file  '" + filename + "': " + nc_strerror(status), AT);

}

} //namespace
