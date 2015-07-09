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
#include "NetCDFIO.h"
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/Timer.h>
#include <meteoio/MathOptim.h>
#include <meteoio/plugins/libncpp.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page netcdf NetCDF
 * @section netcdf_format Format
 * In order to promote creation, access and sharing of scientific data, the NetCDF format has been
 * created as a machine-independent format. NetCDF (network Common Data Form) is therefore an interface
 * for array-oriented data access and a library that provides an implementation of the interface. The
 * <A HREF="http://www.unidata.ucar.edu/downloads/netcdf/index.jsp">NetCDF software</A> was developed
 * at the <A HREF="http://www.unidata.ucar.edu/">Unidata Program Center</A> in Boulder, Colorado.
 * In order to graphicaly explore the content and structure of NetCDF files, you can use the
 * <A HREF="http://www.epic.noaa.gov/java/ncBrowse/">ncBrowse</A> java software.
 *
 * The NetCDF format does not impose a specific set of metadata and therefore in order to easily exchange data
 * within a given field, it is a good idea to standardize the metadata. Several such metadata schema can be used
 * by this plugin:
 * - CF1 - the <A HREF="http://cfconventions.org">conventions</A> for climate and forecast (CF) metadata;
 * - ECMWF - from the <A HREF="http://www.ecmwf.int/">European Centre for Medium-Range Weather Forecasts</A>;
 * - CNRM - from the <A HREF="http://www.cnrm.meteo.fr/">National Centre for Meteorological Research</A>.
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
 * - NETCDF_SCHEMA: the schema to use (either CF1 or CNRM or ECMWF); [Input] and [Output] section
 *
 * @section netcdf_example Example use
 * @code
 * [Input]
 * DEM     = NETCDF
 * DEMFILE = ./input/Aster_tile.nc
 * @endcode
 *
 * @section netcdf_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 */

const double NetCDFIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value
const double NetCDFIO::epsilon = 1.0e-10; //when comparing timestamps

const std::string NetCDFIO::cf_time = "time";
const std::string NetCDFIO::cf_latitude = "lat";
const std::string NetCDFIO::cf_longitude = "lon";
const std::string NetCDFIO::cf_altitude = "z";

NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), in_attributes(), out_attributes(), coordin(), coordinparam(), coordout(), coordoutparam(),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), in_attributes(), out_attributes(), coordin(), coordinparam(), coordout(), coordoutparam(),
                                              in_dflt_TZ(0.), out_dflt_TZ(0.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::~NetCDFIO() throw() {}

void NetCDFIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
	
	initAttributesMap(cfg.get("NETCDF_SCHEMA", "Input", IOUtils::nothrow), in_attributes);
	initAttributesMap(cfg.get("NETCDF_SCHEMA", "Output", IOUtils::nothrow), out_attributes);
	/*std::cout << "In in_attributes map:\n";
	for (std::map<MeteoGrids::Parameters, NetCDFIO::attributes>::const_iterator it = in_attributes.begin(); it != in_attributes.end(); ++it){
		std::cout << MeteoGrids::getParameterName((*it).first) << " = { " << (*it).second.var << " , " << (*it).second.standard_name << " , " <<  (*it).second.long_name << " , " << (*it).second.units << "}\n";
	}
	exit(0);*/
}

void NetCDFIO::initAttributesMap(std::string schema, std::map<MeteoGrids::Parameters, NetCDFIO::attributes> &attr)
{
	if (schema.empty()) return;
	IOUtils::toUpper(schema);
	
	if (schema=="CF1") {
		attr[MeteoGrids::DEM] = attributes("z", "altitude", "height above mean sea level", "m", IOUtils::nodata);
		attr[MeteoGrids::TA] = attributes("temperature", "air_temperature", "near surface air temperature", "K", IOUtils::nodata);
		attr[MeteoGrids::RH] = attributes("humidity", "relative humidity", "relative humidity", "fraction", IOUtils::nodata);
		attr[MeteoGrids::P] = attributes("pressure", "air_pressure", "near surface air pressure", "Pa", IOUtils::nodata);
	} else if (schema=="CNRM") {
		attr[MeteoGrids::DEM] = attributes("ZS", "", "altitude", "m", IOUtils::nodata);
		attr[MeteoGrids::SLOPE] = attributes("slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata);
		attr[MeteoGrids::AZI] = attributes("aspect", "", "slope aspect", "degrees from north", IOUtils::nodata);
		attr[MeteoGrids::TA] = attributes("Tair", "", "Near Surface Air Temperature", "K", IOUtils::nodata);
		attr[MeteoGrids::RH] = attributes("HUMREL", "", "Relative Humidity", "%", IOUtils::nodata);
		attr[MeteoGrids::VW] = attributes("Wind", "", "Wind Speed", "m/s", IOUtils::nodata);
		attr[MeteoGrids::DW] = attributes("Wind_DIR", "", "Wind Direction", "deg", IOUtils::nodata);
		attr[MeteoGrids::QI] = attributes("Qair", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::HNW_L] = attributes("Rainf", "", "Rainfall Rate", "kg/m2/s", IOUtils::nodata);
		attr[MeteoGrids::HNW_S] = attributes("Snowf", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::ISW_DIR] = attributes("DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", IOUtils::nodata);
		attr[MeteoGrids::ISW_DIFF] = attributes("SCA_SWdown", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::P] = attributes("PSurf", "", "Surface Pressure", "Pa", IOUtils::nodata);
		attr[MeteoGrids::ILWR] = attributes("LWdown", "", "Surface Incident Longwave Radiation", "W/m2", IOUtils::nodata);
	} else if (schema=="ECMWF") {
		attr[MeteoGrids::TA] = attributes("t2m", "", "", "", 2.);
		attr[MeteoGrids::P] = attributes("sp", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::ISWR] = attributes("ssrd", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::ILWR] = attributes("strd", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::HNW] = attributes("tp", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::TD] = attributes("d2m", "", "", "", 2.);
		attr[MeteoGrids::U] = attributes("u10m", "", "", "", 10.);
		attr[MeteoGrids::V] = attributes("v10m", "", "", "", 10.);
	} else
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema+"\"", AT);
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
	const string filename = cfg.get("GRID2DFILE", "Input");
	
	string varname="";
	read2DGrid_internal(grid_out, filename, varname, date);
	

}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	//HACK
	const string filename = cfg.get("DEMFILE", "Input");
	const string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	if (!varname.empty()) {
		if (!read2DGrid_internal(dem_out, filename, varname))
			throw InvalidArgumentException("Variable \'"+varname+"\' not found in file \'"+filename+"\'", AT);
	} else {
		const string dem_var = in_attributes[MeteoGrids::DEM].var;
		if (!dem_var.empty() && read2DGrid_internal(dem_out, filename, dem_var)) return;
		if (read2DGrid_internal(dem_out, filename, "Band1")) return; //ASTER naming
		if (read2DGrid_internal(dem_out, filename, "z")) return; //GDAL naming
		
		throw InvalidArgumentException("The variable containing the DEM could not be found. Please specify it using the DEMVAR key.", AT);
	}
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

void NetCDFIO::readStationData(const Date&, std::vector<StationData>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readMeteoData(const Date&, const Date&, std::vector< std::vector<MeteoData> >&, const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&, const std::string&)
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
	// arguments is a string of the format filname:varname
	vector<string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() != 2)
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::write2DGrid is filename:varname", AT);

	write2DGrid_internal(grid_in, vec_argument[0], vec_argument[1]);
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const string filename = cfg.get("GRID2DFILE", "Output");
	//const string varname = get_varname(parameter);
	const string varname = "";
	write2DGrid_internal(grid_in, filename, varname, date);
}

bool NetCDFIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const std::string& varname, const Date& date)
{
	int ncid, varid;
	vector<int> dimid, dim_varid;
	vector<string> dimname;
	vector<size_t> dimlen;

	ncpp::open_file(filename, NC_NOWRITE, ncid);
	if (!ncpp::check_variable(ncid, varname)) return false;
	ncpp::get_variable(ncid, varname, varid);
	ncpp::get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

	const bool is_record = (date != Date()); //HACK: other possibility: if the file only contains 1 grid and no date -> is_record=false
	size_t lat_index = 0, lon_index = 1;
	if (is_record) { // In case we're reading a record the first index is always the record index
		lat_index = 1;
		lon_index = 2;

		if (dimid.size()!=3 || dimlen[0]<1 || dimlen[lat_index]<2 || dimlen[lon_index]<2)
			throw IOException("Variable '" + varname + "' may only have three dimensions, all have to at least have length 1", AT);
	} else if (dimid.size()==3 && dimlen[0]==1) { //in case the variable is associated with a 1 element time dimension
		lat_index = 1;
		lon_index = 2;

		if (dimlen[lat_index]<2 || dimlen[lon_index]<2)
			throw IOException("All dimensions for variable '" + varname + "' have to at least have length 1", AT);
	} else if (dimid.size()!=2 || dimlen[lat_index]<2 || dimlen[lon_index]<2) {
		throw IOException("Variable '" + varname + "' may only have two dimensions and both have to have length >1", AT);
	}

	//read latitude and longitude vectors
	double *lat = new double[dimlen[lat_index]];
	double *lon = new double[dimlen[lon_index]];
	ncpp::read_data(ncid, dimname[lat_index], dim_varid[lat_index], lat);
	ncpp::read_data(ncid, dimname[lon_index], dim_varid[lon_index], lon);

	//read gridded data
	double *grid = new double[dimlen[lat_index]*dimlen[lon_index]];
	if (is_record) {
		const size_t pos = ncpp::find_record(ncid, NetCDFIO::cf_time, dimid[0], date.getModifiedJulianDate());
		if (pos == IOUtils::npos) throw IOException("No record for date " + date.toString(Date::ISO), AT);
		ncpp::read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
	} else {
		ncpp::read_data(ncid, varname, varid, grid);
	}

	//read nodata value
	double missing_value=plugin_nodata;
	if (ncpp::check_attribute(ncid, varid, "missing_value")) ncpp::get_attribute(ncid, varname, varid, "missing_value", missing_value);

	//fill our Grid2DObject with all the data that has been read
	ncpp::copy_grid(coordin, coordinparam, dimlen[lat_index], dimlen[lon_index], lat, lon, grid, missing_value, grid_out);
	delete[] lat; delete[] lon; delete[] grid;

	//handle data packing if necessary
	if (ncpp::check_attribute(ncid, varid, "scale_factor")) {
		double scale_factor=1.;
		ncpp::get_attribute(ncid, varname, varid, "scale_factor", scale_factor);
		grid_out.grid2D *= scale_factor;
	}
	if (ncpp::check_attribute(ncid, varid, "add_offset")) {
		double add_offset=0.;
		ncpp::get_attribute(ncid, varname, varid, "add_offset", add_offset);
		grid_out.grid2D += add_offset;
	}

	ncpp::close_file(filename, ncid);
	return true;
}

void NetCDFIO::write2DGrid_internal(const Grid2DObject& grid_in, const std::string& filename, const std::string& varname, const Date& date)
{
	const bool is_record = (date != Date());
	const bool exists = IOUtils::fileExists(filename);

	double *lat_array = new double[grid_in.getNy()];
	double *lon_array = new double[grid_in.getNx()];
	double *data = new double[grid_in.getNy() * grid_in.getNx()];

	ncpp::calculate_dimensions(grid_in, lat_array, lon_array);
	ncpp::fill_grid_data(grid_in, data);

	int ncid, did_lat, did_lon, did_time, vid_lat, vid_lon, vid_var, vid_time;
	bool create_dimensions(false), create_variable(false), create_time(false);

	if (exists) {
		ncpp::open_file(filename, NC_WRITE, ncid);

		//check of lat/lon are defined and consistent
		if (ncpp::check_dim_var(ncid, cf_latitude) && ncpp::check_dim_var(ncid, cf_longitude)) {
			check_consistency(ncid, grid_in, lat_array, lon_array, did_lat, did_lon, vid_lat, vid_lon);
		} else {
			create_dimensions = true;
		}

		if (is_record) {
			//check if a time dimension/variable already exists
			if (ncpp::check_dim_var(ncid, NetCDFIO::cf_time)) {
				ncpp::get_dimension(ncid, NetCDFIO::cf_time, did_time);
				ncpp::get_variable(ncid, NetCDFIO::cf_time, vid_time);
			} else {
				create_time = true;
			}
		}

		if (ncpp::check_variable(ncid, varname)) { // variable exists
			ncpp::get_variable(ncid, varname, vid_var);

			vector<int> dimid, dim_varid;
			vector<string> dimname;
			vector<size_t> dimlen;

			ncpp::get_dimension(ncid, varname, vid_var, dimid, dim_varid, dimname, dimlen);

			if (is_record) {
				if ((dimname.size() != 3) || (dimname[0] != cf_time) || (dimname[1] != cf_latitude) || (dimname[2] != cf_longitude) || (dimlen[1]!=grid_in.getNy()) || (dimlen[2]!=grid_in.getNx()))
					throw IOException("Variable '" + varname  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			} else {
				if ((dimname[0] != cf_latitude) || (dimname[1] != cf_longitude) || (dimlen[0]!=grid_in.getNy()) || (dimlen[1]!=grid_in.getNx()))
					throw IOException("Variable '" + varname  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			}
		} else {
			create_variable = true;
		}

		ncpp::start_definitions(filename, ncid);
	} else {
		ncpp::create_file(filename, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");

		create_variable = create_dimensions = true;

		if (is_record) create_time = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);
	if (create_time) create_time_dimension(ncid, did_time, vid_time);

	if (is_record && create_variable) {
		ncpp::add_3D_variable(ncid, varname, NC_DOUBLE, did_time, did_lat, did_lon, vid_var);
		//add_attributes_for_variable(ncid, vid_var, varname);
	} else if (create_variable) {
		ncpp::add_2D_variable(ncid, varname, NC_DOUBLE, did_lat, did_lon, vid_var);
		//add_attributes_for_variable(ncid, vid_var, varname);
	}

	ncpp::end_definitions(filename, ncid);

	if (create_dimensions) {
		ncpp::write_data(ncid, cf_latitude, vid_lat, lat_array);
		ncpp::write_data(ncid, cf_longitude, vid_lon, lon_array);
	}

	if (is_record) {
		size_t pos_start = ncpp::add_record(ncid, NetCDFIO::cf_time, vid_time, date.getModifiedJulianDate());
		ncpp::write_data(ncid, varname, vid_var, grid_in.getNy(), grid_in.getNx(), pos_start, data);
	} else {
		ncpp::write_data(ncid, varname, vid_var, data);
	}

	ncpp::close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void NetCDFIO::create_latlon_dimensions(const int& ncid, const Grid2DObject& grid_in, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	ncpp::add_dimension(ncid, cf_latitude, grid_in.getNy(), did_lat);
	ncpp::add_1D_variable(ncid, cf_latitude, NC_DOUBLE, did_lat, vid_lat);
	//add_attributes_for_variable(ncid, vid_lat, cf_latitude);

	ncpp::add_dimension(ncid, cf_longitude, grid_in.getNx(), did_lon);
	ncpp::add_1D_variable(ncid, cf_longitude, NC_DOUBLE, did_lon, vid_lon);
	//add_attributes_for_variable(ncid, vid_lon, cf_longitude);
}

void NetCDFIO::create_time_dimension(const int& ncid, int& did_time, int& vid_time)
{
	ncpp::add_dimension(ncid, NetCDFIO::cf_time, NC_UNLIMITED, did_time);
	ncpp::add_1D_variable(ncid, NetCDFIO::cf_time, NC_DOUBLE, did_time, vid_time); // julian day
	//add_attributes_for_variable(ncid, vid_time, NetCDFIO::cf_time);
}

/*void NetCDFIO::add_attributes_for_variable(const int& ncid, const int& varid, const MeteoGrids& parameter)
{
	ncpp::add_attribute(ncid, varid, "standard_name", getAttribute(parameter, schema, "standard_name"));
	ncpp::add_attribute(ncid, varid, "long_name", attribute(parameter, schema, "long_name"));
	ncpp::add_attribute(ncid, varid, "units", attribute(parameter, schema, "units"));
	if (parameter==MeteoGrids::DEM) {
		ncpp::add_attribute(ncid, varid, "positive", "up");
		ncpp::add_attribute(ncid, varid, "axis", "Z");
	}
}*/

void NetCDFIO::check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
                                 int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	size_t latlen, lonlen;

	ncpp::get_dimension(ncid, cf_latitude, did_lat, latlen);
	ncpp::get_dimension(ncid, cf_longitude, did_lon, lonlen);

	ncpp::get_variable(ncid, cf_latitude, vid_lat);
	ncpp::get_variable(ncid, cf_longitude, vid_lon);

	if ((latlen != grid.getNy()) || (lonlen != grid.getNx()))
		throw IOException("Error while writing grid - grid size and lat/lon coordinates are inconsistent", AT);

	double *lat = new double[grid.getNy()];
	double *lon = new double[grid.getNx()];

	ncpp::read_data(ncid, cf_latitude, vid_lat, lat);
	ncpp::read_data(ncid, cf_longitude, vid_lon, lon);

	for (size_t ii=0; ii<latlen; ++ii) {
		if (lat_array[ii] != lat[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	for (size_t ii=0; ii<lonlen; ++ii) {
		if (lon_array[ii] != lon[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	delete[] lat; delete[] lon;
}

} //namespace
