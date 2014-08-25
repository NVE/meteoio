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
using namespace ncpp; // wrappers for libnetcdf

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
 * <A REF="http://www.epic.noaa.gov/java/ncBrowse/">ncBrowse></A> java software.
 *
 * The <A HREF="http://cfconventions.org/1.6.html">conventions</A> for climate and forecast (CF) metadata
 * are designed to promote the processing and sharing of netCDF files. The conventions define metadata
 * that provide a definitive description of what the data represents, and the spatial and temporal properties of the data.
 * This plugin follows such conventions as well as the naming extensions defined by the
 * <A HREF="http://www.cnrm.meteo.fr/">CNRM</A>.
 *
 * *Put here the more informations about the standard format that is implemented*
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
 * - METEOFILE: the NetCDF file which shall be used for the meteo parameter input/output; [Input] and [Output] section
 * - GRID2DFILE: the NetCDF file which shall be used for gridded input/output; [Input] and [Output] section
 * - STRICTFORMAT: Whether the NetCDF file should be strictly compliant with the CNRM standard; Parameters not present
 *                 in the specification will be omitted; [Input] and [Output] section
 *
 * @section netcdf_example Example use
 * @code
 * [Input]
 * DEM     = NETCDF
 * DEMFILE = ./input/Aster_tile.nc
 * DEMVAR  = z
 * @endcode
 *
 * @section netcdf_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 */

const double NetCDFIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value
const double NetCDFIO::epsilon = 1.0e-10; //when comparing timestamps

const std::string NetCDFIO::cf_time = "time";
const std::string NetCDFIO::cf_units = "units";
const std::string NetCDFIO::cf_days = "days since ";
const std::string NetCDFIO::cf_hours = "hours since ";
const std::string NetCDFIO::cf_seconds = "seconds since ";
const std::string NetCDFIO::cf_latitude = "lat";
const std::string NetCDFIO::cf_longitude = "lon";
const std::string NetCDFIO::cf_altitude = "z";
const std::string NetCDFIO::cf_ta = "temperature";
const std::string NetCDFIO::cf_rh = "humidity";
const std::string NetCDFIO::cf_p = "pressure";

const std::string NetCDFIO::cnrm_points = "Number_of_points";
const std::string NetCDFIO::cnrm_latitude = "LAT";
const std::string NetCDFIO::cnrm_longitude = "LON";
const std::string NetCDFIO::cnrm_altitude = "ZS";
const std::string NetCDFIO::cnrm_aspect = "aspect";
const std::string NetCDFIO::cnrm_slope = "slope";
const std::string NetCDFIO::cnrm_ta = "Tair";
const std::string NetCDFIO::cnrm_rh = "HUMREL";
const std::string NetCDFIO::cnrm_vw = "Wind";
const std::string NetCDFIO::cnrm_dw = "Wind_DIR";
const std::string NetCDFIO::cnrm_qair = "Qair";
const std::string NetCDFIO::cnrm_co2air = "CO2air";
const std::string NetCDFIO::cnrm_theorsw = "theorSW";
const std::string NetCDFIO::cnrm_neb = "NEB";
const std::string NetCDFIO::cnrm_hnw = "Rainf";
const std::string NetCDFIO::cnrm_snowf = "Snowf";
const std::string NetCDFIO::cnrm_swr_direct = "DIR_SWdown";
const std::string NetCDFIO::cnrm_swr_diffuse = "SCA_SWdown";
const std::string NetCDFIO::cnrm_p = "PSurf";
const std::string NetCDFIO::cnrm_ilwr = "LWdown";
const std::string NetCDFIO::cnrm_timestep = "FRC_TIME_STP";

const std::string NetCDFIO::ecmwf_ta = "t2m";
const std::string NetCDFIO::ecmwf_p = "sp";
const std::string NetCDFIO::ecmwf_iswr = "ssrd";
const std::string NetCDFIO::ecmwf_ilwr = "strd";
const std::string NetCDFIO::ecmwf_hnw = "tp";
const std::string NetCDFIO::ecmwf_td = "d2m";
const std::string NetCDFIO::ecmwf_u10 = "u10";
const std::string NetCDFIO::ecmwf_v10 = "v10";

std::map<std::string, size_t> NetCDFIO::paramname;
std::map<std::string, std::string> NetCDFIO::map_name;
const bool NetCDFIO::__init = NetCDFIO::initStaticData();

bool NetCDFIO::initStaticData()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	paramname[cnrm_ta] = MeteoData::TA;
	paramname[cnrm_qair] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_co2air] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_neb] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_theorsw] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_rh] = MeteoData::RH;
	paramname[cnrm_vw] = MeteoData::VW;
	paramname[cnrm_dw] = MeteoData::DW;
	paramname[cnrm_hnw] = IOUtils::npos;
	paramname[cnrm_snowf] = IOUtils::npos;
	paramname[cnrm_swr_direct] = IOUtils::npos;
	paramname[cnrm_swr_diffuse] = IOUtils::npos;
	paramname[cnrm_p] = MeteoData::P;
	paramname[cnrm_ilwr] = MeteoData::ILWR;

	map_name["TA"] = cnrm_ta;
	map_name["RH"] = cnrm_rh;
	map_name["ILWR"] = cnrm_ilwr;
	map_name["P"] = cnrm_p;
	map_name["VW"] = cnrm_vw;
	map_name["DW"] = cnrm_dw;
	map_name["ISWR"] = cnrm_swr_direct;
	map_name["HNW"] = cnrm_hnw;
	map_name[cnrm_co2air] = cnrm_co2air;
	map_name[cnrm_qair] = cnrm_qair;
	map_name[cnrm_theorsw] = cnrm_theorsw;
	map_name[cnrm_neb] = cnrm_neb;

	return true;
}

NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(),
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

	cfg.getValue("STRICTFORMAT", "Input", in_strict, IOUtils::nothrow);
	cfg.getValue("STRICTFORMAT", "Output", out_strict, IOUtils::nothrow);
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
	string filename;
	cfg.getValue("GRID2DFILE", "Input", filename);

	const string varname = get_varname(parameter);

	read2DGrid_internal(grid_out, filename, varname, date);
}

void NetCDFIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const std::string& varname, const Date& date)
{
	const bool is_record = (date != Date());
	size_t lat_index = 0, lon_index = 1;

	int ncid, varid;
	vector<int> dimid, dim_varid;
	vector<string> dimname;
	vector<size_t> dimlen;

	open_file(filename, NC_NOWRITE, ncid);
	get_variable(ncid, varname, varid);
	get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

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

	double *lat = new double[dimlen[lat_index]];
	double *lon = new double[dimlen[lon_index]];
	double *grid = new double[dimlen[lat_index]*dimlen[lon_index]];

	read_data(ncid, dimname[lat_index], dim_varid[lat_index], lat);
	read_data(ncid, dimname[lon_index], dim_varid[lon_index], lon);

	if (is_record) {
		const size_t pos = find_record(ncid, NetCDFIO::cf_time, dimid[0], date.getModifiedJulianDate());
		if (pos == IOUtils::npos)
			throw IOException("No record for date " + date.toString(Date::ISO), AT);

		read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
	} else {
		read_data(ncid, varname, varid, grid);
	}

	double missing_value=plugin_nodata;
	if (ncpp::check_attribute(ncid, varid, "missing_value"))
		ncpp::get_attribute(ncid, varname, varid, "missing_value", missing_value);

	copy_grid(dimlen[lat_index], dimlen[lon_index], lat, lon, grid, missing_value, grid_out);

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

	close_file(filename, ncid);
	delete[] lat; delete[] lon; delete[] grid;
}

void NetCDFIO::copy_grid(const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
                         const double * const grid, const double& nodata, Grid2DObject& grid_out)
{
	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = calculate_cellsize(latlen, lonlen, lat, lon, resampling_factor_x, resampling_factor_y);

	const double cntr_lat = .5*(lat[0]+lat[latlen-1]);
	const double cntr_lon = .5*(lon[0]+lon[lonlen-1]);

	//computing lower left corner by using the center point as reference
	Coords cntr(coordin, coordinparam);
	cntr.setLatLon(cntr_lat, cntr_lon, IOUtils::nodata);
	cntr.moveByXY(-.5*(double)(lonlen-1)*cellsize, -.5*(double)(latlen-1)*cellsize);

	grid_out.set(lonlen, latlen, cellsize, cntr);

	//Handle the case of llcorner/urcorner swapped
	if (lat[0]<=lat[latlen-1]) {
		for (size_t kk=0; kk < latlen; kk++) {
			const size_t row = kk*lonlen;
			if (lon[0]<=lon[lonlen-1]) {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = IOUtils::standardizeNodata(grid[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = IOUtils::standardizeNodata(grid[row + (lonlen -1) - ll], nodata);
			}
		}
	} else {
		for (size_t kk=0; kk < latlen; kk++) {
			const size_t row = ((latlen-1) - kk)*lonlen;
			if (lon[0]<=lon[lonlen-1]) {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = IOUtils::standardizeNodata(grid[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < lonlen; ll++)
					grid_out(ll, kk) = IOUtils::standardizeNodata(grid[row + (lonlen -1) - ll], nodata);
			}
		}
	}

	if (resampling_factor_x != IOUtils::nodata) {
		grid_out.grid2D = ResamplingAlgorithms2D::BilinearResampling(grid_out.grid2D, resampling_factor_x, resampling_factor_y);
	}
}

/* The Grid2DObject holds data and meta data for quadratic cells. However the NetCDF file
 * stores the grid as discrete latitude and longitude values. It is necessary to calculate
 * the distance between the edges of the grid and determine the cellsize. This cellsize may
 * be different for X and Y directions. We then choose one cellsize for our grid and
 * determine a factor that will be used for resampling the grid to likewise consist of
 * quadratic cells.
 */
double NetCDFIO::calculate_cellsize(const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
                                    double& factor_x, double& factor_y)
{
	const double cntr_lat = .5*(lat[0]+lat[latlen-1]);
	const double cntr_lon = .5*(lon[0]+lon[lonlen-1]);
	double alpha;

	const double distanceX = Coords::VincentyDistance(cntr_lat, lon[0], cntr_lat, lon[lonlen-1], alpha);
	const double distanceY = Coords::VincentyDistance(lat[0], cntr_lon, lat[latlen-1], cntr_lon, alpha);

	// lonlen, latlen are decremented by 1; n linearly connected points have (n-1) connections
	const double cellsize_x = distanceX / (lonlen-1);
	const double cellsize_y = distanceY / (latlen-1);

	// round to 1cm precision for numerical stability
	const double cellsize = Optim::round( std::min(cellsize_x, cellsize_y)*100. ) / 100.;

	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
		factor_x =  cellsize_x / cellsize;
		factor_y =  cellsize_y / cellsize;

		return cellsize;
	}
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	string filename, varname;

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
	if (!vecMetaData.empty()) { // We already have meta data
		vecStation = vecMetaData;
		return;
	}

	string filename;
	cfg.getValue("METEOFILE", "Input", filename);

	int ncid;
	open_file(filename, NC_NOWRITE, ncid);
	readMetaData(ncid, vecMetaData);
	close_file(filename, ncid);

	vecStation = vecMetaData;
}

void NetCDFIO::readMetaData(const int& ncid, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	int dimid;
	size_t dimlen;
	map<string, int> map_vid;

	get_dimension(ncid, cnrm_points, dimid, dimlen);
	if (dimlen == 0) return; // There are no stations

	get_meta_data_ids(ncid, map_vid);

	double *alt = new double[dimlen];
	double *lat = new double[dimlen];
	double *lon = new double[dimlen];
	double *aspect = new double[dimlen];
	double *slope = new double[dimlen];

	read_data(ncid, cnrm_altitude, map_vid[cnrm_altitude], alt);
	read_data(ncid, cnrm_latitude, map_vid[cnrm_latitude], lat);
	read_data(ncid, cnrm_longitude, map_vid[cnrm_longitude], lon);
	read_data(ncid, cnrm_aspect, map_vid[cnrm_aspect], aspect);
	read_data(ncid, cnrm_slope, map_vid[cnrm_slope], slope);

	//Parse to StationData objects
	Coords location(coordin, coordinparam);
	ostringstream ss;
	for (size_t ii=0; ii<dimlen; ii++) {
		location.setLatLon(lat[ii], lon[ii], alt[ii]);

		ss << (ii+1);
		const string id( ss.str() );
		ss.str("");

		ss << "Station " << (ii +1);
		const string name( ss.str() );
		ss.str("");

		StationData tmp(location, id, name);
		const double aspect_bearing = (aspect[ii] < 0) ? 0 : aspect[ii]; // aspect allowed to be -1 in CNRM format...
		tmp.setSlope(slope[ii], aspect_bearing);
		vecStation.push_back(tmp);
	}

	delete[] alt; delete[] lat; delete[] lon; delete[] aspect; delete[] slope;
}

void NetCDFIO::get_meta_data_ids(const int& ncid, std::map<std::string, int>& map_vid)
{
	const string names[] = {cnrm_altitude, cnrm_latitude, cnrm_longitude, cnrm_aspect, cnrm_slope};
	vector<string> varname(names, names + sizeof(names) / sizeof(names[0]));

	vector<string> dimensions;
	dimensions.push_back(cnrm_points); // All variables have to have the dimension cnrm_points

	for (vector<string>::const_iterator it = varname.begin(); it != varname.end(); ++it) {
		int varid;
		const string& name = *it;

		get_variable(ncid, name, varid);
		if (!check_dimensions(ncid, name, varid, dimensions))
			throw IOException("Variable '" + name  + "' fails dimension check", AT);

		map_vid[name] = varid;
	}
}

void NetCDFIO::readMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, const size_t&)
{
	vecMeteo.clear();

	string filename;
	cfg.getValue("METEOFILE", "Input", filename);

	int ncid;
	open_file(filename, NC_NOWRITE, ncid);

	if (vecMetaData.empty()) readMetaData(ncid, vecMetaData);

	if (!vecMetaData.empty()) { //at least one station exists
		size_t index_start, index_end;
		vector<Date> vec_date;
		get_indices(ncid, dateStart, dateEnd, index_start, index_end, vec_date); //get indices for dateStart and dateEnd

		MeteoData meteo_data; //the template MeteoData object
		if ((index_start != IOUtils::npos) && (index_end != IOUtils::npos)) {
			map<string, size_t> map_parameters;
			get_parameters(ncid, map_parameters, meteo_data); //get a list of parameters present an render the template

			readData(ncid, index_start, vec_date, map_parameters, meteo_data, vecMeteo);
		}
	}

	close_file(filename, ncid);
}

void NetCDFIO::readData(const int& ncid, const size_t& index_start, const std::vector<Date>& vec_date,
                        const std::map<std::string, size_t>& map_parameters, const MeteoData& meteo_data, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	const size_t number_of_stations = vecMetaData.size();
	const size_t number_of_records = vec_date.size();

	// Allocate all the MeteoData objects based on the template meteo_data
	vector<MeteoData> tmp_vec(number_of_records, meteo_data);
	for (size_t jj=0; jj<number_of_records; jj++) tmp_vec[jj].date = vec_date[jj]; //set correct date for every record

	for (size_t ii=0; ii<number_of_stations; ii++) {
		for (size_t jj=0; jj<number_of_records; jj++) tmp_vec[jj].meta = vecMetaData[ii]; //adapt meta data
		vecMeteo.push_back(tmp_vec);
	}

	// Allocate enough linear space for each parameter and read the data from NetCDF
	map<string, double*> map_data;
	for (map<string, size_t>::const_iterator it = map_parameters.begin(); it != map_parameters.end(); ++it) {
		double* data = new double[number_of_stations*number_of_records];
		const string& varname = it->first;

		map_data[varname] = data;

		int varid;
		get_variable(ncid, varname, varid);
		read_data_2D(ncid, varname, varid, index_start, number_of_records, number_of_stations, data);
	}

	copy_data(ncid, map_parameters, map_data, number_of_stations, number_of_records, vecMeteo);

	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); ++it) {
		delete[] it->second;
	}
}

// The copying of data into vecMeteo is a process consisting of:
// 1. A check what the relation between MeteoIO parameters and CNRM parameters present is, check map_parameters
// 2. If there is no direct association between the parameters present and the meteo_data parameters we might
//    have to deal with the parameter in a more complex way: e.g., HNW or SWR measurements
// 3. Once we know how to deal with the parameter we loop through all stations and all parameters and copy them
//    into the appropriate places. All unit conversion have been accomplished at that point.
void NetCDFIO::copy_data(const int& ncid, const std::map<std::string, size_t>& map_parameters, const std::map<std::string, double*> map_data,
                         const size_t& number_of_stations, const size_t& number_of_records, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); ++it) {
		const string& varname = it->first;

		//find correct handling for each parameter
		bool simple_copy = false, mutiply_copy = false, hnw_measurement = false, sw_measurement = false;
		double multiplier = IOUtils::nodata;
		const size_t param = map_parameters.find(varname)->second; //must exist, at this point we know it does

		if (param == IOUtils::npos) {
			if ((varname == cnrm_snowf) || (varname == cnrm_hnw)) {
				int varid;
				get_variable(ncid, cnrm_timestep, varid);
				read_value(ncid, cnrm_timestep, varid, multiplier);

				if (multiplier <= 0) throw InvalidArgumentException("The variable '" + cnrm_timestep + "' is invalid", AT);

				hnw_measurement = true;
			} else if ((varname == cnrm_swr_diffuse) || (varname == cnrm_swr_direct)) {
				sw_measurement = true;
			} else {
				throw IOException("Don't know how to deal with parameter " + varname, AT);
			}
		} else {
			if (varname == cnrm_rh) {
				mutiply_copy = true;
				multiplier = 0.01;
			} else {
				simple_copy = true;
			}
		}

		// Loop through all times and all stations
		for (size_t jj=0; jj<number_of_records; jj++) {
			for (size_t ii=0; ii<number_of_stations; ii++) {
				double& value = (it->second)[jj*number_of_stations + ii];
				bool nodata = false;

				if (value == plugin_nodata) {
					nodata = true;
					value = IOUtils::nodata;
				}

				if (simple_copy) {
					vecMeteo[ii][jj](param) = value;
				} else if (mutiply_copy) {
					if (nodata) {
						vecMeteo[ii][jj](param) = value;
					} else {
						vecMeteo[ii][jj](param) = value * multiplier;
					}
				} else if (hnw_measurement) {
					if (!nodata) {
						double& hnw = vecMeteo[ii][jj](MeteoData::HNW);
						if (hnw == IOUtils::nodata) hnw = 0.0;
						hnw += value * multiplier;
					}
				} else if (sw_measurement) {
					if (!nodata) {
						double& iswr = vecMeteo[ii][jj](MeteoData::ISWR);
						if (iswr == IOUtils::nodata) iswr = 0.0;
						iswr += value;
					}
				}
			}
		}
	}
}

// Go through all variables present in the NetCDF dataset that have the correct dimensions. A map called
// map_parameters will associate all parameters present with MeteoData parameters or IOUtils::npos). If
// the CNRM parameter does not have a corresponding parameter in the meteo_data object we can add a new
// parameter (e.g. cnrm_theorsw) or if the situation is more complex (e.g. rainfall is measured with two
// parameters) we deal with the situation in copy_data().
void NetCDFIO::get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data)
{
	vector<string> dimensions;
	dimensions.push_back(cf_time);
	dimensions.push_back(cnrm_points);

	vector<string> present_parameters;
	get_variables(ncid, dimensions, present_parameters);

	for (vector<string>::const_iterator it = present_parameters.begin(); it != present_parameters.end(); ++it) {
		const string& name = *it;
		//cout << "Found parameter: " << name << endl;

		// Check if parameter exists in paramname, which holds strict CNRM parameters
		const map<string, size_t>::const_iterator strict_it = paramname.find(name);
		if (strict_it != paramname.end()) { // parameter is a part of the CNRM specification
			size_t index = strict_it->second;

			if ((name == cnrm_theorsw) || (name == cnrm_qair) || (name == cnrm_co2air) || (name == cnrm_neb)) {
			 	index = meteo_data.addParameter(name);
			}

			map_parameters[name] = index;
		} else if (!in_strict) { // parameter will be read anyway
			size_t index = IOUtils::npos;

			if (meteo_data.param_exists(name)) {
				index = meteo_data.getParameterIndex(name);
			} else {
			 	index = meteo_data.addParameter(name);
			}

			map_parameters[name] = index;
		}
	}
}

// The CNRM format stores timestamps as doubles (either seconds or days counted from a start date)
// This method takes the dateStart and dateEnd requested and looks for the corresponding indices
// of the time variable indexStart and indexEnd.
// Furthermore the timestamps are converted to mio::Date objects and stored in vecDate
void NetCDFIO::get_indices(const int& ncid, const Date& dateStart, const Date& dateEnd, size_t& indexStart, size_t& indexEnd, std::vector<Date>& vecDate)
{
	indexStart = indexEnd = IOUtils::npos;

	int varid, dimid;
	size_t dimlen;
	get_dimension(ncid, NetCDFIO::cf_time, dimid, dimlen);
	get_variable(ncid, NetCDFIO::cf_time, varid);

	// Get the units attribute and calculate the offset date
	string units_str;
	NetCDFIO::TimeUnit unit_type;
	Date offset;
	get_attribute(ncid, NetCDFIO::cf_time, varid, cf_units, units_str);
	calculate_offset(units_str, unit_type, offset);

	double *time = new double[dimlen];
	read_data(ncid, NetCDFIO::cf_time, varid, time);

	// Firstly, check whether search makes any sense, that is dateStart and dateEnd overlap with the times present
	bool search = true;
	if (dimlen > 0) {
		Date time_start(offset), time_end(offset);

		double start = time[0];
		double end = time[dimlen-1];

		if (unit_type == seconds) {
			start /= 86400;
			end   /= 86400;
		}
		if (unit_type == hours) {
			start /= 24;
			end   /= 24;
		}
		time_start += Date(start, 0.0);
		time_end += Date(end, 0.0);

		if (time_start > dateEnd) search = false;
		if (time_end < dateStart) search = false;
	}

	// If search is feasible then loop through the existent timestamps and find the relevant indices
	bool start_found = false;
	if (search) {
		for (size_t ii=0; ii<dimlen; ii++) {
			if (unit_type == seconds) {
				time[ii] /= 86400;
			}
			if (unit_type == hours) {
				time[ii] /= 24;
			}

			const Date tmp_date = offset + Date(time[ii], 0.0);

			if (!start_found && (dateStart <= tmp_date && tmp_date <= dateEnd)) {
				start_found = true;
				indexStart = ii;
			} else if (start_found && (tmp_date > dateEnd)) {
				indexEnd = ii-1;
				break;
			}

			if (start_found) vecDate.push_back(tmp_date);
		}

		if (start_found && (indexEnd == IOUtils::npos)) {
			indexEnd = dimlen-1;
		}
	}

	delete[] time;
}

// The CNRM timestamps have an offset that is saved in the units attribute of
// the time variable - this method retrieves that offset
void NetCDFIO::calculate_offset(const std::string& units, NetCDFIO::TimeUnit& time_unit, Date& offset)
{
	string tmp(units);
	const size_t found_sec = units.find(NetCDFIO::cf_seconds);
	const size_t found_hour = units.find(NetCDFIO::cf_hours);
	const size_t found_day = units.find(NetCDFIO::cf_days);

	if (found_sec != string::npos) {
		time_unit = seconds;
		tmp = tmp.substr(found_sec + NetCDFIO::cf_seconds.size());
	} else if (found_hour != string::npos) {
		time_unit = hours;
		tmp = tmp.substr(found_hour + NetCDFIO::cf_hours.size());
	} else if (found_day != string::npos) {
		time_unit = days;
		tmp = tmp.substr(found_day + NetCDFIO::cf_days.size());
	} else {
		throw InvalidFormatException("Variable '"+NetCDFIO::cf_time+"' has no valid attribute '" + cf_units + "'" , AT);
	}

	const bool success = IOUtils::convertString(offset, tmp, in_dflt_TZ);
	if (!success) throw InvalidFormatException("Cannot parse time: " + tmp, AT);
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	const size_t number_of_stations = vecMeteo.size();
	if (number_of_stations == 0) return; //Nothing to write

	const size_t number_of_records = vecMeteo[0].size();

	string filename;
	cfg.getValue("METEOFILE", "Output", filename);

	int ncid, did_time, vid_time, did_points;
	bool create_time = false, create_points = false, create_locations = false, create_variables = false;

	const bool exists = IOUtils::fileExists(filename);
	if (exists) remove(filename.c_str()); // NOTE: file is deleted if it exists

	double* dates;
	map<string, double*> map_data; // holds a pointer for every C array to be written
	map_data[cnrm_latitude] = new double[number_of_stations];
	map_data[cnrm_longitude] = new double[number_of_stations];
	map_data[cnrm_altitude] = new double[number_of_stations];
	map_data[cnrm_aspect] = new double[number_of_stations];
	map_data[cnrm_slope] = new double[number_of_stations];

	map<string, int> varid;
	map<size_t, string> map_param_name;

	get_parameters(vecMeteo, map_param_name, map_data, dates);

	create_file(filename, NC_CLASSIC_MODEL, ncid);
	create_time = create_points = create_locations = create_variables = true;

	if (create_time) create_time_dimension(ncid, did_time, vid_time);
	if (create_points) add_dimension(ncid, cnrm_points, number_of_stations, did_points);
	if (create_locations) create_meta_data(ncid, did_points, map_data, varid);
	if (create_variables) create_parameters(ncid, did_time, did_points, number_of_records, number_of_stations, map_param_name, map_data, varid);

	end_definitions(filename, ncid);

	copy_data(number_of_stations, number_of_records, vecMeteo, map_param_name, map_data);

	write_record(ncid, NetCDFIO::cf_time, vid_time, 0, number_of_records, dates);
	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); ++it) {
		const string& varname = it->first;
		write_data(ncid, varname, varid[varname], map_data[varname]);
		delete[] it->second;
	}

	close_file(filename, ncid);

	delete[] dates;
}

// Copy the data from the MeteoData objects into C arrays, perform all necessary
// conversions (multiplications) and set plugin_nodata values where required.
// A loop over all parameters present is performed.
void NetCDFIO::copy_data(const size_t& number_of_stations, const size_t& number_of_records, const std::vector< std::vector<MeteoData> >& vecMeteo,
                         const std::map<size_t, std::string>& map_param_name, std::map<std::string, double*>& map_data_2D)
{
	for (map<size_t, string>::const_iterator it = map_param_name.begin(); it != map_param_name.end(); ++it) {
		const size_t param = it->first;
		const string varname = it->second;

		bool simple_copy = false, multiply_copy = false;
		double multiplier = IOUtils::nodata;

		double* data = map_data_2D[varname];

		if (param == MeteoData::RH) {
			multiplier = 100.;
			multiply_copy = true;
		} else if (param == MeteoData::HNW) {
			multiply_copy = true;
			multiplier = 1./3600.;
		} else {
			simple_copy = true;
		}

		for (size_t ii=0; ii<number_of_stations; ++ii) {
			for (size_t jj=0; jj<number_of_records; ++jj) {
				const double& value = vecMeteo[ii][jj](param);

				if (value == IOUtils::nodata) {
					data[jj*number_of_stations + ii] = plugin_nodata;
				} else if (simple_copy) {
					data[jj*number_of_stations + ii] = value;
				} else if (multiply_copy) {
					data[jj*number_of_stations + ii] = value * multiplier;
				}
			}
		}
	}
}

// Create meta data variables in the NetCDF dataset
void NetCDFIO::create_meta_data(const int& ncid, const int& did, std::map<std::string, double*>& map_data_1D, std::map<std::string, int>& varid)
{
	for (map<string, double*>::const_iterator it = map_data_1D.begin(); it != map_data_1D.end(); ++it) {
		int vid;
		const string& varname = it->first;

		if (varname == cnrm_timestep) {
			add_0D_variable(ncid, cnrm_timestep, NC_DOUBLE, vid);
		} else {
			add_1D_variable(ncid, varname, NC_DOUBLE, did, vid);
		}
		add_attribute(ncid, vid, "_FillValue", plugin_nodata);
		add_attributes_for_variable(ncid, vid, varname);

		varid[varname] = vid;
	}
}

// Create the parameter variables in the NetCDF dataset, allocate memory for the
// respective C arrays and store the variable ids in the varid map.
void NetCDFIO::create_parameters(const int& ncid, const int& did_time, const int& did_points, const size_t& number_of_records,
						   const size_t& number_of_stations, std::map<size_t, std::string>& map_param_name,
                                 std::map<std::string, double*>& map_data_2D, std::map<std::string, int>& varid)
{
	// At this point map_param_name holds all parameters that have values different from nodata
	for (map<size_t, string>::iterator it = map_param_name.begin(); it != map_param_name.end();) {
		bool create = false;
		string& varname = it->second;

		const map<string, string>::const_iterator it_cnrm = map_name.find(varname);
		if (it_cnrm != map_name.end()) {
			varname = it_cnrm->second; // the offical CNRM name for the parameter
			create = true;
			++it;
		} else {
			if (out_strict) {
				// ignore any parameters not defined in the CNRM standard:
				// if a parameter in map_param_name has no equivalent in the map_name map
				// it is deleted from map_param_name and henceforth ignored.
				map_param_name.erase(it++);
			} else {
				create = true;
				++it;
			}
		}

		if (create) {
			int vid;

			double* data = new double[number_of_records*number_of_stations];
			map_data_2D[varname] = data;

			add_2D_variable(ncid, varname, NC_DOUBLE, did_time, did_points, vid);
			add_attribute(ncid, vid, "_FillValue", plugin_nodata);
			add_attributes_for_variable(ncid, vid, varname);

			varid[varname] = vid;
		}
	}
}

// Retrieve the parameters in use (parameters, that are different from nodata
// for at least one timestamp for at least one station) and store them in
// map_param_name. map_param_name associates a MeteoData parameter index with a
// string name, that is the CNRM name for the parameter to use in the NetCDF
// file. Furthermore this method copies the meta data into the appropriate C
// arrays. The timestep interval is also calculated and added to the map_data_1D
void NetCDFIO::get_parameters(const std::vector< std::vector<MeteoData> >& vecMeteo, std::map<size_t, std::string>& map_param_name,
                              std::map<std::string, double*>& map_data_1D, double*& dates)
{
	const size_t number_of_records = vecMeteo[0].size();
	dates = new double[number_of_records];

	double interval = 0;
	for (size_t ii=0; ii<number_of_records; ii++) {
		dates[ii] = vecMeteo[0][ii].date.getModifiedJulianDate();
		if (ii == 1) interval = Optim::round((dates[ii] - dates[ii-1]) * 86400.);
	}

	const size_t nr_of_parameters = (!vecMeteo[0].empty())? vecMeteo[0][0].getNrOfParameters() : 0 ;
	vector<bool> vec_param_in_use(nr_of_parameters, false);
	vector<string> vec_param_name(nr_of_parameters, "");

	//Check consistency, dates must be existent everywhere
	bool inconsistent = false;
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (number_of_records != vecMeteo[ii].size()) inconsistent = true;
		for (size_t jj=0; jj<vecMeteo[ii].size(); ++jj) {
			const MeteoData& meteo_data = vecMeteo[ii][jj];

			if (!IOUtils::checkEpsilonEquality(dates[jj], meteo_data.date.getModifiedJulianDate(), NetCDFIO::epsilon))
				inconsistent = true;

			if (jj == 0) {
				map_data_1D[cnrm_latitude][ii] = meteo_data.meta.position.getLat();
				map_data_1D[cnrm_longitude][ii] = meteo_data.meta.position.getLon();
				map_data_1D[cnrm_altitude][ii] = meteo_data.meta.position.getAltitude();
				map_data_1D[cnrm_slope][ii] = meteo_data.meta.getSlopeAngle();
				map_data_1D[cnrm_aspect][ii] = meteo_data.meta.getAzimuth();
			}

			//Check which parameters are in use
			for (size_t kk=0; kk<nr_of_parameters; ++kk) {
				if (!vec_param_in_use[kk]){
					if (meteo_data(kk) != IOUtils::nodata){
						vec_param_in_use[kk] = true;
						vec_param_name[kk] = meteo_data.getNameForParameter(kk);
					}
				}
			}
		}
	}

	if (inconsistent) throw IOException("Inconsistent dates in vecMeteo between different stations", AT);

	for (size_t kk=0; kk<nr_of_parameters; ++kk) {
		if (vec_param_in_use[kk])
			map_param_name[kk] = vec_param_name[kk];
	}


	double* timestep = new double[1];
	*timestep = interval;
	map_data_1D[cnrm_timestep] = timestep;
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
	string filename;
	cfg.getValue("GRID2DFILE", "Output", filename);

	const string varname = get_varname(parameter);

	write2DGrid_internal(grid_in, filename, varname, date);
}

void NetCDFIO::write2DGrid_internal(const Grid2DObject& grid_in, const std::string& filename, const std::string& varname, const Date& date)
{
	const bool is_record = (date != Date());
	const bool exists = IOUtils::fileExists(filename);

	double *lat_array = new double[grid_in.getNy()];
	double *lon_array = new double[grid_in.getNx()];
	double *data = new double[grid_in.getNy() * grid_in.getNx()];

	calculate_dimensions(grid_in, lat_array, lon_array);
	fill_data(grid_in, data);

	int ncid, did_lat, did_lon, did_time, vid_lat, vid_lon, vid_var, vid_time;
	bool create_dimensions(false), create_variable(false), create_time(false);

	if (exists) {
		open_file(filename, NC_WRITE, ncid);

		//check of lat/lon are defined and consistent
		if (check_dim_var(ncid, cf_latitude) && check_dim_var(ncid, cf_longitude)) {
			check_consistency(ncid, grid_in, lat_array, lon_array, did_lat, did_lon, vid_lat, vid_lon);
		} else {
			create_dimensions = true;
		}

		if (is_record) {
			//check if a time dimension/variable already exists
			if (check_dim_var(ncid, NetCDFIO::cf_time)) {
				get_dimension(ncid, NetCDFIO::cf_time, did_time);
				get_variable(ncid, NetCDFIO::cf_time, vid_time);
			} else {
				create_time = true;
			}
		}

		if (check_variable(ncid, varname)) { // variable exists
			get_variable(ncid, varname, vid_var);

			vector<int> dimid, dim_varid;
			vector<string> dimname;
			vector<size_t> dimlen;

			get_dimension(ncid, varname, vid_var, dimid, dim_varid, dimname, dimlen);

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

		start_definitions(filename, ncid);
	} else {
		create_file(filename, NC_CLASSIC_MODEL, ncid);
		add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");

		create_variable = create_dimensions = true;

		if (is_record) create_time = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);
	if (create_time) create_time_dimension(ncid, did_time, vid_time);

	if (is_record && create_variable) {
		add_3D_variable(ncid, varname, NC_DOUBLE, did_time, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, varname);
	} else if (create_variable) {
		add_2D_variable(ncid, varname, NC_DOUBLE, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, varname);
	}

	end_definitions(filename, ncid);

	if (create_dimensions) {
		write_data(ncid, cf_latitude, vid_lat, lat_array);
		write_data(ncid, cf_longitude, vid_lon, lon_array);
	}

	if (is_record) {
		size_t pos_start = add_record(ncid, NetCDFIO::cf_time, vid_time, date.getModifiedJulianDate());
		write_data(ncid, varname, vid_var, grid_in.getNy(), grid_in.getNx(), pos_start, data);
	} else {
		write_data(ncid, varname, vid_var, data);
	}

	close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void NetCDFIO::create_latlon_dimensions(const int& ncid, const Grid2DObject& grid_in, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	add_dimension(ncid, cf_latitude, grid_in.getNy(), did_lat);
	add_1D_variable(ncid, cf_latitude, NC_DOUBLE, did_lat, vid_lat);
	add_attributes_for_variable(ncid, vid_lat, cf_latitude);

	add_dimension(ncid, cf_longitude, grid_in.getNx(), did_lon);
	add_1D_variable(ncid, cf_longitude, NC_DOUBLE, did_lon, vid_lon);
	add_attributes_for_variable(ncid, vid_lon, cf_longitude);
}

void NetCDFIO::create_time_dimension(const int& ncid, int& did_time, int& vid_time)
{
	add_dimension(ncid, NetCDFIO::cf_time, NC_UNLIMITED, did_time);
	add_1D_variable(ncid, NetCDFIO::cf_time, NC_DOUBLE, did_time, vid_time); // julian day
	add_attributes_for_variable(ncid, vid_time, NetCDFIO::cf_time);
}

void NetCDFIO::fill_data(const Grid2DObject& grid, double*& data)
{
	for (size_t kk=0; kk<grid.getNy(); ++kk) {
		for (size_t ll=0; ll<grid.getNx(); ++ll) {
			data[kk*grid.getNx() + ll] = grid.grid2D(ll,kk);
		}
	}
}

// When reading or writing gridded variables we should have a consistent naming
// scheme: http://cfconventions.org/1.6.html
std::string NetCDFIO::get_varname(const MeteoGrids::Parameters& parameter)
{
	string varname = MeteoGrids::getParameterName(parameter);

	const Naming naming_scheme=full;

	if (naming_scheme==full) {
		if (parameter == MeteoGrids::TA) varname = cf_ta;
		else if (parameter == MeteoGrids::RH) varname = cf_rh;
		else if (parameter == MeteoGrids::DEM) varname = cf_altitude;
		else if (parameter == MeteoGrids::P) varname = cf_p;
	} else if (naming_scheme==cnrm) {
		if (parameter == MeteoGrids::TA) varname = cnrm_ta;
		else if (parameter == MeteoGrids::RH) varname = cnrm_rh;
		else if (parameter == MeteoGrids::DEM) varname = cnrm_altitude;
		else if (parameter == MeteoGrids::P) varname = cnrm_p;
		else if (parameter == MeteoGrids::VW) varname = cnrm_vw;
		else if (parameter == MeteoGrids::DW) varname = cnrm_dw;
		else if (parameter == MeteoGrids::ILWR) varname = cnrm_ilwr;
		else if (parameter == MeteoGrids::HNW) varname = cnrm_hnw; //HACK this should add snowf!
		else if (parameter == MeteoGrids::SLOPE) varname = cnrm_slope;
		else if (parameter == MeteoGrids::AZI) varname = cnrm_aspect;
		//HACK: iswr=dir+diff
		//HACK: U, V from vw, dw
	} else if (naming_scheme==ecmwf) {
		if (parameter == MeteoGrids::TA) varname = ecmwf_ta;
		else if (parameter == MeteoGrids::P) varname = ecmwf_p;
		else if (parameter == MeteoGrids::ISWR) varname = ecmwf_iswr;
		else if (parameter == MeteoGrids::ILWR) varname = ecmwf_ilwr;
		else if (parameter == MeteoGrids::HNW) varname = ecmwf_hnw;
		else if (parameter == MeteoGrids::U) varname = ecmwf_u10; //HACK: height correction
		else if (parameter == MeteoGrids::V) varname = ecmwf_v10; //HACK: height correction
		//HACK RH from Td, VW, DW from U, V
	}

	return varname;
}

void NetCDFIO::add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname)
{
	if (varname == cf_latitude) {
		add_attribute(ncid, varid, "standard_name", "latitude");
		add_attribute(ncid, varid, "long_name", "latitude");
		add_attribute(ncid, varid, "units", "degrees_north");
	} else if (varname == cf_longitude) {
		add_attribute(ncid, varid, "standard_name", "longitude");
		add_attribute(ncid, varid, "long_name", "longitude");
		add_attribute(ncid, varid, "units", "degrees_east");
	} else if (varname == cf_altitude) {
		add_attribute(ncid, varid, "standard_name", "altitude");
		add_attribute(ncid, varid, "long_name", "height above mean sea level");
		add_attribute(ncid, varid, "units", "m");
		add_attribute(ncid, varid, "positive", "up");
		add_attribute(ncid, varid, "axis", "Z");
	} else if (varname == cf_p) {
		add_attribute(ncid, varid, "standard_name", "air_pressure");
		add_attribute(ncid, varid, "long_name", "near surface air pressure");
		add_attribute(ncid, varid, "units", "Pa");
	} else if (varname == cf_ta) {
		add_attribute(ncid, varid, "standard_name", "air_temperature");
		add_attribute(ncid, varid, "long_name", "near surface air temperature");
		add_attribute(ncid, varid, "units", "K");
	} else if (varname == cf_rh) {
		add_attribute(ncid, varid, "standard_name", "relative humidity");
		add_attribute(ncid, varid, "long_name", "relative humidity");
		add_attribute(ncid, varid, "units", "fraction");
	} else if (varname == cf_time) {
		add_attribute(ncid, varid, "standard_name", NetCDFIO::cf_time);
		add_attribute(ncid, varid, "long_name", NetCDFIO::cf_time);
		add_attribute(ncid, varid, "units", "days since 1858-11-17 00:00:00");
	} else if (varname == NetCDFIO::cnrm_altitude) {
		add_attribute(ncid, varid, "long_name", "altitude");
		add_attribute(ncid, varid, "units", "m");
	} else if (varname == NetCDFIO::cnrm_aspect) {
		add_attribute(ncid, varid, "long_name", "slope aspect");
		add_attribute(ncid, varid, "units", "degrees from north");
	} else if (varname == NetCDFIO::cnrm_slope) {
		add_attribute(ncid, varid, "long_name", "slope angle");
		add_attribute(ncid, varid, "units", "degrees from horizontal");
	} else if (varname == NetCDFIO::cnrm_latitude) {
		add_attribute(ncid, varid, "long_name", "latitude");
		add_attribute(ncid, varid, "units", "degrees_north");
	} else if (varname == NetCDFIO::cnrm_longitude) {
		add_attribute(ncid, varid, "long_name", "longitude");
		add_attribute(ncid, varid, "units", "degrees_east");
	} else if (varname == NetCDFIO::cnrm_ta) {
		add_attribute(ncid, varid, "long_name", "Near Surface Air Temperature");
		add_attribute(ncid, varid, "units", "K");
	} else if (varname == NetCDFIO::cnrm_timestep) {
		add_attribute(ncid, varid, "long_name", "Forcing_Time_Step");
		add_attribute(ncid, varid, "units", "s");
	} else if (varname == NetCDFIO::cnrm_vw) {
		add_attribute(ncid, varid, "long_name", "Wind Speed");
		add_attribute(ncid, varid, "units", "m/s");
	} else if (varname == NetCDFIO::cnrm_dw) {
		add_attribute(ncid, varid, "long_name", "Wind Direction");
		add_attribute(ncid, varid, "units", "deg");
	} else if (varname == NetCDFIO::cnrm_swr_direct) {
		add_attribute(ncid, varid, "long_name", "Surface Incident Direct Shortwave Radiation");
		add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == NetCDFIO::cnrm_hnw) {
		add_attribute(ncid, varid, "long_name", "Rainfall Rate");
		add_attribute(ncid, varid, "units", "kg/m2/s");
	} else if (varname == NetCDFIO::cnrm_rh) {
		add_attribute(ncid, varid, "long_name", "Relative Humidity");
		add_attribute(ncid, varid, "units", "%");
	} else if (varname == NetCDFIO::cnrm_ilwr) {
		add_attribute(ncid, varid, "long_name", "Surface Incident Longwave Radiation");
		add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == NetCDFIO::cnrm_p) {
		add_attribute(ncid, varid, "long_name", "Surface Pressure");
		add_attribute(ncid, varid, "units", "Pa");
	}
}

void NetCDFIO::calculate_dimensions(const Grid2DObject& grid, double*& lat_array, double*& lon_array)
{
	lat_array[0] = grid.llcorner.getLat();
	lon_array[0] = grid.llcorner.getLon();

	// The idea is to use the difference in coordinates of the upper right and the lower left
	// corner to calculate the lat/lon intervals between cells
	Coords urcorner(grid.llcorner);
	urcorner.setGridIndex(grid.getNx()-1, grid.getNy()-1, IOUtils::nodata, true);
	grid.gridify(urcorner);

	const double lat_interval = (urcorner.getLat() - lat_array[0]) / (grid.getNy()-1);
	const double lon_interval = (urcorner.getLon() - lon_array[0]) / (grid.getNx()-1);

	// The method to use interval*ii is consistent with the corresponding
	// calculation of the Grid2DObject::gridify method -> numerical stability
	for (size_t ii=1; ii<grid.getNy(); ++ii) {
		lat_array[ii] = lat_array[0] + lat_interval*ii;
	}

	for (size_t ii=1; ii<grid.getNx(); ++ii) {
		lon_array[ii] = lon_array[0] + lon_interval*ii;
	}
}

void NetCDFIO::check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
                                 int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	size_t latlen, lonlen;

	get_dimension(ncid, cf_latitude, did_lat, latlen);
	get_dimension(ncid, cf_longitude, did_lon, lonlen);

	get_variable(ncid, cf_latitude, vid_lat);
	get_variable(ncid, cf_longitude, vid_lon);

	if ((latlen != grid.getNy()) || (lonlen != grid.getNx()))
		throw IOException("Error while writing grid - grid size and lat/lon coordinates are inconsistent", AT);

	double *lat = new double[grid.getNy()];
	double *lon = new double[grid.getNx()];

	read_data(ncid, cf_latitude, vid_lat, lat);
	read_data(ncid, cf_longitude, vid_lon, lon);

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
