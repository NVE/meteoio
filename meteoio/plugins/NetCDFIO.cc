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
#include <meteoio/MathOptim.h>

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
 *
 * @section example Example use
 * @code
 * [Input]
 * DEM     = NETCDF
 * DEMFILE = ./input/Aster_tile.nc
 * DEMVAR  = z
 * @endcode
 *
 * @section Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 */

const double NetCDFIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value

const std::string NetCDFIO::cf_time = "time";
const std::string NetCDFIO::cf_units = "units";
const std::string NetCDFIO::cf_days = "days since ";
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

std::map<std::string, size_t> NetCDFIO::paramname;
std::map<std::string, std::string> NetCDFIO::map_name;
const bool NetCDFIO::__init = NetCDFIO::initStaticData();

bool NetCDFIO::initStaticData()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	paramname[cnrm_ta] = MeteoData::TA;
	//paramname[cnrm_qair] = IOUtils::npos; // not a standard MeteoIO parameter
	//paramname[cnrm_co2air] = IOUtils::npos; // not a standard MeteoIO parameter
	//paramname[cnrm_neb] = IOUtils::npos; // not a standard MeteoIO parameter
	//paramname[cnrm_theorsw] = IOUtils::npos; // not a standard MeteoIO parameter
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

	return true;
}

NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), coordin(""), coordinparam(""), coordout(""), coordoutparam(""),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), coordin(""), coordinparam(""), coordout(""), coordoutparam(""),
                                              in_dflt_TZ(0.), out_dflt_TZ(0.), vecMetaData()
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

	copy_grid(dimlen[lat_index], dimlen[lon_index], lat, lon, grid, grid_out);

	close_file(filename, ncid);

	delete[] lat; delete[] lon; delete[] grid;
}

void NetCDFIO::copy_grid(const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
                         const double * const grid, Grid2DObject& grid_out)
{
	Coords location(coordin, coordinparam);
	location.setLatLon(lat[0], lon[0], grid[0]);

	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = calculate_cellsize(latlen, lonlen, lat, lon, resampling_factor_x, resampling_factor_y);

	grid_out.set(lonlen, latlen, cellsize, location);

	for (size_t kk=0; kk < latlen; kk++) {
		for (size_t ll=0; ll < lonlen; ll++) {
			grid_out(ll, kk) = IOUtils::standardizeNodata(grid[kk*lonlen + ll], plugin_nodata);
		}
	}

	if (resampling_factor_x != IOUtils::nodata) {
		grid_out.grid2D = ResamplingAlgorithms2D::BilinearResampling(grid_out.grid2D, resampling_factor_x, resampling_factor_y);
		grid_out.ncols = grid_out.grid2D.getNx();
		grid_out.nrows = grid_out.grid2D.getNy();
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
	Coords llcorner(coordin, coordinparam);
	llcorner.setLatLon(lat[0], lon[0], IOUtils::nodata);

	Coords urcorner(coordin, coordinparam);
	urcorner.setLatLon(lat[latlen-1], lon[lonlen-1], IOUtils::nodata);

	const double ll_easting=llcorner.getEasting(), ll_northing=llcorner.getNorthing();
	const double ur_easting=urcorner.getEasting(), ur_northing=urcorner.getNorthing();

	const double distanceX = ur_easting - ll_easting;
	const double distanceY = ur_northing - ll_northing;
	if(distanceX<0 || distanceY<0) {
		ostringstream ss;
		ss << "Can not compute cellsize: this is most probably due to an inappropriate input coordinate system (COORDSYS).";
		ss << "Please configure one that can accomodate (" << llcorner.getLat() << "," << llcorner.getLon() << ") - ";
		ss << "(" << urcorner.getLat() << "," << urcorner.getLon() << ")";
		throw InvalidArgumentException(ss.str(), AT);
	}

	// lonlen, latlen are decremented by 1; n linearly connected points have (n-1) connections
	const double cellsize_x = distanceX / (lonlen-1);
	const double cellsize_y = distanceY / (latlen-1);

	// We're using a precision for the cellsize that is equal to 1cm, more
	// precision makes ensuing calculations numerically instable
	const double value =  max(cellsize_x, cellsize_y) * 100.0;
	const double cellsize = floor(value) / 100.0;

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
	if (!vecMetaData.empty()) { // We already have meta data
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
		check_dimensions(ncid, name, varid, dimensions);

		map_vid[name] = varid;
	}
}

void NetCDFIO::readMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, const size_t&)
{
	vecMeteo.clear();

	string filename("");
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

// Go through all known CNRM parameters defined in the map paramname and check which ones are present
// in the current NetCDF dataset. A map called map_parameters will associate all parameters present
// with MeteoData parameters or IOUtils::npos). If the CNRM parameter does not have a corresponding
// parameter in the meteo_data object we can add a new parameter (e.g. cnrm_theorsw) or if the situation
// is more complex (e.g. rainfall is measured with two parameters) we deal with the situation in copy_data().
// Furthermore the dimensions of each present parameter are checked.
void NetCDFIO::get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data)
{
	vector<string> dimensions;
	dimensions.push_back(cnrm_points);
	dimensions.push_back(cf_time);

	for (map<string, size_t>::const_iterator it = paramname.begin(); it != paramname.end(); ++it) {
		if (check_variable(ncid, it->first)) {
			const string& name = it->first;
			size_t index = it->second;

			//cout << "Found parameter: " << name << endl;
			if ((name == cnrm_theorsw) || (name == cnrm_qair) || (name == cnrm_co2air) || (name == cnrm_neb)) {
			 	index = meteo_data.addParameter(name);
			}

			map_parameters[it->first] = index;

			// Now check the dimensions of the current variable
			int varid;
			get_variable(ncid, name, varid);
			check_dimensions(ncid, name, varid, dimensions);
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
	const size_t found_day = units.find(NetCDFIO::cf_days);

	if (found_sec != string::npos) {
		time_unit = seconds;
		tmp = tmp.substr(found_sec + NetCDFIO::cf_seconds.size());
	} else if (found_day != string::npos) {
		time_unit = days;
		tmp = tmp.substr(found_day+ + NetCDFIO::cf_days.size());
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

	string filename("");
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
		const size_t& param = it->first;
		const string& varname = it->second;

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

		for (size_t ii=0; ii<number_of_stations; ii++) {
			for (size_t jj=0; jj<number_of_records; jj++) {
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
// NOTE: if a parameter in map_param_name has no equivalent in the map_name map
//       it is deleted from map_param_name and henceforth ignored.
void NetCDFIO::create_parameters(const int& ncid, const int& did_time, const int& did_points, const size_t& number_of_records,
						   const size_t& number_of_stations, std::map<size_t, std::string>& map_param_name,
                                 std::map<std::string, double*>& map_data_2D, std::map<std::string, int>& varid)
{
	map<string, string>::const_iterator it_cnrm;

	for (map<size_t, string>::iterator it = map_param_name.begin(); it != map_param_name.end();) {
		string& varname = it->second;

		it_cnrm = map_name.find(varname);
		if (it_cnrm != map_name.end()) {
			const string& cnrm_name = it_cnrm->second;
			varname = cnrm_name;

			int vid;

			double* data = new double[number_of_records*number_of_stations];
			map_data_2D[cnrm_name] = data;

			add_2D_variable(ncid, cnrm_name, NC_DOUBLE, did_time, did_points, vid);
			add_attribute(ncid, vid, "_FillValue", plugin_nodata);
			add_attributes_for_variable(ncid, vid, varname);

			varid[varname] = vid;
			++it;
		} else {
			map_param_name.erase(it++);
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

	size_t nr_of_parameters = 0;
	if (!vecMeteo[0].empty()) nr_of_parameters = vecMeteo[0][0].getNrOfParameters();

	vector<bool> vec_param_in_use(nr_of_parameters, false);
	vector<string> vec_param_name(nr_of_parameters, "");

	//Check consistency, dates must be existent everywhere
	bool inconsistent = false;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (number_of_records != vecMeteo[ii].size()) inconsistent = true;
		for (size_t jj=0; jj<vecMeteo[ii].size(); jj++) {
			const MeteoData& meteo_data = vecMeteo[ii][jj];

			if (dates[jj] != meteo_data.date.getModifiedJulianDate()) inconsistent = true;

			if (jj == 0) {
				map_data_1D[cnrm_latitude][ii] = meteo_data.meta.position.getLat();
				map_data_1D[cnrm_longitude][ii] = meteo_data.meta.position.getLon();
				map_data_1D[cnrm_altitude][ii] = meteo_data.meta.position.getAltitude();
				map_data_1D[cnrm_slope][ii] = meteo_data.meta.getSlopeAngle();
				map_data_1D[cnrm_aspect][ii] = meteo_data.meta.getAzimuth();
			}

			//Check which parameters are in use
			for (size_t kk=0; kk<nr_of_parameters; kk++) {
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

	for (size_t kk=0; kk<nr_of_parameters; kk++) {
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
	string filename("");
	cfg.getValue("GRID2DFILE", "Output", filename);

	const string varname = get_varname(parameter);

	write2DGrid_internal(grid_in, filename, varname, date);
}

void NetCDFIO::write2DGrid_internal(const Grid2DObject& grid_in, const std::string& filename, const std::string& varname, const Date& date)
{
	const bool is_record = (date != Date());
	const bool exists = IOUtils::fileExists(filename);

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
				if ((dimname.size() != 3) || (dimname[0] != cf_time) || (dimname[1] != cf_latitude) || (dimname[2] != cf_longitude) || (dimlen[1]!=grid_in.nrows) || (dimlen[2]!=grid_in.ncols))
					throw IOException("Variable '" + varname  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			} else {
				if ((dimname[0] != cf_latitude) || (dimname[1] != cf_longitude) || (dimlen[0]!=grid_in.nrows) || (dimlen[1]!=grid_in.ncols))
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
		write_data(ncid, varname, vid_var, grid_in, pos_start, data);
	} else {
		write_data(ncid, varname, vid_var, data);
	}

	close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void NetCDFIO::create_latlon_dimensions(const int& ncid, const Grid2DObject& grid_in, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	add_dimension(ncid, cf_latitude, grid_in.nrows, did_lat);
	add_1D_variable(ncid, cf_latitude, NC_DOUBLE, did_lat, vid_lat);
	add_attributes_for_variable(ncid, vid_lat, cf_latitude);

	add_dimension(ncid, cf_longitude, grid_in.ncols, did_lon);
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
	for (size_t kk=0; kk<grid.nrows; kk++) {
		for (size_t ll=0; ll<grid.ncols; ll++) {
			data[kk*grid.ncols + ll] = grid.grid2D(ll,kk);
		}
	}
}

// When reading or writing gridded variables we should have a consistent naming
// scheme: http://cfconventions.org/1.6.html
std::string NetCDFIO::get_varname(const MeteoGrids::Parameters& parameter)
{
	string varname = MeteoGrids::getParameterName(parameter);

	if (parameter == MeteoGrids::TA) varname = cf_ta;
	else if (parameter == MeteoGrids::RH) varname = cf_rh;
	else if (parameter == MeteoGrids::DEM) varname = cf_altitude;
	else if (parameter == MeteoGrids::P) varname = cf_p;

	//TODO: complete mapping

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
	urcorner.setGridIndex(grid.ncols-1, grid.nrows-1, IOUtils::nodata, true);
	grid.gridify(urcorner);

	const double lat_interval = (urcorner.getLat() - lat_array[0]) / (grid.nrows-1);
	const double lon_interval = (urcorner.getLon() - lon_array[0]) / (grid.ncols-1);

	// The method to use interval*ii is consistent with the corresponding
	// calculation of the Grid2DObject::gridify method -> numerical stability
	for (size_t ii=1; ii<grid.nrows; ii++) {
		lat_array[ii] = lat_array[0] + lat_interval*ii;
	}

	for (size_t ii=1; ii<grid.ncols; ii++) {
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

	if ((latlen != grid.nrows) || (lonlen != grid.ncols))
		throw IOException("Error while writing grid - grid size and lat/lon coordinates are inconsistent", AT);

	double *lat = new double[grid.nrows];
	double *lon = new double[grid.ncols];

	read_data(ncid, cf_latitude, vid_lat, lat);
	read_data(ncid, cf_longitude, vid_lon, lon);

	for (size_t ii=0; ii<latlen; ii++) {
		if (lat_array[ii] != lat[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	for (size_t ii=0; ii<lonlen; ii++) {
		if (lon_array[ii] != lon[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	delete[] lat; delete[] lon;
}

//
// NetCDF C Library wrappers
//
void NetCDFIO::open_file(const std::string& filename, const int& omode, int& ncid)
{
	const int status = nc_open(filename.c_str(), omode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void NetCDFIO::create_file(const std::string& filename, const int& cmode, int& ncid)
{
	const int status = nc_create(filename.c_str(), cmode, &ncid);
	if (status != NC_NOERR)
		throw IOException("Could not create netcdf file '" + filename + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_variable(const int& ncid, const std::string& varname, int& varid)
{
	const int status = nc_inq_varid(ncid, varname.c_str(), &varid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve varid for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_dimension(const int& ncid, const std::string& dimname, int& dimid)
{
	const int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimid for dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_dimension(const int& ncid, const std::string& dimname, int& dimid, size_t& dimlen)
{
	get_dimension(ncid, dimname, dimid);

	const int status = nc_inq_dimlen(ncid, dimid, &dimlen);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve length for dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, std::string& attr_value)
{
	size_t attr_len;

     int status = nc_inq_attlen (ncid, varid, attr_name.c_str(), &attr_len);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve attribute '" + attr_name + "'for var '" + varname + "': " + nc_strerror(status), AT);

     char* value = new char[attr_len + 1]; // +1 for trailing null

     status = nc_get_att_text(ncid, varid, attr_name.c_str(), value);
	if (status != NC_NOERR)
		throw IOException("Could not read attribute '" + attr_name + "'for var '" + varname + "': " + nc_strerror(status), AT);

     value[attr_len] = '\0';
	attr_value = string(value);

	delete[] value;
}

bool NetCDFIO::check_variable(const int& ncid, const std::string& varname)
{
	int varid;
	const int status = nc_inq_varid(ncid, varname.c_str(), &varid);

	if (status != NC_NOERR) return false;

	return true;
}

bool NetCDFIO::check_dim_var(const int& ncid, const std::string& dimname)
{
	int dimid;
	const int status = nc_inq_dimid(ncid, dimname.c_str(), &dimid);
	if (status != NC_NOERR) return false;

	return check_variable(ncid, dimname);
}

void NetCDFIO::check_dimensions(const int& ncid, const std::string& varname, const int& varid, const std::vector<std::string>& names)
{
	int dimids[NC_MAX_VAR_DIMS], ndimsp;

	const int status = nc_inq_var(ncid, varid, NULL, NULL, &ndimsp, dimids, NULL);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve dimensions for variable '" + varname + "': " + nc_strerror(status), AT);

	if ((int)names.size() != ndimsp)
		throw IOException("Variable '" + varname  + "' fails dimension check", AT);

	for (int ii=0; ii<ndimsp; ii++) {
		char name[NC_MAX_NAME+1];

		const int stat = nc_inq_dimname(ncid, dimids[ii], name);
		if (stat != NC_NOERR) throw IOException(nc_strerror(stat), AT);

		const string dimname = string(name);
		const bool exists = (find(names.begin(), names.end(), dimname) != names.end());

		if (!exists)
			throw IOException("Variable '" + varname  + "' fails dimension check", AT);
	}
}

void NetCDFIO::get_dimension(const int& ncid, const std::string& varname, const int& varid,
                             std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen)
{
	dimid.clear(); dim_varid.clear(); dimname.clear(); dimlen.clear();

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

void NetCDFIO::read_data_2D(const int& ncid, const std::string& varname, const int& varid,
                            const size_t& record, const size_t& nr_of_records, const size_t& length, double*& data)
{
	size_t start[] = {record, 0};
	size_t count[] = {nr_of_records, length};

	const int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::read_value(const int& ncid, const std::string& varname, const int& varid, double& data)
{
	read_value(ncid, varname, varid, 0, data);
}

void NetCDFIO::read_value(const int& ncid, const std::string& varname, const int& varid, const size_t& pos, double& data)
{
	size_t index[] = {pos};

	const int status = nc_get_var1_double(ncid, varid, index, &data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);

}

void NetCDFIO::read_data(const int& ncid, const std::string& varname, const int& varid,
                         const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data)
{
	size_t start[] = {pos, 0, 0};
	size_t count[] = {1, latlen, lonlen};

	const int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::read_data(const int& ncid, const std::string& varname, const int& varid, double*& data)
{
	const int status = nc_get_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::write_data(const int& ncid, const std::string& varname, const int& varid, const double * const data)
{
	const int status = nc_put_var_double(ncid, varid, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::write_data(const int& ncid, const std::string& varname, const int& varid, const Grid2DObject& grid,
                          const size_t& pos_start, const double * const data)
{
	size_t start[] = {pos_start, 0, 0};
	size_t count[] = {1, grid.nrows, grid.ncols};

	const int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR) {
		throw IOException("Could not write variable '" + varname + "': " + string(nc_strerror(status)), AT);
	}
}

// Adding a record value (e.g. timestamp), in case it doesn't already exist and
// that the value is greater than the last record variable value. For example,
// timestamps have to be strictly monotonically increasing or already existent.
size_t NetCDFIO::add_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);

	//check if record already exists
	if (dimlen > 0) {
		double last_value = IOUtils::nodata;
		read_value(ncid, varname, varid, dimlen-1, last_value);

		if (last_value == data) return (dimlen - 1); //The timestamp already exists

		if (last_value > data) {
			size_t pos = find_record(ncid, varname, dimid, data); // Search for a possible match

			if (pos != IOUtils::npos) {
				return pos;
			} else {
				throw IOException("The variable '" + varname + "' has to be linearly increasing", AT);
			}
		}
	}

	write_record(ncid, varname, varid, dimlen, 1, &data);
	return dimlen;
}

// Finding a certain record variable value (e.g. timestamp) by retrieving all
// record values and then performing a linear search
size_t NetCDFIO::find_record(const int& ncid, const std::string& varname, const int& varid, const double& data)
{
	int dimid;
	size_t dimlen;

	get_dimension(ncid, varname, dimid, dimlen);

	//check if record already exists
	if (dimlen > 0) {
		double *record_value = new double[dimlen];
		read_data(ncid, varname, varid, record_value);

		for (size_t ii=0; ii<dimlen; ii++) {
			if (record_value[ii] == data) {
				delete[] record_value;
				return ii;
			}
		}

		delete[] record_value;
	}

	return IOUtils::npos; // data not found
}

// In case the dimension length of the record variable is less than start_pos
// values will be added (containing the _FillValue) until a length of start_pos-1
// has been reached. Finally the length amount elements from start_pos and on
// will be added.
void NetCDFIO::write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& start_pos, const size_t& length, const double * const data)
{
	size_t start[] = {start_pos};
	size_t count[] = {length};

	const int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for record variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid)
{
	const int status = nc_def_dim(ncid, dimname.c_str(), length, &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not define dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value)
{
	const int status = nc_put_att_double(ncid, varid, attr_name.c_str(), NC_DOUBLE, 1, &attr_value);
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value)
{
	const int status = nc_put_att_text(ncid, varid, attr_name.c_str(), attr_value.size(), attr_value.c_str());
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_0D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, int& varid)
{
	int dimid;
	const int status = nc_def_var(ncid, varname.c_str(), xtype, 0, &dimid, &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid)
{
	const int status = nc_def_var(ncid, varname.c_str(), xtype, 1, &dimid, &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid)
{
	vector<int> dimids;
	dimids.push_back(dimid1);
	dimids.push_back(dimid2);

	const int status = nc_def_var(ncid, varname.c_str(), xtype, 2, &dimids[0], &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_3D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record, const int& dimid1, const int& dimid2, int& varid)
{
	vector<int> dimids;
	dimids.push_back(dimid_record); // has to be the first one, the slowest changing index
	dimids.push_back(dimid1);
	dimids.push_back(dimid2);


	const int status = nc_def_var(ncid, varname.c_str(), xtype, 3, &dimids[0], &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::start_definitions(const std::string& filename, const int& ncid)
{
	const int status = nc_redef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not open define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void NetCDFIO::end_definitions(const std::string& filename, const int& ncid)
{
	const int status = nc_enddef(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close define mode for file '" + filename + "': " + nc_strerror(status), AT);

}

void NetCDFIO::close_file(const std::string& filename, const int& ncid)
{
	const int status = nc_close(ncid);
	if (status != NC_NOERR)
		throw IOException("Could not close netcdf file  '" + filename + "': " + nc_strerror(status), AT);

}

} //namespace
