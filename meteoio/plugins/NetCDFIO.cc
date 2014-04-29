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
 * - METEOFILE: the NetCDF file which shall be used for the meteo parameter input/output; [Input] and [Output] section
 * - GRID2DFILE: the NetCDF file which shall be used for gridded input/output; [Input] and [Output] section
 */

const double NetCDFIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value
const std::string NetCDFIO::lat_str = "lat";
const std::string NetCDFIO::lon_str = "lon";
const std::string NetCDFIO::z_str = "z";
const std::string NetCDFIO::ta_str = "temperature";
const std::string NetCDFIO::rh_str = "humidity";

const std::string NetCDFIO::cf_time = "time";
const std::string NetCDFIO::cf_units = "units";
const std::string NetCDFIO::cf_days = "days since ";
const std::string NetCDFIO::cf_seconds = "seconds since ";

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

NetCDFIO::~NetCDFIO() throw()
{

}

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

	size_t pos = find_record(ncid, NetCDFIO::cf_time, dimid[0], date.getModifiedJulianDate());
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
	int vid_alt, vid_lat, vid_lon, vid_aspect, vid_slope, dimid;
	size_t dimlen;

	//HACK: could check dimension of all the vars, must be 'Number_of_points'
	get_dimension(ncid, cnrm_points, dimid, dimlen);
	get_variable(ncid, cnrm_altitude, vid_alt);
	get_variable(ncid, IOUtils::strToUpper(NetCDFIO::lat_str), vid_lat);
	get_variable(ncid, IOUtils::strToUpper(NetCDFIO::lon_str), vid_lon);
	get_variable(ncid, cnrm_aspect, vid_aspect);
	get_variable(ncid, cnrm_slope, vid_slope);

	double *alt = new double[dimlen];
	double *lat = new double[dimlen];
	double *lon = new double[dimlen];
	double *aspect = new double[dimlen];
	double *slope = new double[dimlen];

	read_data(ncid, "ZS", vid_alt, alt);
	read_data(ncid, NetCDFIO::lat_str, vid_lat, lat);
	read_data(ncid, NetCDFIO::lon_str, vid_lon, lon);
	read_data(ncid, NetCDFIO::lat_str, vid_aspect, aspect);
	read_data(ncid, NetCDFIO::lat_str, vid_slope, slope);

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
		double aspect_bearing = (aspect[ii] < 0) ? 0 : aspect[ii]; // aspect allowed to be -1... HACK
		tmp.setSlope(slope[ii], aspect_bearing);
		vecStation.push_back(tmp);

		//cout << "Station " << (ii+1) << "  alt: " << alt[ii] << "  aspect: " << aspect_bearing << "  slope: " << slope[ii] << endl;
	}

	delete[] alt; delete[] lat; delete[] lon; delete[] aspect; delete[] slope;
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
			get_parameters(ncid, map_parameters, meteo_data); //get a list of parameters present
			
			readData(ncid, index_start, vec_date, map_parameters, meteo_data, vecMeteo);
		}
	}

	close_file(filename, ncid);
}

void NetCDFIO::readData(const int& ncid, const size_t& index_start, const std::vector<Date>& vec_date, 
                        const std::map<std::string, size_t>& map_parameters, const MeteoData& meteo_data, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	size_t number_of_stations = vecMetaData.size();
	size_t number_of_records = vec_date.size();

	//Allocate all the MeteoData objects
	vector<MeteoData> tmp_vec(number_of_records, meteo_data);
	for (size_t jj=0; jj<number_of_records; jj++) tmp_vec[jj].date = vec_date[jj]; //set correct date for every record

	for (size_t ii=0; ii<number_of_stations; ii++) {
		for (size_t jj=0; jj<number_of_records; jj++) tmp_vec[jj].meta = vecMetaData[ii]; //adapt meta data
		vecMeteo.push_back(tmp_vec);
	}

	//allocate enough linear space for each parameter
	map<string, double*> map_data;
	for (map<string, size_t>::const_iterator it = map_parameters.begin(); it != map_parameters.end(); it++) {
		double* data = new double[number_of_stations*number_of_records];
		const string& varname = it->first;

		map_data[varname] = data;

		int varid; 
		get_variable(ncid, varname, varid);
		read_data_2D(ncid, varname, varid, index_start, number_of_records, number_of_stations, data);
	}

	copy_data(ncid, map_parameters, map_data, number_of_stations, number_of_records, vecMeteo); 

	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); it++) {
		delete[] it->second;
	}
}

void NetCDFIO::copy_data(const int& ncid, const std::map<std::string, size_t>& map_parameters, const std::map<std::string, double*> map_data, 
                         const size_t& number_of_stations, const size_t& number_of_records, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); it++) {
		const string& varname = it->first;

		//find correct handling for each parameter
		bool simple_copy = false, mutiply_copy = false, hnw_measurement = false, sw_measurement = false;
		double multiplier = IOUtils::nodata;
		size_t param = map_parameters.find(varname)->second; //must exist, at this point we know it does

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

void NetCDFIO::get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data)
{
	for (map<string, size_t>::const_iterator it = paramname.begin(); it != paramname.end(); it++) {
		if (check_variable(ncid, it->first)) {
			const string& name = it->first;
			size_t index = it->second;

			//cout << "Found parameter: " << name << endl; 
			if ((name == cnrm_theorsw) || (name == cnrm_qair) || (name == cnrm_co2air) || (name == cnrm_neb)) {
			 	index = meteo_data.addParameter(name);
			}

			map_parameters[it->first] = index;
		}
	}

	//TODO: check dimensions?
}

void NetCDFIO::get_indices(const int& ncid, const Date& dateStart, const Date& dateEnd, size_t& indexStart, size_t& indexEnd, std::vector<Date>& vecDate)
{
	//get dimid for time
	int varid, dimid;
	size_t dimlen;
	get_dimension(ncid, NetCDFIO::cf_time, dimid, dimlen);
	get_variable(ncid, NetCDFIO::cf_time, varid);

	//get attributes, calculate offset date
	string units_str;
	NetCDFIO::TimeUnit unit_type;
	Date offset;
	get_attribute(ncid, NetCDFIO::cf_time, varid, NetCDFIO::cf_units, units_str);
	calculate_offset(units_str, unit_type, offset); //HACK: should only be exctracted once

	//read values, find indices
	double *time = new double[dimlen];
	read_data(ncid, NetCDFIO::cf_time, varid, time);

	bool start_found = false;
	indexStart = indexEnd = IOUtils::npos;

	//check whether search makes any sense
	bool search = true;
	if (dimlen > 0) {
		Date time_start(offset), time_end(offset);
	
		double start = time[0];
		double end = time[dimlen-1];

		if (unit_type == seconds) {
			start /= 86400;
			end   /= 86400;
		}
		time_start += Date(start, in_dflt_TZ);
		time_end += Date(end, in_dflt_TZ);

		if (time_start > dateEnd) search = false;
		if (time_end < dateStart) search = false;
	}

	if (search) {
		for (size_t ii=0; ii<dimlen; ii++) {
			if (unit_type == seconds) {
				time[ii] /= 86400;
			}

			Date tmp_date = offset + Date(time[ii], in_dflt_TZ);
			
			//cout << ii << "  julian: " << time[ii];
			//cout << "\t" << tmp_date.toString(Date::ISO) << endl;
			
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
	/*
	cout << "vecDate:" << vecDate.size() << endl;
	vector<Date>::iterator it;
	for (it = vecDate.begin(); it != vecDate.end(); it++) {
		cout << (*it).toString() << endl;
	}
	*/
	//cout << dateStart.toString(Date::ISO) << "  -  " << dateEnd.toString(Date::ISO) << endl;
	//cout << "indexStart: " << indexStart << "  indexEnd: " << indexEnd << endl;

	delete[] time;
}

void NetCDFIO::calculate_offset(const std::string& units, NetCDFIO::TimeUnit& time_unit, Date& offset)
{
	string tmp(units);
	size_t found_sec = units.find(NetCDFIO::cf_seconds);
	size_t found_day = units.find(NetCDFIO::cf_days);

	if (found_sec != string::npos) {
		time_unit = seconds;
		tmp = tmp.substr(found_sec + NetCDFIO::cf_seconds.size());
	} else if (found_day != string::npos) {
		time_unit = days;
		tmp = tmp.substr(found_day+ + NetCDFIO::cf_days.size());
	} else {
		throw InvalidFormatException("Variable '"+NetCDFIO::cf_time+"' has no valid attribute 'units'" , AT);
	}
	
	bool success = IOUtils::convertString(offset, tmp, in_dflt_TZ);
	if (!success) throw InvalidFormatException("Cannot parse time: " + tmp, AT);

	//cout << "Parsing : " << tmp << endl;
	//cout << offset.toString() << endl;
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	size_t number_of_stations = vecMeteo.size();
	if (number_of_stations == 0) return; //Nothing to write

	size_t number_of_records = vecMeteo[0].size();

	string filename("");
	cfg.getValue("METEOFILE", "Output", filename);

	int ncid, did_time, vid_time, did_points;
	bool create_time = false, create_points = false, create_locations = false, create_variables = false;
	bool exists = IOUtils::fileExists(filename);

	double* dates;
	map<string, double*> map_data;
	map_data[IOUtils::strToUpper(NetCDFIO::lat_str)] = new double[number_of_stations];
	map_data[IOUtils::strToUpper(NetCDFIO::lon_str)] = new double[number_of_stations];
	map_data[cnrm_altitude] = new double[number_of_stations];
	map_data[cnrm_aspect] = new double[number_of_stations];
	map_data[cnrm_slope] = new double[number_of_stations];

	map<string, int> varid;
	map<size_t, string> map_param_name;

	get_parameters(vecMeteo, map_param_name, map_data, dates);

	if (exists) {
		open_file(filename, NC_WRITE, ncid);
		start_definitions(filename, ncid);
	} else {
		create_file(filename, NC_CLASSIC_MODEL, ncid);
		create_time = create_points = create_locations = create_variables = true;
	}

	if (create_time) create_time_dimension(ncid, did_time, vid_time);
	if (create_points) add_dimension(ncid, cnrm_points, number_of_stations, did_points);
	if (create_locations) create_meta_data(ncid, did_points, map_data, varid);
	if (create_variables) create_parameters(ncid, did_time, did_points, number_of_records, number_of_stations, map_param_name, map_data, varid);

	end_definitions(filename, ncid);

	copy_data(number_of_stations, number_of_records, vecMeteo, map_param_name, map_data);

	write_record(ncid, NetCDFIO::cf_time, vid_time, number_of_records, dates);
	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); it++) {
		const string& varname = it->first;
		write_data(ncid, varname, varid[varname], map_data[varname]);
		delete[] it->second;
	}

	close_file(filename, ncid);

	delete[] dates;
}

void NetCDFIO::copy_data(const size_t& number_of_stations, const size_t& number_of_records, const std::vector< std::vector<MeteoData> >& vecMeteo,
                         const std::map<size_t, std::string>& map_param_name, std::map<std::string, double*>& map_data_2D)
{
	for (map<size_t, string>::const_iterator it = map_param_name.begin(); it != map_param_name.end(); it++) {
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
					// do nothing, rely on the _FillValue
				} else if (simple_copy) {
					data[jj*number_of_stations + ii] = value;
				} else if (multiply_copy) {
					data[jj*number_of_stations + ii] = value * multiplier;
				}
			} 
		}
	}
}

void NetCDFIO::create_meta_data(const int& ncid, const int& did, std::map<std::string, double*>& map_data_1D, std::map<std::string, int>& varid)
{
	for (map<string, double*>::const_iterator it = map_data_1D.begin(); it != map_data_1D.end(); it++) {
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
			it++;
		} else {
			map_param_name.erase(it++);
		}
	}	
}

void NetCDFIO::get_parameters(const std::vector< std::vector<MeteoData> >& vecMeteo, std::map<size_t, std::string>& map_param_name, std::map<std::string, double*>& map_data_1D, double*& dates)
{
	size_t number_of_records = vecMeteo[0].size();
	dates = new double[number_of_records];

	double interval = 0;
	for (size_t ii=0; ii<number_of_records; ii++) {
		dates[ii] = vecMeteo[0][ii].date.getModifiedJulianDate();

		if (ii == 1) interval = round((dates[ii] - dates[ii-1]) * 86400.);
	}

	size_t nr_of_parameters = 0;
	if (vecMeteo[0].size() > 0) nr_of_parameters = vecMeteo[0][0].getNrOfParameters();

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
				map_data_1D[IOUtils::strToUpper(NetCDFIO::lat_str)][ii] = meteo_data.meta.position.getLat();
				map_data_1D[IOUtils::strToUpper(NetCDFIO::lon_str)][ii] = meteo_data.meta.position.getLon();
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
		if (check_dim_var(ncid, NetCDFIO::cf_time)) {
			get_dimension(ncid, NetCDFIO::cf_time, did_time);
			get_variable(ncid, NetCDFIO::cf_time, vid_time);			
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
		add_3D_variable(ncid, varname, NC_DOUBLE, did_time, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, varname);
	}

	end_definitions(filename, ncid);

	if (create_dimensions) {
		write_data(ncid, NetCDFIO::lat_str, vid_lat, lat_array);
		write_data(ncid, NetCDFIO::lon_str, vid_lon, lon_array);
	}

	size_t pos_start = append_record(ncid, NetCDFIO::cf_time, vid_time, date.getModifiedJulianDate());
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
	add_dimension(ncid, NetCDFIO::lat_str, grid_in.nrows, did_lat);
	add_1D_variable(ncid, NetCDFIO::lat_str, NC_DOUBLE, did_lat, vid_lat);
	add_attributes_for_variable(ncid, vid_lat, NetCDFIO::lat_str);

	add_dimension(ncid, NetCDFIO::lon_str, grid_in.ncols, did_lon);
	add_1D_variable(ncid, NetCDFIO::lon_str, NC_DOUBLE, did_lon, vid_lon);
	add_attributes_for_variable(ncid, vid_lon, NetCDFIO::lon_str);
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

void NetCDFIO::add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname)
{
	if (varname == NetCDFIO::lat_str) {
		add_attribute(ncid, varid, "standard_name", "latitude");
		add_attribute(ncid, varid, "long_name", "latitude");
		add_attribute(ncid, varid, "units", "degrees_north");
	} else if (varname == NetCDFIO::lon_str) {
		add_attribute(ncid, varid, "standard_name", "longitude");
		add_attribute(ncid, varid, "long_name", "longitude");
		add_attribute(ncid, varid, "units", "degrees_east");
	} else if (varname == NetCDFIO::z_str) {
		add_attribute(ncid, varid, "standard_name", "altitude");
		add_attribute(ncid, varid, "long_name", "height above mean sea level");
		add_attribute(ncid, varid, "units", "m");
		add_attribute(ncid, varid, "positive", "up");
		add_attribute(ncid, varid, "axis", "Z");
	} else if (varname == NetCDFIO::cf_time) {
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

void NetCDFIO::read_data_2D(const int& ncid, const std::string& varname, const int& varid,
                            const size_t& record, const size_t& nr_of_records, const size_t& length, double*& data)
{
	size_t start[] = {record, 0};
     size_t count[] = {nr_of_records, length};

	int status = nc_get_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);	
}

void NetCDFIO::read_value(const int& ncid, const std::string& varname, const int& varid, double& data)
{
	size_t index[] = {0};

	int status = nc_get_var1_double(ncid, varid, index, &data);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve data for variable '" + varname + "': " + nc_strerror(status), AT);	

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

void NetCDFIO::write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& length, double*& data)
{
	size_t start[] = {0};
     size_t count[] = {length};

	int status = nc_put_vara_double(ncid, varid, start, count, data);
	if (status != NC_NOERR)
		throw IOException("Could not write data for variable '" + varname + "': " + nc_strerror(status), AT);	
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

void NetCDFIO::add_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid)
{
	int status = nc_def_dim(ncid, dimname.c_str(), length, &dimid);
	if (status != NC_NOERR)
		throw IOException("Could not define dimension '" + dimname + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value)
{
	int status = nc_put_att_double(ncid, varid, attr_name.c_str(), NC_DOUBLE, 1, &attr_value);
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value)
{
	int status = nc_put_att_text(ncid, varid, attr_name.c_str(), attr_value.size(), attr_value.c_str());
	if (status != NC_NOERR)
		throw IOException("Could not add attribute '" + attr_name + "': " + nc_strerror(status), AT);
}

void NetCDFIO::add_0D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, int& varid)
{
	int dimid;
	int status = nc_def_var(ncid, varname.c_str(), xtype, 0, &dimid, &varid);
	if (status != NC_NOERR)
		throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);
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

void NetCDFIO::add_3D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record, const int& dimid1, const int& dimid2, int& varid)
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
