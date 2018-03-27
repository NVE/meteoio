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
#include <meteoio/plugins/NetCDFIO.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/FileUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/plugins/libncpp.h>
#include <meteoio/meteoStats/libresampling2D.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

namespace mio {
//helper function to sort the cache of grid files
inline bool sort_cache_grids(const std::pair<std::pair<Date,Date>,ncParameters> &left, const std::pair<std::pair<Date,Date>,ncParameters> &right) {
	if (left.first.first < right.first.first) return true;
	if (left.first.first > right.first.first) return false;
	return left.first.second < right.first.second; //date_start equallity case
}

NetCDFIO::NetCDFIO(const std::string& configfile) 
         : cfg(configfile), cache_grid_files(), available_params(), in_schema("ECMWF"), out_schema("ECMWF"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), out_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), dem_altimeter(false), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) 
         : cfg(cfgreader), cache_grid_files(), available_params(), in_schema("ECMWF"), out_schema("ECMWF"), in_grid2d_path(), in_nc_ext(".nc"), out_grid2d_path(), out_file(), in_dflt_TZ(0.), out_dflt_TZ(0.), dem_altimeter(false), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

void NetCDFIO::parseInputOutputSection()
{
	std::string in_grid2d, out_grid2d;
	cfg.getValue("GRID2D", "Input", in_grid2d, IOUtils::nothrow);
	cfg.getValue("GRID2D", "Output", out_grid2d, IOUtils::nothrow);
	if (in_grid2d=="NETCDF") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
		cfg.getValue("NETCDF_SCHEMA", "Input", in_schema, IOUtils::nothrow); IOUtils::toUpper(in_schema);
		cfg.getValue("GRID2DPATH", "Input", in_grid2d_path);
		cfg.getValue("NC_EXT", "INPUT", in_nc_ext, IOUtils::nothrow);
		cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow);
		cfg.getValue("NC_DEBUG", "INPUT", debug, IOUtils::nothrow);
	}
	
	if (out_grid2d=="NETCDF") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
		cfg.getValue("NETCDF_SCHEMA", "Output", out_schema, IOUtils::nothrow); IOUtils::toUpper(out_schema);
		cfg.getValue("GRID2DPATH", "Output", out_grid2d_path);
		cfg.getValue("NC_FILE", "Output", out_file);
	}
}

void NetCDFIO::scanMeteoPath(const std::string& in_path, const std::string& nc_ext, std::vector< std::pair<std::pair<Date,Date>, ncParameters> > &meteo_files)
{
	meteo_files.clear();
	std::list<std::string> dirlist( FileUtils::readDirectory(in_path, nc_ext) );
	if (dirlist.empty()) return; //nothing to do if the directory is empty, we will transparently swap to using GRID2DFILE
	dirlist.sort();

	//Check date range in every filename and cache it
	std::list<std::string>::const_iterator it = dirlist.begin();
	while ((it != dirlist.end())) {
		const std::string filename( in_path + "/" + *it );
		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		const ncParameters ncFile = ncParameters(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		meteo_files.push_back( make_pair(ncFile.getDateRange(), ncFile) );
		it++;
	}

	std::sort(meteo_files.begin(), meteo_files.end(), &sort_cache_grids);

	//now handle overlaping files: truncate the end date of the file starting earlier
	for (size_t ii=0; ii<(meteo_files.size()-1); ii++) {
		if (meteo_files[ii].first.second > meteo_files[ii+1].first.first)
			meteo_files[ii].first.second = meteo_files[ii+1].first.first;
	}
}

bool NetCDFIO::list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& list)
{
	if (cache_grid_files.empty()) scanMeteoPath(in_grid2d_path, in_nc_ext, cache_grid_files);
	if (cache_grid_files.empty()) return true; //there are no grids to read
	
	//HACK handle the case of file_start & file_end are undef() (example: DEM)
	for (size_t ii=0; ii<cache_grid_files.size(); ii++) {
		const Date file_start( cache_grid_files[ii].first.first );
		const Date file_end( cache_grid_files[ii].first.second );
		
		if (file_start > end) return true; //no more files to process (since the files are sorted in cache_grid_files)
		if (file_end < start) continue;
		
		//we consider that the exact same parameters are available at all time steps in the current file
		const std::set<size_t> params_set( cache_grid_files[ii].second.getParams() );
		const std::vector<Date> ts( cache_grid_files[ii].second.getTimestamps() );
		for (size_t jj=0; jj<ts.size(); jj++) {
			if (ts[jj]>end) break; //no more timestamps in the right range
			if (ts[jj]<start) continue;
			list[ ts[jj] ] = params_set;
		}
	}

	return true;
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const std::string& arguments)
{
	std::vector<std::string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() == 2) {
		const ncParameters ncFile(vec_argument[0], ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		grid_out = ncFile.read2DGrid(vec_argument[1]);
	} else {
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::read2DGrid is filename:varname", AT);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (cache_grid_files.empty()) scanMeteoPath(in_grid2d_path, in_nc_ext, cache_grid_files);
	
	if (!cache_grid_files.empty()) {
		for (size_t ii=0; ii<cache_grid_files.size(); ii++) {
			const Date file_start( cache_grid_files[ii].first.first );
			const Date file_end( cache_grid_files[ii].first.second );
			if (file_start > date) return;
			if (file_end < date) continue;
			
			if (date>=file_start && date<=file_end) {
				grid_out = cache_grid_files[ii].second.read2DGrid(parameter, date);
				return;
			}
		}
		//the date was not found
		throw InvalidArgumentException("No Gridded data found for "+date.toString(Date::ISO)+" in '"+in_grid2d_path+"'", AT);
	} else {
		const std::string filename = cfg.get("GRID2DFILE", "Input");
		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		const ncParameters ncFile(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
		grid_out = ncFile.read2DGrid(parameter, date);
	}
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	const std::string filename = cfg.get("DEMFILE", "Input");
	const std::string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	const ncParameters ncFile(filename, ncParameters::READ, cfg, in_schema, in_dflt_TZ, debug);
	const Grid2DObject grid = (varname.empty())? ncFile.readDEM() : ncFile.read2DGrid(varname);
	dem_out = DEMObject( grid ); //we can not directly assign a Grid2DObject to a DEMObject
}

void NetCDFIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& arguments)
{
	throw IOException("Not implemented yet!", AT);
	// arguments is a string of the format filname:varname
	std::vector<std::string> vec_argument;
	if (IOUtils::readLineToVec(arguments, vec_argument, ':')  != 2)
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::write2DGrid is filename:varname", AT);

	const ncParameters::var_attr attr(-1, vec_argument[1], IOUtils::nodata);
	ncParameters::nc_variable tmp_var(attr);
	//write2DGrid(grid_in, vec_argument[0], tmp_var);
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const std::string file_and_path( out_grid2d_path + out_file );
	ncParameters ncFile(file_and_path, ncParameters::WRITE, cfg, out_schema, out_dflt_TZ, debug);
	if (parameter==MeteoGrids::DEM || parameter==MeteoGrids::SHADE || parameter==MeteoGrids::SLOPE || parameter==MeteoGrids::AZI)
		ncFile.write2DGrid(grid_in, parameter, Date()); //do not assign a date to a DEM?
	else
		ncFile.write2DGrid(grid_in, parameter, date);
}


///////////////////////////////////////////////////// Now the ncParameters class starts //////////////////////////////////////////
std::vector<std::string> ncParameters::dimnames( initDimensionNames() );
std::map< std::string, std::vector<ncParameters::nc_dimension> > ncParameters::schemas_dims( initSchemasDims() );
std::map< std::string, std::vector<ncParameters::var_attr> > ncParameters::schemas_vars( initSchemasVars() );

std::vector<std::string> ncParameters::initDimensionNames()
{
	//the order must be the same as in the enum Dimensions
	std::vector<std::string> tmp;
	tmp.push_back("NONE"); tmp.push_back("TIME"); 
	tmp.push_back("LATITUDE"); tmp.push_back("LONGITUDE"); tmp.push_back("ALTITUDE");
	tmp.push_back("NORTHING"); tmp.push_back("EASTING");
	
	return tmp;
}

std::map< std::string, std::vector<ncParameters::nc_dimension> > ncParameters::initSchemasDims()
{
	std::map< std::string, std::vector<ncParameters::nc_dimension> > results;
	std::vector<ncParameters::nc_dimension> tmp;
	
	//CF1 schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "time") );
	tmp.push_back( nc_dimension(LATITUDE, "lat") );
	tmp.push_back( nc_dimension(LONGITUDE, "lon") );
	results["CF1"] = tmp;
	
	//CNRM schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "time") );
	tmp.push_back( nc_dimension(LATITUDE, "latitude") );
	tmp.push_back( nc_dimension(LONGITUDE, "longitude") );
	results["CNRM"] = tmp;
	
	//ECMWF schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "time") );
	tmp.push_back( nc_dimension(LATITUDE, "latitude") );
	tmp.push_back( nc_dimension(LONGITUDE, "longitude") );
	results["ECMWF"] = tmp;
	
	//WRF schema
	tmp.clear();
	tmp.push_back( nc_dimension(TIME, "Times") );
	tmp.push_back( nc_dimension(LATITUDE, "south_north") );
	tmp.push_back( nc_dimension(LONGITUDE, "west_east") );
	results["WRF"] = tmp;
	
	return results;
}

std::map< std::string, std::vector<ncParameters::var_attr> > ncParameters::initSchemasVars()
{
	std::map< std::string, std::vector<ncParameters::var_attr> > results;
	std::vector<ncParameters::var_attr> tmp;

	//CF1 schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "time", "time", "time", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "lat", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "lon", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "height", "altitude", "height above mean sea level", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "temperature", "air_temperature", "near surface air temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RH, "humidity", "relative humidity", "relative humidity", "fraction", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "pressure", "air_pressure", "near surface air pressure", "Pa", IOUtils::nodata) );
	results["CF1"] = tmp;

	//CNRM schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "time", "time", "time", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "ZS", "", "altitude", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::SLOPE, "slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::AZI, "aspect", "", "slope aspect", "degrees from north", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "Tair", "", "Near Surface Air Temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RH, "HUMREL", "", "Relative Humidity", "%", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::QI, "Qair", "", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::VW, "Wind", "", "Wind Speed", "m/s", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DW, "Wind_DIR", "", "Wind Direction", "deg", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::PSUM_L, "Rainf", "", "Rainfall Rate", "kg/m2/s", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::PSUM_S, "Snowf", "", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR_DIR, "DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR_DIFF, "SCA_SWdown", "", "", "", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "PSurf", "", "Surface Pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "LWdown", "", "Surface Incident Longwave Radiation", "W/m2", IOUtils::nodata) );
	results["CNRM"] = tmp;

	//ECMWF schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "time", "time", "time", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "latitude", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "longitude", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "z", "geopotential_height", "geopotential_height", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "t2m", "", "2 metre temperature", "K", 2.) );
	tmp.push_back( var_attr(MeteoGrids::TD, "d2m", "", "2 metre dewpoint temperature", "K", 2.) );
	tmp.push_back( var_attr(MeteoGrids::P, "sp", "surface_air_pressure", "Surface pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P_SEA, "msl", "air_pressure_at_sea_level", "Mean sea level pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ISWR, "ssrd", "surface_downwelling_shortwave_flux_in_air", "Surface solar radiation downwards", "J m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "strd", "", "Surface thermal radiation downwards", "J m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::PSUM, "tp", "", "Total precipitation", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::U, "u10", "", "10 metre U wind component", "m s**-1", 10.) );
	tmp.push_back( var_attr(MeteoGrids::V, "v10", "", "10 metre V wind component", "m s**-1", 10.) );
	tmp.push_back( var_attr(MeteoGrids::SWE, "sd", "lwe_thickness_of_surface_snow_amount", "Snow depth", "m of water equivalent", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TSS, "skt", "", "Skin temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TSG, "stl1", "surface_temperature", "Soil temperature level 1", "K", IOUtils::nodata) ); //this is from 0 to -7cm
	tmp.push_back( var_attr(MeteoGrids::ALB, "al", "surface_albedo", "Albedo", "(0 - 1)", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ALB, "fal", "", "Forecast albedo", "(0 - 1)", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RSNO, "rsn", "", "Snow density", "kg m**-3", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ROT, "ro", "", "Runoff", "m", IOUtils::nodata) );
	results["ECMWF"] = tmp;

	//WRF schema
	tmp.clear();
	tmp.push_back( var_attr(TIME, "Times", "Times", "Times", "", IOUtils::nodata) );
	tmp.push_back( var_attr(LATITUDE, "XLAT", "latitude", "latitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(LONGITUDE, "XLONG", "longitude", "longitude", "degrees", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::DEM, "HGT", "Terrain Height", "Terrain Height", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "PSFC", "Surface pressure", "Surface pressure", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "T2", "2-meter temperature", "2-meter temperature", "K", 2.) );
	tmp.push_back( var_attr(MeteoGrids::QI, "Q2", "2-meter specific humidity", "2-meter specific humidity", "kg kg-1", 2) );
	tmp.push_back( var_attr(MeteoGrids::ISWR, "ACSWDNB", "Downward SW surface radiation", "Downward SW surface radiation", "W m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RSWR, "ACSWUPB", "Upwelling Surface Shortwave Radiation", "Upwelling Surface Shortwave Radiation", "W m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ILWR, "ACLWDNB", "Downward LW surface radiation", "Downward LW surface radiation", "W m**-2", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::ROT, "SFROFF", "Surface runoff ", "Surface runoff ", "kg*m2*s-1", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::HS, "SNOWH", "Snow depth", "Snow depth", "Pa", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TSS, "TSK", "Surface skin temperature", "Surface skin temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::U, "U10", "10-meter wind speed", "10 metre U wind component", "m s**-1", 10.) );
	tmp.push_back( var_attr(MeteoGrids::V, "V10", "10-meter wind speed", "10 metre V wind component", "m s**-1", 10.) );
	results["WRF"] = tmp;
	
	return results;
}

//The user can provide his own variables properties as NETCDF::{param} = {name}
std::vector<ncParameters::var_attr> ncParameters::initUserSchemas(const Config& i_cfg)
{
	std::vector<ncParameters::var_attr> results;
	
	//TODO allow the user to provide the units, long name, etc?
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string meteo_grid( custom_attr[ii].substr(found+1) );
		const std::string netcdf_param = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = getParameterIndex(meteo_grid);
		if (param_index==IOUtils::npos)
			throw InvalidArgumentException("Parameter '"+meteo_grid+"' is not a valid MeteoGrid! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( var_attr(param_index, netcdf_param, IOUtils::nodata) );
	}
	
	return results;
}

//The user can provide his own variables properties as NETCDF_DIM::{dimension_param} = {name_in_current_file}
std::vector<ncParameters::nc_dimension> ncParameters::initUserDimensions(const Config& i_cfg)
{
	std::vector<ncParameters::nc_dimension> results;
	
	//TODO allow the user to provide the units, long name, etc?
	const std::vector<std::string> custom_attr( i_cfg.getKeys("NETCDF_DIM::", "Input") );
	const size_t nrOfCustoms = custom_attr.size();
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string dim_str( custom_attr[ii].substr(found+1) );
		const std::string netcdf_dim = i_cfg.get(custom_attr[ii], "Input");
		const size_t param_index = getParameterIndex(dim_str);
		if (param_index==IOUtils::npos || param_index<firstdimension || param_index>lastdimension)
			throw InvalidArgumentException("Dimension '"+dim_str+"' is not a valid dimension! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( nc_dimension( static_cast<Dimensions>(param_index), netcdf_dim) );
	}
	
	return results;
}

ncParameters::ncParameters(const std::string& filename, const Mode& mode, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug)
             : user_schemas( initUserSchemas(cfg) ), user_dimensions( initUserDimensions(cfg) ), vars(), unknown_vars(), vecTime(), vecLat(), vecLon(), dimensions_map(), file_and_path(filename), coordin(), coordinparam(), TZ(tz_in), wrf_hacks(schema=="WRF"), debug(i_debug)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	
	if (mode==WRITE) {
		initFromSchema(schema);
		if (FileUtils::fileExists(filename)) initFromFile(filename, schema);
	} else if (mode==READ)
		initFromFile(filename, schema);
	
	if (debug) {
		std::cout << filename << ":\n";
		std::cout << "\tDimensions:\n";
		for (std::map<size_t, nc_dimension>::const_iterator it = dimensions_map.begin(); it!=dimensions_map.end(); ++it)
			std::cout << "\t\t" << it->second.toString() << "\n";
		if (!vecTime.empty()) std::cout << "\ttime range: [" << vecTime.front().toString(Date::ISO) << " - " << vecTime.back().toString(Date::ISO) << "]\n";
		std::cout << "\tVariables:\n";
		for (std::map<size_t, nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
			std::cout << "\t\t" << getParameterName( it->first ) << " -> " << it->second.toString() << "\n";
		std::cout << "\tUnrecognized variables:\n";
		for (std::map<std::string, nc_variable>::const_iterator it=unknown_vars.begin(); it!=unknown_vars.end(); ++it)
			std::cout << "\t\t" << it->first << " -> " << it->second.toString() << "\n";
	}
}

//populate the dimensions_map from the selected schema
void ncParameters::initFromSchema(const std::string& schema)
{
	for (size_t ii=0; ii<schemas_dims[schema].size(); ii++) {
		dimensions_map[ schemas_dims[schema][ii].type ] = schemas_dims[schema][ii];
	}
	if (dimensions_map.count(TIME)==0) throw IOException("No TIME dimension in schema '"+schema+"'", AT);
	dimensions_map[ TIME ].isUnlimited = true;
	
	for (size_t ii=0; ii<schemas_vars[schema].size(); ii++) {
		vars[ schemas_vars[schema][ii].param ] = nc_variable( schemas_vars[schema][ii] );
	}
}

//populate the dimensions_map and vars and unknown_vars from the file
void ncParameters::initFromFile(const std::string& filename, const std::string& schema)
{
	if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
	
	int ncid;
	ncpp::open_file(filename, NC_NOWRITE, ncid);
	
	//read the dimensions and variables
	initDimensionsFromFile(ncid, schema);
	initVariablesFromFile(ncid, schema);
	
	if (dimensions_map.count(TIME)!=0) vecTime = read_1Dvariable(ncid);
	if (dimensions_map.count(LATITUDE)!=0) vecLat = read_1Dvariable(ncid, LATITUDE);
	if (dimensions_map.count(LONGITUDE)!=0) vecLon = read_1Dvariable(ncid, LONGITUDE);
	
	ncpp::close_file(filename, ncid);
}

std::pair<Date, Date> ncParameters::getDateRange() const
{
	if (vecTime.empty()) return make_pair( Date(), Date() );
	return make_pair( vecTime.front(), vecTime.back() );
}

std::set<size_t> ncParameters::getParams() const 
{
	std::set<size_t> available_params;
	for (std::map<size_t, nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
		available_params.insert( it->first );
	
	return available_params;
}

Grid2DObject ncParameters::readDEM() const
{
	if (vars.count(MeteoGrids::DEM)!=0) 
		return read2DGrid(MeteoGrids::DEM, Date());
	
	for (std::map<std::string, nc_variable>::const_iterator it = unknown_vars.begin(); it!=unknown_vars.end(); ++it) {
		if (it->first=="Band1") return read2DGrid(it->second, IOUtils::npos); //ASTER naming
		if (it->first=="z") return read2DGrid(it->second, IOUtils::npos); //GDAL naming
		if (it->first=="height") return read2DGrid(it->second, IOUtils::npos); //MeteoCH naming
		if (it->first=="HGT") return read2DGrid(it->second, IOUtils::npos); //WRF naming
	}
	
	throw NotFoundException("No DEM could not be found in file '"+file_and_path+"'", AT);
}

Grid2DObject ncParameters::read2DGrid(const std::string& varname) const
{
	for (std::map<size_t, nc_variable>::const_iterator it = vars.begin(); it!=vars.end(); ++it) {
		if (it->second.attributes.name==varname)
			return read2DGrid(it->second, IOUtils::npos);
	}
	
	for (std::map<std::string, nc_variable>::const_iterator it = unknown_vars.begin(); it!=unknown_vars.end(); ++it) {
		if (it->first==varname)
			return read2DGrid(it->second, IOUtils::npos);
	}
	
	throw NotFoundException("The variable '"+varname+"' could not be found in file '"+file_and_path+"'", AT);
}

Grid2DObject ncParameters::read2DGrid(const size_t& param, const Date& date) const
{
	const std::map <size_t, nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end()) NoDataException("No "+MeteoGrids::getParameterName( param )+" grid in file "+file_and_path, AT);
	
	size_t time_pos = IOUtils::npos;
	if (!date.isUndef()) {
		const std::vector<Date>::const_iterator low = std::lower_bound(vecTime.begin(), vecTime.end(), date);
		if (*low!=date) throw NoDataException("No data at "+date.toString(Date::ISO)+" in file "+file_and_path, AT);
		time_pos = static_cast<size_t>( std::distance(vecTime.begin(), low) );
	} else { 
		//no date has been provided, check if this parameter depends on time
		const std::map<size_t, nc_dimension>::const_iterator it2 = dimensions_map.find(TIME);
		const int time_id = (it2!=dimensions_map.end())? it2->second.dimid : -1;
		const bool depend_on_time = (std::find(it->second.dimids.begin(), it->second.dimids.end(), time_id) != it->second.dimids.end());
		if (depend_on_time && !vecTime.empty() && (param!=MeteoGrids::DEM && param!=MeteoGrids::SLOPE && param!=MeteoGrids::AZI))
			throw InvalidFormatException("No time requirement has been provided for a file that contains multiple timestamps", AT);
	}
	
	const bool isPrecip = (param==MeteoGrids::PSUM || param==MeteoGrids::PSUM_L || param==MeteoGrids::PSUM_S);
	const bool isRad = (param==MeteoGrids::ISWR || param==MeteoGrids::RSWR || param==MeteoGrids::ISWR_DIFF || param==MeteoGrids::ISWR_DIR);
	return read2DGrid(it->second, time_pos, isPrecip, (isPrecip || isRad));
}

Grid2DObject ncParameters::read2DGrid(const nc_variable& var, const size_t& time_pos, const bool& m2mm, const bool& reZero) const
{
	//define the results grid
	mio::Coords llcorner(coordin, coordinparam);
	llcorner.setLatLon( std::min(vecLat.front(), vecLat.back()), std::min(vecLon.front(), vecLon.back()), IOUtils::nodata);
	double resampling_factor_x = IOUtils::nodata, resampling_factor_y=IOUtils::nodata;
	const double cellsize = calculate_cellsize(resampling_factor_x, resampling_factor_y); //HACK expand the definition of Grid2DObject to support lat/lon grids and reproject in GridsManager
	Grid2DObject grid(vecLon.size(), vecLat.size(), cellsize, llcorner);
	
	//read the raw data, copy it into the Grid2DObject
	int ncid;
	ncpp::open_file(file_and_path, NC_NOWRITE, ncid);
	
	double *data = new double[vecLat.size()*vecLon.size()]; //HACK do it with variable/dimensions dependencies
	if (time_pos!=IOUtils::npos)
		ncpp::read_data(ncid, var.attributes.name, var.varid, time_pos, vecLat.size(), vecLon.size(), data);
	else
		ncpp::read_data(ncid, var.attributes.name, var.varid, data);
	fill2DGrid(grid, data, var.nodata);
	delete[] data;
	ncpp::close_file(file_and_path, ncid);
	
	//handle data packing and units, if necessary
	if (var.scale!=1.) grid *= var.scale;
	if (var.offset!=0.) grid += var.offset;
	const std::string units( var.attributes.units );
	if (!units.empty()) {
		if (units=="m**2 s**-2") grid /= Cst::gravity;
		else if (units=="%") grid /= 100.;
		else if (units=="J m**-2") grid /= (3600.*3.); //HACK: this is the ECMWF sampling rate
		else if (m2mm && units=="m") grid *= 1000.;
		
		if (reZero) {//reset very low values to zero
			for (size_t ii=0; ii<grid.size(); ii++)
				if (grid(ii)<1e-6 && grid(ii)!=mio::IOUtils::nodata) grid(ii)=0.;
		}
	}
	
	return grid;
}

void ncParameters::write2DGrid(Grid2DObject grid_in, const size_t& param, const Date& date)
{
	if (vars.count(param)>0) {
		write2DGrid(grid_in, vars[param], date);
	} else {
		const std::string param_name( MeteoGrids::getParameterName(param) );
		const var_attr tmp_attr(param, param_name, IOUtils::nodata);
		nc_variable tmp_var(tmp_attr);
		write2DGrid(grid_in, tmp_var, date);
	}
}

void ncParameters::write2DGrid(Grid2DObject grid_in, nc_variable& var, const Date& date)
{
	const bool is_record = (!date.isUndef());
	
	bool create_spatial_dimensions(false), create_variable(false), create_time(false);
	int ncid;
	if ( FileUtils::fileExists(file_and_path) ) {
		ncpp::open_file(file_and_path, NC_WRITE, ncid);
		//HACK chech "somewhere" that all dimensions / variables (TIME, LAT, LON) are available in dimensions_map (ie once and for all)
		if (is_record) create_time = (dimensions_map[TIME].dimid == -1);
		const bool hasLatitude = (dimensions_map[LATITUDE].dimid != -1) && (dimensions_map[LATITUDE].length == grid_in.getNy());
		const bool hasLongitude = (dimensions_map[LONGITUDE].dimid != -1) && (dimensions_map[LONGITUDE].length == grid_in.getNy());
		create_spatial_dimensions = (!hasLatitude && !hasLongitude);
		ncpp::start_definitions(file_and_path, ncid); //call nc_redef
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.6");
		create_variable = create_spatial_dimensions = true;
		if (is_record) create_time = true;
	}
	
	if (create_time) {
		std::string date_str( Date(date.getYear(), 1, 1, 0, 0, TZ).toString(Date::ISO) );
		IOUtils::replace_all(date_str, "T", " ");
		vars[TIME].attributes.units = "hours since " + date_str;
		create_dimension_and_variable(ncid, dimensions_map[TIME], vars[TIME]);
	}
	if (create_spatial_dimensions) {
		dimensions_map[LATITUDE].length = grid_in.getNy();
		create_dimension_and_variable(ncid, dimensions_map[LATITUDE], vars[LATITUDE]);
		
		dimensions_map[LONGITUDE].length = grid_in.getNx();
		create_dimension_and_variable(ncid, dimensions_map[LONGITUDE], vars[LONGITUDE]);
	}
	if (create_variable) {
		if (is_record) var.dimids.push_back(dimensions_map[TIME].dimid);
		var.dimids.push_back(dimensions_map[LATITUDE].dimid);
		var.dimids.push_back(dimensions_map[LONGITUDE].dimid);
		ncParameters::create_variable(ncid, var);
	}
	ncpp::end_definitions(file_and_path, ncid);
	
	//write the dimensions' data
	if (create_spatial_dimensions) {
		double *lat_array = new double[grid_in.getNy()];
		double *lon_array = new double[grid_in.getNx()];
		ncpp::calculate_dimensions(grid_in, lat_array, lon_array);
		ncpp::write_data(ncid, vars[LATITUDE].attributes.name, vars[LATITUDE].varid, lat_array);
		ncpp::write_data(ncid, vars[LONGITUDE].attributes.name, vars[LONGITUDE].varid, lon_array);
		delete[] lat_array; delete[] lon_array; 
	}
	
	//now write the data
	int *data = new int[grid_in.getNy() * grid_in.getNx()];
	ncpp::fill_grid_data(grid_in, IOUtils::nodata, data);
	if (is_record) {
		if (dimensions_map[TIME].length>0) {
			double last_value = IOUtils::nodata;
			const size_t index = dimensions_map[TIME].length - 1;
			const int status = nc_get_var1_double(ncid, vars[TIME].varid, &index, &last_value);
			if (status != NC_NOERR) throw IOException("Could not retrieve data for variable '" + vars[TIME].attributes.name + "': " + nc_strerror(status), AT);

			//const Date last_date( read_variableAtPos(ncid, index) );
			
			//if (last_value == data) return (dimlen - 1); //The timestamp already exists

			/*if (last_value > data) {
				const size_t pos = find_record(ncid, varname, dimid, data); // Search for a possible match

				if (pos != IOUtils::npos) {
					return pos;
				} else {
					throw IOException("The variable '" + varname + "' has to be linearly increasing", AT);
				}
			}*/
			std::cout << "Last time value: " << last_value << "\n";
		}
		
		const size_t start[] = {dimensions_map[TIME].length};
		const size_t count[] = {1};

		double dt = date.getJulian();
		const int status = nc_put_vara_double(ncid, vars[TIME].varid, start, count, &dt);
		if (status != NC_NOERR) throw IOException("Could not write data for record variable '" + vars[TIME].attributes.name + "': " + nc_strerror(status), AT);
		
		
		size_t pos_start = dimensions_map[TIME].length;
		
		ncpp::write_data(ncid, var.attributes.name, var.varid, grid_in.getNy(), grid_in.getNx(), pos_start, data);
	} else {
		ncpp::write_data(ncid, var.attributes.name, var.varid, data);
	}
	delete[] data;
	
	ncpp::close_file(file_and_path, ncid);
}

void ncParameters::create_dimension_and_variable(const int& ncid, nc_dimension &dim, nc_variable& var)
{
	create_dimension(ncid, dim);
	var.dimids.push_back( dim.dimid );
	ncParameters::create_variable(ncid, var);
}

void ncParameters::create_dimension(const int& ncid, nc_dimension &dim)
{
	const nc_type length = (dim.isUnlimited)? NC_UNLIMITED : static_cast<int>(dim.length);
	int status = nc_def_dim(ncid, dim.name.c_str(), length, &dim.dimid);
	if (status != NC_NOERR) throw IOException("Could not define dimension '" + dim.name + "': " + nc_strerror(status), AT);
}

//write the variable's attributes into the file
void ncParameters::create_variable(const int& ncid, nc_variable& var)
{
	const int ndims = static_cast<int>( var.dimids.size() );
	const int status = nc_def_var(ncid, var.attributes.name.c_str(), NC_DOUBLE, ndims, &var.dimids[0], &var.varid);
	if (status != NC_NOERR) throw IOException("Could not define variable '" + var.attributes.name + "': " + nc_strerror(status), AT);
	
	if (!var.attributes.standard_name.empty()) ncpp::add_attribute(ncid, var.varid, "standard_name", var.attributes.standard_name);
	if (!var.attributes.long_name.empty()) ncpp::add_attribute(ncid, var.varid, "long_name", var.attributes.long_name);
	if (!var.attributes.units.empty()) ncpp::add_attribute(ncid, var.varid, "units", var.attributes.units);
	
	if (var.attributes.param==TIME) ncpp::add_attribute(ncid, var.varid, "calendar", "gregorian");
	if (var.attributes.param==MeteoGrids::DEM) {
		ncpp::add_attribute(ncid, var.varid, "positive", "up");
		ncpp::add_attribute(ncid, var.varid, "axis", "Z");
	}
}

void ncParameters::initDimensionsFromFile(const int& ncid, const std::string& schema_name)
{
	int status;
	int ndims;
	status = nc_inq_ndims(ncid, &ndims);
	if (status != NC_NOERR) throw IOException("Could not retrieve number of dimensions: " + std::string(nc_strerror(status)), AT);
	
	int unlim_id;
	status = nc_inq_unlimdim(ncid, &unlim_id);
	if (status != NC_NOERR) throw IOException("Could not retrieve unlimited dimension: " + std::string(nc_strerror(status)), AT);
	
	int *dimids = (int*)malloc(ndims * sizeof(int));
	status = nc_inq_dimids(ncid, &ndims, dimids, 0);
	if (status != NC_NOERR) throw IOException("Could not retrieve dimensions IDs: " + std::string(nc_strerror(status)), AT);
	
	for (int idx=0; idx<ndims; idx++) {
		char name[NC_MAX_NAME+1];
		status = nc_inq_dimname(ncid, dimids[idx], name);
		if (status != NC_NOERR) throw IOException("Could not retrieve dimension name: " + std::string(nc_strerror(status)), AT);
		
		nc_dimension tmp_dim( getSchemaDimension(name, schema_name) ); //set name and type
		tmp_dim.dimid = idx;
		tmp_dim.isUnlimited = (idx==unlim_id);
		status = nc_inq_dimlen(ncid, dimids[idx], &tmp_dim.length);
		if (status != NC_NOERR) throw IOException("Could not retrieve dimension lenght: " + std::string(nc_strerror(status)), AT);
		
		dimensions_map[ tmp_dim.type ] = tmp_dim;
	}
	
	free( dimids );
}

void ncParameters::initVariablesFromFile(const int& ncid, const std::string& schema_name)
{
	int nr_of_variables = -1;
	int status = nc_inq_nvars(ncid, &nr_of_variables);
	if (status != NC_NOERR)
		throw IOException("Could not retrieve variables for dataset: " + string(nc_strerror(status)), AT);
	
	// Variable IDs in a NetCDF file are consecutive integers starting with 0
	for (int ii=0; ii<nr_of_variables; ++ii) {
		int nrdims, varid = -1;
		int dimids[NC_MAX_VAR_DIMS];
		char name[NC_MAX_NAME+1];
		
		status = nc_inq_varname(ncid, ii, name);
		status = nc_inq_varid (ncid, name, &varid);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
		status = nc_inq_var(ncid, varid, NULL, NULL, &nrdims, dimids, NULL);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
		
		const std::string varname( name );
		nc_variable tmp_var(getSchemaAttributes(varname, schema_name), 1., 0., IOUtils::nodata, varid);
		getAttribute(ncid, tmp_var.varid, varname, "missing_value", tmp_var.nodata);
		getAttribute(ncid, tmp_var.varid, varname, "scale_factor", tmp_var.scale);
		getAttribute(ncid, tmp_var.varid, varname, "add_offset", tmp_var.offset);
		getAttribute(ncid, tmp_var.varid, varname, "units", tmp_var.attributes.units);
		tmp_var.dimids.assign(dimids, dimids+nrdims);
		
		if (tmp_var.attributes.param!=IOUtils::npos) {
			if (tmp_var.attributes.param==TIME && !wrf_hacks) 
				getTimeTransform(tmp_var.attributes.units, TZ, tmp_var.offset, tmp_var.scale);
			vars[ tmp_var.attributes.param ] = tmp_var;
		} else
			unknown_vars[ varname ] = tmp_var;
	}
}

Date ncParameters::read_variableAtPos(const int& ncid, const size_t pos) const
{
	const std::map<size_t, nc_variable>::const_iterator it = vars.find( TIME );
	if (it==vars.end()) throw InvalidArgumentException("Could not find parameter \"TIME\" in file \""+file_and_path+"\"", AT);
	
	if (!wrf_hacks) {
		double value;
		const int status = nc_get_var1_double(ncid, it->second.varid, &pos, &value);
		if (status != NC_NOERR) throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);
		return Date(value*it->second.scale + it->second.offset, TZ);
	} else {
		static const size_t DateStrLen = 19; //HACK DateStrLen = 19, defined in Dimensions
		char *data = (char*)calloc(1, sizeof(char)*DateStrLen);
		const int status = nc_get_var1_text(ncid, it->second.varid, &pos, data);
		if (status != NC_NOERR) throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);
		std::string tmp( data );
		IOUtils::replace_all(tmp, "_", "T");
		Date result;
		IOUtils::convertString(result, tmp, TZ);
		return result;
	}
}

std::vector<Date> ncParameters::read_1Dvariable(const int& ncid) const
{
	if (!wrf_hacks) {
		const std::vector<double> tmp_results( read_1Dvariable(ncid, TIME) );
		const std::map<size_t, nc_variable>::const_iterator it = vars.find( TIME ); //it exists since it has been read above
		std::vector<Date> results(tmp_results.size());
		for (size_t ii=0; ii<tmp_results.size(); ii++)
			results[ii].setDate(tmp_results[ii]*it->second.scale + it->second.offset, TZ);
		return results;
	} else {
		static const size_t DateStrLen = 19; //HACK DateStrLen = 19, defined in Dimensions
		const std::map<size_t, nc_variable>::const_iterator it = vars.find( TIME );
		if (it==vars.end()) throw InvalidArgumentException("Could not find parameter \"TIME\" in file \""+file_and_path+"\"", AT);
		const size_t length = read_1DvariableLength(it->second);
		
		char *data = (char*)calloc(length, sizeof(char)*DateStrLen);
		const int status = nc_get_var_text(ncid, it->second.varid, data);
		if (status != NC_NOERR) throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);

		std::vector<Date> results(length);
		for(size_t ii=0; ii<length; ii++) {
			std::string tmp(DateStrLen, '\0');
			for(size_t jj=0; jj<DateStrLen; jj++) {
				const char c = data[ii*DateStrLen+jj];
				tmp[jj] = (c!='_')? c : 'T';
			}
			IOUtils::convertString(results[ii], tmp, TZ);
		}
		free( data );
		return results;
	}
}

std::vector<double> ncParameters::read_1Dvariable(const int& ncid, const size_t& param) const
{
	const std::map<size_t, nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end()) throw InvalidArgumentException("Could not find parameter \""+getParameterName(param)+"\" in file \""+file_and_path+"\"", AT);
	const size_t length = read_1DvariableLength(it->second);
	
	std::vector<double> results( length );
	double *data = new double[ length ];
	ncpp::read_data(ncid, it->second.attributes.name, it->second.varid, data);
	std::copy(data, data+length, results.begin());
	delete[] data;
	return results;
}

size_t ncParameters::read_1DvariableLength(const nc_variable& var) const
{
	if (var.dimids.size()!=1) throw InvalidArgumentException("Parameter \""+getParameterName(var.attributes.param)+"\" in file \""+file_and_path+"\" is not a 1D variable", AT);
	
	const int dimid = var.dimids[0];
	std::map<size_t, nc_dimension>::const_iterator it = dimensions_map.begin();
	for (; it!=dimensions_map.end(); ++it) {
		if (it->second.dimid==dimid) break;
	}
	if (it==dimensions_map.end()) throw InvalidArgumentException("Could not find a dimension in file \""+file_and_path+"\"", AT);
	
	return it->second.length;
}

double ncParameters::calculate_cellsize(double& factor_x, double& factor_y) const
{
	//in order to handle swapped llcorner/urcorner, we use "fabs" everywhere
	double alpha;
	const double cntr_lat = .5*fabs(vecLat.front()+vecLat.back());
	const double cntr_lon = .5*fabs(vecLon.front()+vecLon.back());
	const double distanceX = CoordsAlgorithms::VincentyDistance(cntr_lat, vecLon.front(), cntr_lat, vecLon.back(), alpha);
	const double distanceY = CoordsAlgorithms::VincentyDistance(vecLat.front(), cntr_lon, vecLat.back(), cntr_lon, alpha);

	//round to 1cm precision for numerical stability
	const double cellsize_x = static_cast<double>(Optim::round( distanceX / static_cast<double>(vecLon.size())*100. )) / 100.;
	const double cellsize_y = static_cast<double>(Optim::round( distanceY / static_cast<double>(vecLat.size())*100. )) / 100.;
	if (cellsize_x == cellsize_y) {
		return cellsize_x;
	} else {
		const double cellsize = std::min(cellsize_x, cellsize_y);
		factor_x =  cellsize_x / cellsize;
		factor_y =  cellsize_y / cellsize;
		return cellsize;
	}
}

//populate the results grid and handle the case of llcorner/urcorner swapped
void ncParameters::fill2DGrid(Grid2DObject& grid, const double data[], const double& nodata) const
{
	if (vecLat.front()<=vecLat.back()) {
		for (size_t kk=0; kk < vecLat.size(); kk++) {
			const size_t row = kk*vecLon.size();
			if (vecLon.front()<=vecLon.back()) {
				for (size_t ll=0; ll < vecLon.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < vecLon.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + (vecLon.size() -1) - ll], nodata);
			}
		}
	} else {
		for (size_t kk=0; kk < vecLat.size(); kk++) {
			const size_t row = ((vecLat.size()-1) - kk)*vecLon.size();
			if (vecLon.front()<=vecLon.back()) {
				for (size_t ll=0; ll < vecLon.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + ll], nodata);
			} else {
				for (size_t ll=0; ll < vecLon.size(); ll++)
					grid(ll, kk) = mio::IOUtils::standardizeNodata(data[row + (vecLon.size() -1) - ll], nodata);
			}
		}
	}
}

void ncParameters::getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier)
{
	static const double equinox_year = 365.242198781; //definition used by the NetCDF Udunits package
	
	std::vector<std::string> vecString;
	const size_t nrWords = IOUtils::readLineToVec(time_units, vecString);
	if (nrWords<3 || nrWords>4) throw InvalidArgumentException("Invalid format for time units: \'"+time_units+"\'", AT);
	
	if (vecString[0]=="years") o_time_multiplier = equinox_year;
	else if (vecString[0]=="months") o_time_multiplier = equinox_year/12.;
	else if (vecString[0]=="days") o_time_multiplier = 1.;
	else if (vecString[0]=="hours") o_time_multiplier = 1./24.;
	else if (vecString[0]=="minutes") o_time_multiplier = 1./(24.*60.);
	else if (vecString[0]=="seconds") o_time_multiplier = 1./(24.*3600);
	else throw InvalidArgumentException("Unknown time unit \'"+vecString[0]+"\'", AT);
	
	const std::string ref_date_str = (nrWords==3)? vecString[2] : vecString[2]+"T"+vecString[3];
	Date refDate;
	if (!IOUtils::convertString(refDate, ref_date_str, i_TZ))
		throw InvalidArgumentException("Invalid reference date \'"+ref_date_str+"\'", AT);
	
	o_time_offset = refDate.getJulian();
}

const ncParameters::var_attr ncParameters::getSchemaAttributes(const std::string& var, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_schemas.size(); ii++) {
		if (user_schemas[ii].name==var) return user_schemas[ii];
	}
	
	std::map< std::string, std::vector<ncParameters::var_attr> >::const_iterator it = schemas_vars.find( schema_name );
	if (it==schemas_vars.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==var) return it->second[ii];
	}
	
	return var_attr();
}

const ncParameters::nc_dimension ncParameters::getSchemaDimension(const std::string& dimname, const std::string& schema_name) const
{
	//the user defined schema has priority
	for (size_t ii=0; ii<user_dimensions.size(); ii++) {
		if (user_dimensions[ii].name==dimname) return user_dimensions[ii];
	}
	
	std::map< std::string, std::vector<ncParameters::nc_dimension> >::const_iterator it = schemas_dims.find( schema_name );
	if (it==schemas_dims.end())
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema_name+"\"", AT);
	
	for (size_t ii=0; ii<it->second.size(); ii++) {
		if (it->second[ii].name==dimname) return it->second[ii];
	}
	
	return nc_dimension();
}

//if the attribute is not found, an empty string is returned
void ncParameters::getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, std::string& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, value_id, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		char* value = new char[attr_len + 1]; // +1 for trailing null
		status = nc_get_att_text(ncid, value_id, attr_name.c_str(), value);
		if (status != NC_NOERR) throw IOException("Could not read attribute '" + attr_name + "' for '" + value_name + "': " + nc_strerror(status), AT);

		value[attr_len] = '\0';
		attr_value = value;
		delete[] value;
	}
}

//if the attribute is not found, the attr_value is not changed
void ncParameters::getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, double& attr_value)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, value_id, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		status = nc_get_att_double(ncid, value_id, attr_name.c_str(), &attr_value);
		if (status != NC_NOERR) throw IOException("Could not read attribute '" + attr_name + "' for '" + value_name + "': " + nc_strerror(status), AT);
	}
}

//Since we had to extend MeteoGrids::Parameters, we must redefine this method
std::string ncParameters::getParameterName(const size_t& param)
{
	if (param==IOUtils::npos) return "";
	
	if (param>=NONE) {
		if (param>lastdimension) 
			throw IndexOutOfBoundsException("Trying to get name for a dimension that does not exist", AT);
		return dimnames[ param - firstdimension ];
	}
	
	return MeteoGrids::getParameterName( param );
}

//Since we had to extend MeteoGrids::Parameters, we must redefine this method
size_t ncParameters::getParameterIndex(const std::string& param)
{
	for (size_t ii=firstdimension; ii<=lastdimension; ii++) {
		if (dimnames[ii]==param) return ii;
	}
	
	return MeteoGrids::getParameterIndex( param );
}

} //namespace
