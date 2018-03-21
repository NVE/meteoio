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
         : cfg(configfile), cache_grid_files(), available_params(), in_schema(), out_schema(), in_grid2d_path(), in_nc_ext(".nc"), in_dflt_TZ(0.), out_dflt_TZ(0.), dem_altimeter(false), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) 
         : cfg(cfgreader), cache_grid_files(), available_params(), in_schema(), out_schema(), in_grid2d_path(), in_nc_ext(".nc"), in_dflt_TZ(0.), out_dflt_TZ(0.), dem_altimeter(false), debug(false)
{
	//IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

void NetCDFIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow);
	
	cfg.getValue("NETCDF_SCHEMA", "Input", in_schema, IOUtils::nothrow); IOUtils::toUpper(in_schema);
	cfg.getValue("NETCDF_SCHEMA", "Output", out_schema, IOUtils::nothrow); IOUtils::toUpper(out_schema);
	
	cfg.getValue("GRID2DPATH", "Input", in_grid2d_path, IOUtils::nothrow);
	cfg.getValue("NC_EXT", "INPUT", in_nc_ext, IOUtils::nothrow);
	
	cfg.getValue("NC_DEBUG", "INPUT", debug, IOUtils::nothrow);
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
		const ncParameters ncFile = ncParameters(filename, cfg, in_schema, in_dflt_TZ, debug);
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
		const ncParameters ncFile(vec_argument[0], cfg, in_schema, in_dflt_TZ, debug);
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
		const ncParameters ncFile(filename, cfg, in_schema, in_dflt_TZ, debug);
		grid_out = ncFile.read2DGrid(parameter, date);
	}
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	const std::string filename = cfg.get("DEMFILE", "Input");
	const std::string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	const ncParameters ncFile(filename, cfg, in_schema, in_dflt_TZ, debug);
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

	const std::string name( vec_argument[1] );
	/*const attributes attr(name, name, name, "", IOUtils::nodata);
	write2DGrid_internal(grid_in, vec_argument[0], attr);*/
}

void NetCDFIO::write2DGrid(const Grid2DObject& /*grid_in*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	throw IOException("Not implemented yet!", AT);
}


///////////////////////////////////////////////////// Now the ncParameters class starts //////////////////////////////////////////
std::map< std::string, std::vector<ncParameters::dim_attributes> > ncParameters::schemas_dims( initSchemasDims() );
std::map< std::string, std::vector<ncParameters::var_attr> > ncParameters::schemas_vars( initSchemasVars() );

std::map< std::string, std::vector<ncParameters::dim_attributes> > ncParameters::initSchemasDims()
{
	std::map< std::string, std::vector<ncParameters::dim_attributes> > results;
	std::vector<ncParameters::dim_attributes> tmp;
	
	//CF1 schema
	tmp.clear();
	tmp.push_back( dim_attributes(TIME, "time", "time", "years since 1900-01-01 00:00:00") );
	tmp.push_back( dim_attributes(LATITUDE, "lat", "latitude", "degrees") );
	tmp.push_back( dim_attributes(LONGITUDE, "lon", "longitude", "degrees") );
	results["CF1"] = tmp;
	
	//CNRM schema
	tmp.clear();
	tmp.push_back( dim_attributes(TIME, "time", "time", "") );
	tmp.push_back( dim_attributes(LATITUDE, "latitude", "Latitude", "degrees_north") );
	tmp.push_back( dim_attributes(LONGITUDE, "longitude", "Longitude", "degrees_east") );
	results["CNRM"] = tmp;
	
	//ECMWF schema
	tmp.clear();
	tmp.push_back( dim_attributes(TIME, "time", "time", "") );
	tmp.push_back( dim_attributes(LATITUDE, "latitude", "latitude", "degrees_north") );
	tmp.push_back( dim_attributes(LONGITUDE, "longitude", "longitude", "degrees_east") );
	results["ECMWF"] = tmp;
	
	//WRF schema
	tmp.clear();
	tmp.push_back( dim_attributes(TIME, "Times", "Times", "") );
	tmp.push_back( dim_attributes(LATITUDE, "south_north", "latitude", "degrees") );
	tmp.push_back( dim_attributes(LONGITUDE, "west_east", "longitude", "degrees") );
	results["WRF"] = tmp;
	
	return results;
}

std::map< std::string, std::vector<ncParameters::var_attr> > ncParameters::initSchemasVars()
{
	std::map< std::string, std::vector<ncParameters::var_attr> > results;
	std::vector<ncParameters::var_attr> tmp;

	//CF1 schema
	tmp.clear();
	tmp.push_back( var_attr(MeteoGrids::DEM, "z", "altitude", "height above mean sea level", "m", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::TA, "temperature", "air_temperature", "near surface air temperature", "K", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::RH, "humidity", "relative humidity", "relative humidity", "fraction", IOUtils::nodata) );
	tmp.push_back( var_attr(MeteoGrids::P, "pressure", "air_pressure", "near surface air pressure", "Pa", IOUtils::nodata) );
	results["CF1"] = tmp;

	//CNRM schema
	tmp.clear();
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

std::vector<ncParameters::var_attr>  ncParameters::initUserSchemas(const Config& i_cfg)
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
		const size_t param_index = MeteoGrids::getParameterIndex(meteo_grid);
		if (param_index==IOUtils::npos)
			throw InvalidArgumentException("Parameter '"+meteo_grid+"' is not a valid MeteoGrid! Please correct key '"+custom_attr[ii]+"'", AT);
		
		results.push_back( var_attr(param_index, netcdf_param, "", "", "", IOUtils::nodata) );
	}
	
	return results;
}

ncParameters::ncParameters(const std::string& filename, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug)
             : user_schemas( initUserSchemas(cfg) ), vars(), unknown_vars(), vecTime(), vecLat(), vecLon(), dimensions_map(), file_and_path(filename), coordin(), coordinparam(), TZ(tz_in), max_dimension(-1), wrf_hacks(schema=="WRF"), debug(i_debug)
{
	if (!FileUtils::fileExists(file_and_path)) throw AccessException(file_and_path, AT); //prevent invalid filenames
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	if (debug) std::cout << file_and_path << ":\n";
	
	int ncid;
	int status = nc_open(file_and_path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) throw IOException("Could not open netcdf file '" + file_and_path + "': " + nc_strerror(status), AT);
	
	//read the dimensions
	initDimensions(ncid);
	const std::map<ncParameters::Dimensions, nc_dimension>::const_iterator time = dimensions_map.find( TIME );
	if (time!=dimensions_map.end()) vecTime = readTimeDimension(ncid, time->second);
	const std::map<ncParameters::Dimensions, nc_dimension>::const_iterator latitude = dimensions_map.find( LATITUDE );
	if (latitude!=dimensions_map.end()) vecLat = readDimension(ncid, latitude->second);
	const std::map<ncParameters::Dimensions, nc_dimension>::const_iterator longitude = dimensions_map.find( LONGITUDE );
	if (longitude!=dimensions_map.end()) vecLon = readDimension(ncid, longitude->second);
	
	//read all variables
	initVariables(ncid, schema);
	status = nc_close(ncid);
	if (status != NC_NOERR) throw IOException("Could not close netcdf file  '" + file_and_path + "': " + nc_strerror(status), AT);
	
	if (debug) {
		std::cout << "\tParameters: \n";
		for (std::map<size_t, nc_variable>::const_iterator it=vars.begin(); it!=vars.end(); ++it)
			std::cout << "\t\t" << MeteoGrids::getParameterName( it->first ) << " -> " << it->second.toString() << "\n";
		std::cout << "\tUnrecognized variables: \n";
		for (std::map<std::string, nc_variable>::const_iterator it=unknown_vars.begin(); it!=unknown_vars.end(); ++it)
			std::cout << "\t\t" << it->first << " -> " << it->second.toString() << "\n";
		std::cout << "\ttime range: [" << vecTime.front().toString(Date::ISO) << " - " << vecTime.back().toString(Date::ISO) << "]\n";
	}
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
	size_t time_pos = IOUtils::npos;
	if (!date.isUndef()) {
		const std::vector<Date>::const_iterator low = std::lower_bound(vecTime.begin(), vecTime.end(), date);
		if (*low!=date) throw NoDataException("No data at "+date.toString(Date::ISO)+" in file "+file_and_path, AT);
		time_pos = static_cast<size_t>( std::distance(vecTime.begin(), low) );
	} else {
		if (param!=MeteoGrids::DEM && !vecTime.empty()) //HACK check if param depends on TIME dimension
			throw InvalidFormatException("No time requirement has been provided for a file that contains multiple timestamps", AT);
	}
	const std::map <size_t, nc_variable>::const_iterator it = vars.find( param );
	if (it==vars.end()) NoDataException("No "+MeteoGrids::getParameterName( param )+" grid in file "+file_and_path, AT);
	
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
	const double cellsize = calculate_cellsize(resampling_factor_x, resampling_factor_y);
	Grid2DObject grid(vecLon.size(), vecLat.size(), cellsize, llcorner);
	
	//read the raw data, copy it into the Grid2DObject
	int ncid;
	int status = nc_open(file_and_path.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) throw IOException("Could not open netcdf file '" + file_and_path + "': " + nc_strerror(status), AT);
	
	double *data = new double[vecLat.size()*vecLon.size()];
	if (time_pos!=IOUtils::npos)
		ncpp::read_data(ncid, var.attributes.name, var.varid, time_pos, vecLat.size(), vecLon.size(), data);
	else
		ncpp::read_data(ncid, var.attributes.name, var.varid, data);
	fill2DGrid(grid, data, var.nodata);
	delete[] data;
	status = nc_close(ncid);
	if (status != NC_NOERR) throw IOException("Could not close netcdf file  '" + file_and_path + "': " + nc_strerror(status), AT);
	
	//handle data packing and units, if necessary
	if (var.scale!=1.) grid *= var.scale;
	if (var.offset!=0.) grid += var.offset;
	const std::string units( var.units );
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
	
	//HACK expand the definition of Grid2DObject to support lat/lon grids and reproject in GridsManager
	/*if (resampling_factor_x != mio::IOUtils::nodata || resampling_factor_y != mio::IOUtils::nodata) {
		grid.grid2D = mio::LibResampling2D::Bilinear(grid.grid2D, resampling_factor_x, resampling_factor_y);
	}*/	
	return grid;
}

void ncParameters::write2DGrid(Grid2DObject grid_in, const Date& date, const std::string& schema)
{
	const bool is_record = (!date.isUndef());
	
	//define C-like structures for libnetcdf
	double *lat_array = new double[grid_in.getNy()];
	double *lon_array = new double[grid_in.getNx()];
	int *data = new int[grid_in.getNy() * grid_in.getNx()];
	ncpp::calculate_dimensions(grid_in, lat_array, lon_array);
	ncpp::fill_grid_data(grid_in, IOUtils::nodata, data);
	
	bool create_spatial_dimensions(false), create_variable(false), create_time(false);
	int ncid, status;
	if ( FileUtils::fileExists(file_and_path) ) {
		status = nc_open(file_and_path.c_str(), NC_WRITE, &ncid);
		if (status != NC_NOERR) throw IOException("Could not open netcdf file '" + file_and_path + "': " + nc_strerror(status), AT);
		//HACK TODO
		
		ncpp::start_definitions(file_and_path, ncid);
		status = nc_redef(ncid);
		if (status != NC_NOERR) throw IOException("Could not open define mode for file '" + file_and_path + "': " + nc_strerror(status), AT);
	} else {
		if (!FileUtils::validFileAndPath(file_and_path)) throw InvalidNameException(file_and_path, AT);
		ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.6");
		create_variable = create_spatial_dimensions = true;
		if (is_record) create_time = true;
	}
	
	if (create_time) create_dimension(ncid, TIME, schema);
	if (create_spatial_dimensions) {
		create_dimension(ncid, LATITUDE, schema);
		create_dimension(ncid, LONGITUDE, schema);
	}
	
	if (is_record && create_variable) {
		//ncpp::add_3D_variable(ncid, attr.var, NC_INT, did_time, did_lat, did_lon, vid_var); //NC_DOUBLE or NC_INT or NC_SHORT
		//add_attributes_for_variable(ncid, vid_var, attr, IOUtils::nodata);
	} else if (create_variable) {
		/*ncpp::add_2D_variable(ncid, attr.var, NC_INT, did_lat, did_lon, vid_var); //NC_DOUBLE
		add_attributes_for_variable(ncid, vid_var, attr, IOUtils::nodata);
		
		std::vector<int> dimids;
		dimids.push_back(dimensions_map[LATITUDE].dimid);
		dimids.push_back(dimensions_map[LONGITUDE].dimid);

		const int status = nc_def_var(ncid, varname.c_str(), NC_INT, 2, &dimids[0], &varid);
		if (status != NC_NOERR) throw IOException("Could not define variable '" + varname + "': " + nc_strerror(status), AT);*/
		
	}
	status = nc_enddef(ncid);
	if (status != NC_NOERR) throw IOException("Could not close define mode for file '" + file_and_path + "': " + nc_strerror(status), AT);
	
	status = nc_close(ncid);
	if (status != NC_NOERR) throw IOException("Could not close netcdf file  '" + file_and_path + "': " + nc_strerror(status), AT);
	delete[] lat_array; delete[] lon_array; delete[] data;
	
}

void ncParameters::create_dimension(const int& ncid, const Dimensions& dimType, const std::string& schema)
{
	//HACK should we keep that here, or simply fill dimensions_map from schema somewhere else?
	//make sure we have all the metadata for this dimension, either from previous calls (or reads) or from schema
	if (dimensions_map.count(dimType)==0) { //create the dimension from the staic schemas_dims
		if (schemas_dims.count(schema)==0) 
			throw InvalidNameException("The schema '"+schema+"' provided to write gridded data does not exists", AT);
		
		bool dimension_found = false;
		for (size_t ii=0; ii<schemas_dims[schema].size(); ii++) {
			if (schemas_dims[schema][ii].type==dimType) {
				dimensions_map[ dimType ] = nc_dimension( schemas_dims[schema][ii] );
				dimension_found = true;
				break;
			}
		}
		
		if (!dimension_found)
			throw IOException("The requested dimension could not be found in schema '"+schema+"'", AT);
	}
	if (dimType==TIME) {
		dimensions_map[ dimType ].isUnlimited = true;
		//dimensions_map[ dimType ].attributes.units = "hours since 1970-1-1"; //HACK How to handle this in append?
	}
	
	//now create the dimension into the file
	nc_dimension &dim = dimensions_map[ dimType ];
	
	const nc_type length = (dim.isUnlimited)? NC_UNLIMITED : static_cast<int>(dim.length);
	int status = nc_def_dim(ncid, dim.attributes.name.c_str(), length, &dim.dimid);
	if (status != NC_NOERR) throw IOException("Could not define dimension '" + dim.attributes.name + "': " + nc_strerror(status), AT);
	
	//write the dimension's attributes into the file
	static const int ndims = 1;
	status = nc_def_var(ncid, dim.attributes.name.c_str(), NC_DOUBLE, ndims, &dim.dimid, &dim.varid);
	if (status != NC_NOERR) throw IOException("Could not define variable '" + dim.attributes.name + "': " + nc_strerror(status), AT);
	
	status = nc_put_att_text(ncid, dim.varid, "standard_name", dim.attributes.name.size(), dim.attributes.name.c_str());
	if (status != NC_NOERR) throw IOException( std::string("Could not add attribute 'standard_name': ") + nc_strerror(status), AT);
	
	status = nc_put_att_text(ncid, dim.varid, "long_name", dim.attributes.long_name.size(), dim.attributes.long_name.c_str());
	if (status != NC_NOERR) throw IOException( std::string("Could not add attribute 'long_name': ") + nc_strerror(status), AT);
	
	status = nc_put_att_text(ncid, dim.varid, "units", dim.attributes.units.size(), dim.attributes.units.c_str());
	if (status != NC_NOERR) throw IOException( std::string("Could not add attribute 'units': ") + nc_strerror(status), AT);
	
	if (dim.attributes.type==TIME) {
		status = nc_put_att_text(ncid, dim.varid, "calendar", 9, "gregorian");
		if (status != NC_NOERR) throw IOException( std::string("Could not add attribute 'calendar': ") + nc_strerror(status), AT);
	}
}

void ncParameters::initDimensions(const int& ncid)
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
		
		Dimensions type = NONE;
		std::string dimname( name );
		if (dimname=="Time" || dimname=="Times" || dimname=="time")
			type = TIME;
		else if (dimname=="latitude" || dimname=="lat" || dimname=="south_north")
			type = LATITUDE;
		else if (dimname=="longitude" || dimname=="lon" || dimname=="west_east")
			type = LONGITUDE;
		else if (dimname=="x" || dimname=="easting")
			type = EASTING;
		else if (dimname=="y" || dimname=="northing")
			type = NORTHING;
		
		size_t dimlen;
		status = nc_inq_dimlen(ncid, dimids[idx], &dimlen);
		if (status != NC_NOERR) throw IOException("Could not retrieve dimension lenght: " + std::string(nc_strerror(status)), AT);
		
		const std::string units( getAttribute(ncid, dimids[idx], dimname, "units") );
		const std::string long_name( getAttribute(ncid, dimids[idx], dimname, "long_name") );
		if (wrf_hacks) { //wrf does not use the same names between the dimensions list and the variables dependencies...
			if (type==TIME) dimname = "Times";
			if (type==LATITUDE) dimname = "XLAT";
			if (type==LONGITUDE) dimname = "XLONG";
		}
		
		dimensions_map[ type ] = nc_dimension( dim_attributes(type, dimname, long_name, units), dimlen, idx, -1, (idx==unlim_id));
		if (idx>max_dimension) max_dimension = idx;
	}
	
	free( dimids );
	
	if (debug) {
		std::map<ncParameters::Dimensions, nc_dimension>::const_iterator it;
		for(it=dimensions_map.begin(); it!=dimensions_map.end(); ++it)
			std::cout << "\t" << it->second.toString() << "\n";
	}
}

void ncParameters::initVariables(const int& ncid, const std::string& schema_name)
{
	int nr_of_variables = -1;
	int status = nc_inq_nvars(ncid, &nr_of_variables); //Trick! this also returns the dimensions!
	if (status != NC_NOERR)
		throw IOException("Could not retrieve variables for dataset: " + string(nc_strerror(status)), AT);

	// Variable IDs in a NetCDF file are consecutive integers starting with 0
	for (int ii=0; ii<nr_of_variables; ++ii) {
		char name[NC_MAX_NAME+1];
		status = nc_inq_varname(ncid, ii, name);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
		const std::string varname( name );
		
		int varid = -1;
		status = nc_inq_varid (ncid, name, &varid);
		if (status != NC_NOERR) throw IOException(nc_strerror(status), AT);
		if (varid<=max_dimension) continue;
		
		nc_variable tmp_var(getSchemaAttributes(varname, schema_name), "", 1., 0., IOUtils::nodata, varid);
		getAttribute(ncid, tmp_var.varid, varname, "missing_value", tmp_var.nodata);
		getAttribute(ncid, tmp_var.varid, varname, "scale_factor", tmp_var.scale);
		getAttribute(ncid, tmp_var.varid, varname, "add_offset", tmp_var.offset);
		tmp_var.units = getAttribute(ncid, tmp_var.varid, varname, "units");
		
		if (tmp_var.attributes.param!=IOUtils::npos)
			vars[ tmp_var.attributes.param ] = tmp_var;
		else
			unknown_vars[ varname ] = tmp_var;
	}
}

std::vector<Date> ncParameters::readTimeDimension(const int& ncid, const nc_dimension& dim) const
{
	std::vector<Date> results(dim.length);
	
	if (!wrf_hacks) {
		//get time decoding parameters
		double offset = 0., multiplier = 1.;
		if (!wrf_hacks) getTimeTransform(dim.attributes.units, TZ, offset, multiplier);
		
		double *data = new double[ dim.length ];
		ncpp::read_data(ncid, dim.attributes.name, dim.dimid, data);
		for (size_t ii=0; ii<dim.length; ii++)
			results[ii].setDate(data[ii]*multiplier + offset, TZ);
		delete[] data;
	} else {
		static const size_t DateStrLen = 19; //HACK DateStrLen = 19, defined in Dimensions

		char *record_value = (char*)calloc(dim.length, sizeof(char)*DateStrLen);
		const int status = nc_get_var_text(ncid, dim.dimid, record_value);
		if (status != NC_NOERR)
			throw IOException("Could not retrieve data for Time variable: " + std::string(nc_strerror(status)), AT);

		for(size_t ii=0; ii<dim.length; ii++) {
			std::string tmp(DateStrLen, '\0');
			for(size_t jj=0; jj<DateStrLen; jj++) {
				const char c = record_value[ii*DateStrLen+jj];
				tmp[jj] = (c!='_')? c : 'T';
			}
			IOUtils::convertString(results[ii], tmp, TZ);
		}
		free( record_value );
	}
	
	return results;
}

std::vector<double> ncParameters::readDimension(const int& ncid, const nc_dimension& dim)
{
	std::vector<double> results( dim.length );
	double *data = new double[ dim.length ];
	ncpp::read_data(ncid, dim.attributes.name, dim.dimid, data);
	std::copy(data, data+dim.length, results.begin());
	delete[] data;
	return results;
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

//if the attribute is not found, an empty string is returned
std::string ncParameters::getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name)
{
	size_t attr_len;
	int status = nc_inq_attlen (ncid, value_id, attr_name.c_str(), &attr_len);
	if (status == NC_NOERR) {
		char* value = new char[attr_len + 1]; // +1 for trailing null
		status = nc_get_att_text(ncid, value_id, attr_name.c_str(), value);
		if (status != NC_NOERR) throw IOException("Could not read attribute '" + attr_name + "' for '" + value_name + "': " + nc_strerror(status), AT);

		value[attr_len] = '\0';
		const std::string attr_value(value);
		delete[] value;
		return attr_value;
	}
	return std::string();
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

} //namespace
