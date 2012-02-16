/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include "GRIBIO.h"

#include <meteoio/meteolaws/Atmosphere.h>
#include <meteoio/meteolaws/Meteoconst.h> //for PI
#include <meteoio/DEMObject.h>

#include <errno.h>
#include <grib_api.h>

using namespace std;

namespace mio {
/**
 * @page gribio GRIBIO
 * @section gribio_format Format
 * This plugin reads GRIB (https://en.wikipedia.org/wiki/GRIB) files as produced by meteorological models.
 * Being based on GRIB API (http://www.ecmwf.int/products/data/software/grib_api.html), it should support both version 1 and 2 of the format.
 * This has been developed for reading MeteoSwiss Cosmo data and reads fields based on their marsParam code (see for example
 * http://www-imk.fzk.de/~kouker/mars/param.html even if these don't match excatly MeteoSwiss...)
 *
 * Levels description is available at also http://www.nco.ncep.noaa.gov/pmb/docs/on388/ . no correction for grid rotation is currently performed.
 * read2DGrid(grid, filename) returns the first grid that is found.
 *
 * @section gribio_units Units
 * As specified by WMO.
 *
 * @section gribio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords)
 * - COORDPARAM: extra coordinates parameters (see Coords)
 * - GRID2DPATH: path where to find the grids
 * - GRIB_DEM_UPDATE: recompute slope/azimuth from the elevations when reading a DEM (default=false,
 * that is we use the slope and azimuth included in the GRIB file)
 * - STATION#: coordinates for virtual stations (if using GRIB as METEO plugin). Each station is given by its coordinates and the closest
 * grid point will be chosen. Coordinates are given one one line as "lat lon" or "xcoord ycoord epsg_code"
 *
 */

const double GRIBIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const double GRIBIO::tz_in = 0.; //GRIB time zone, always UTC
const std::string GRIBIO::prefix="laf"; //filename prefix
const std::string GRIBIO::ext=".grb"; //filename extension

GRIBIO::GRIBIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	setOptions();
	indexed = false;
	idx=NULL;
	fp = NULL;
	meteopath_scanned = false;
}

GRIBIO::GRIBIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	setOptions();
	indexed = false;
	idx=NULL;
	fp = NULL;
	meteopath_scanned = false;
}

GRIBIO::GRIBIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	setOptions();
	indexed = false;
	idx=NULL;
	fp = NULL;
	meteopath_scanned = false;
}

GRIBIO::~GRIBIO() throw()
{
	cleanup();
}

void GRIBIO::setOptions()
{
	std::string coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	update_dem = false;

	string tmp="";
	cfg.getValue("GRID2D", "Input", tmp, Config::nothrow);
	if (tmp == "GRIB") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
		cfg.getValue("GRIB_DEM_UPDATE", "Input", update_dem, Config::nothrow);
	}

	tmp="";
	cfg.getValue("METEO", "Input", tmp, Config::nothrow);
	if (tmp == "GRIB") { //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("METEOPATH", "Input", meteopath_in);
		size_t current_stationnr = 1;
		string current_station;
		do {
			current_station = "";
			stringstream ss;
			ss << "STATION" << current_stationnr;
			cfg.getValue(ss.str(), "Input", current_station, Config::nothrow);

			if (current_station != "") {
				addStation(current_station);
				std::cout <<  "\tRead virtual station " << vecPts.back().printLatLon() << "\n";
			}
			current_stationnr++;
		} while (current_station != "");
	}
}

void GRIBIO::addStation(const std::string& coord_spec)
{
	std::istringstream iss(coord_spec);
	double coord1=IOUtils::nodata, coord2=IOUtils::nodata;
	int epsg=IOUtils::inodata;

	iss >> std::skipws >> coord1;
	iss >> std::skipws >> coord2;
	iss >> std::skipws >> epsg;

	if(coord1!=IOUtils::nodata && coord2!=IOUtils::nodata && epsg!=IOUtils::inodata) {
		Coords point;
		point.setEPSG(epsg);
		point.setXY(coord1, coord2, IOUtils::nodata);
		vecPts.push_back(point);
		return;
	}
	if(coord1!=IOUtils::nodata && coord2!=IOUtils::nodata) {
		Coords point(coordin, coordinparam);
		point.setLatLon(coord1, coord2, IOUtils::nodata);
		vecPts.push_back(point);
		return;
	}

	throw InvalidArgumentException("Coordinate specification \""+coord_spec+"\" is invalid!", AT);
}

void GRIBIO::listKeys(grib_handle** h, const std::string& filename)
{
	const unsigned int MAX_VAL_LEN=1024; //max value string length in GRIB
	unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;
	char* name_space=NULL;
	grib_keys_iterator* kiter=grib_keys_iterator_new(*h,key_iterator_filter_flags,name_space);

	if (!kiter) {
		cleanup();
		throw IOException("Unable to create keys iterator for \""+filename+"\"", AT);
	}

	//Iterating over all keys
	while(grib_keys_iterator_next(kiter)) {
		char value[MAX_VAL_LEN];
		size_t vlen=MAX_VAL_LEN;
		const char* name = grib_keys_iterator_get_name(kiter);
		vlen=MAX_VAL_LEN;
		bzero(value,vlen);
		GRIB_CHECK(grib_get_string(*h,name,value,&vlen),name);
		std::cout << name << " = " << value << "\n";
	}

	grib_keys_iterator_delete(kiter);
}

void  GRIBIO::listFields(const std::string& filename)
{
	grib_handle* h=NULL;
	int err=0;

	//loop over the messages
	while((h = grib_handle_new_from_file(0,fp,&err)) != NULL) {
		if(!h) {
			cleanup();
			throw IOException("Unable to create grib handle for \""+filename+"\"", AT);
		}

		double marsParam;
		GRIB_CHECK(grib_get_double(h,"marsParam",&marsParam),0);

		size_t len=500;
		char name[500], shortname[500];
		GRIB_CHECK(grib_get_string(h,"name",name, &len),0);
		len=500; //this has been reset by the previous call... no comments
		GRIB_CHECK(grib_get_string(h,"shortName",shortname, &len),0);
		long levelType;
		GRIB_CHECK(grib_get_long(h,"indicatorOfTypeOfLevel", &levelType),0);
		long level=0;
		if(levelType!=1) GRIB_CHECK(grib_get_long(h,"level", &level),0);
		std::cout << marsParam << " " << shortname << " " << name << " type " << levelType << " level " << level << "\n";
		grib_handle_delete(h);
	}
}

Date GRIBIO::getDate(grib_handle* h) {
	long dataDate, dataTime;
	GRIB_CHECK(grib_get_long(h,"dataDate",&dataDate),0);
	GRIB_CHECK(grib_get_long(h,"dataTime",&dataTime),0);

	const int year=dataDate/10000, month=dataDate/100-year*100, day=dataDate-month*100-year*10000;
	const int hour=dataTime/100, minutes=dataTime-hour*100;
	return Date(year, month, day, hour, minutes, tz_in);
}

Coords GRIBIO::getGeolocalization(grib_handle* h, double &cellsize_x, double &cellsize_y)
{
	//getting transformation parameters
	double angleOfRotationInDegrees;
	GRIB_CHECK(grib_get_double(h,"angleOfRotationInDegrees",&angleOfRotationInDegrees),0);
	if(angleOfRotationInDegrees!=0.) {
		throw InvalidArgumentException("Rotated grids not supported!", AT);
	}
	double latitudeOfSouthernPole, longitudeOfSouthernPole;
	GRIB_CHECK(grib_get_double(h,"latitudeOfSouthernPoleInDegrees",&latitudeOfSouthernPole),0);
	GRIB_CHECK(grib_get_double(h,"longitudeOfSouthernPoleInDegrees",&longitudeOfSouthernPole),0);

	//determining llcorner
	double ll_latitude, ll_longitude;
	GRIB_CHECK(grib_get_double(h,"latitudeOfFirstGridPointInDegrees",&ll_latitude),0);
	GRIB_CHECK(grib_get_double(h,"longitudeOfFirstGridPointInDegrees",&ll_longitude),0);
	Coords llcorner(coordin, coordinparam);
	llcorner.setLatLon( ll_latitude+(90.+latitudeOfSouthernPole), ll_longitude+longitudeOfSouthernPole, 0.);

	//determining cell size
	double ur_latitude, ur_longitude;
	GRIB_CHECK(grib_get_double(h,"latitudeOfLastGridPointInDegrees",&ur_latitude),0);
	GRIB_CHECK(grib_get_double(h,"longitudeOfLastGridPointInDegrees",&ur_longitude),0);
	const double cntr_latitude = .5*(ll_latitude+ur_latitude)+(90.+latitudeOfSouthernPole);

	double d_i, d_j;
	GRIB_CHECK(grib_get_double(h,"jDirectionIncrementInDegrees",&d_j),0);
	GRIB_CHECK(grib_get_double(h,"iDirectionIncrementInDegrees",&d_i),0);

	cellsize_x = Coords::lon_degree_lenght(cntr_latitude)*d_i;
	cellsize_y = Coords::lat_degree_lenght(cntr_latitude)*d_j;

	//checking that cellsize does not vary too much across the grid
	const double cellsize_x_ll = Coords::lon_degree_lenght(ll_latitude)*d_j;
	const double cellsize_x_ur = Coords::lon_degree_lenght(ur_latitude)*d_j;
	if( fabs(cellsize_x_ll-cellsize_x_ur)/cellsize_x > 1./100.) {
		stringstream ss;
		ss << "Cell size varying too much in the x direction between lower left and upper right corner: ";
		ss << cellsize_x_ll << "m to " << cellsize_x_ur << "m";
		throw IOException(ss.str(), AT);
	}

	return llcorner;
}

void GRIBIO::read2Dlevel(grib_handle* h, Grid2DObject& grid_out)
{ //HACK: why is it that we don't have to corect the aspect ratio??
	long Ni, Nj;
	GRIB_CHECK(grib_get_long(h,"Ni",&Ni),0);
	GRIB_CHECK(grib_get_long(h,"Nj",&Nj),0);

	double *values;
	size_t values_len= 0;
	GRIB_CHECK(grib_get_size(h,"values",&values_len),0);
	if(values_len!=(unsigned)(Ni*Nj)) {
		stringstream ss;
		ss << "Declaring grid of size " << Ni << "x" << Nj << "=" << Ni*Nj << " ";
		ss << "but containing " << values_len << " values. This is inconsistent!";
		throw InvalidArgumentException(ss.str(), AT);
	}
	values = (double*)malloc(values_len*sizeof(double));

	GRIB_CHECK(grib_get_double_array(h,"values",values,&values_len),0);
	double cellsize_x, cellsize_y;
	const Coords llcorner = getGeolocalization(h, cellsize_x, cellsize_y);
	grid_out.set(static_cast<unsigned int>(Ni), static_cast<unsigned int>(Nj), cellsize_x, llcorner); //HACK
	int i=0;
	for(unsigned int jj=0; jj<(unsigned)Nj; jj++) {
		for(unsigned int ii=0; ii<(unsigned)Ni; ii++)
			grid_out(ii,jj) = values[i++];
	}
	free(values);
}

bool GRIBIO::read2DGrid_indexed(const double& in_marsParam, const long& i_levelType, const long& i_level, const Date i_date, Grid2DObject& grid_out)
{
	GRIB_CHECK(grib_index_select_double(idx,"marsParam",in_marsParam),0);
	GRIB_CHECK(grib_index_select_long(idx,"indicatorOfTypeOfLevel", i_levelType),0);

	grib_handle* h=NULL;
	int err=0;
	while((h = grib_handle_new_from_index(idx,&err)) != NULL) {
		if(!h) {
			cleanup();
			throw IOException("Unable to create grib handle from index", AT);
		}

		double timeRange=0.;
		GRIB_CHECK(grib_get_double(h,"timeRangeIndicator", &timeRange),0);
		const Date validity_start=getDate(h);
		const Date validity_end=validity_start+timeRange/24.;
		long level=0;
		if(i_level!=0) GRIB_CHECK(grib_get_long(h,"level", &level),0);
		if(level==i_level && (i_date.isUndef() || (i_date>=validity_start && i_date<=validity_end)) ) {
			read2Dlevel(h, grid_out);
			return true;
		}
		grib_handle_delete(h);
	}
	return false;
}

void GRIBIO::read2DGrid(Grid2DObject& grid_out, const std::string& i_name)
{
	const std::string filename = grid2dpath_in+"/"+i_name;
	fp = fopen(filename.c_str(),"r");
	if(fp==NULL) {
		stringstream ss;
		ss << "Error openning file \"" << filename << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	grib_handle* h=NULL;
	int err=0;
	if((h = grib_handle_new_from_file(0,fp,&err)) != NULL) {
		if(!h) {
			cleanup();
			throw IOException("Unable to create grib handle for \""+filename+"\"", AT);
		}

		read2Dlevel(h, grid_out);
		grib_handle_delete(h);
	} else {
		cleanup();
		throw IOException("No grid found in file \""+filename+"\"", AT);
	}

	cleanup();
}

void GRIBIO::indexFile(const std::string& filename)
{
	fp = fopen(filename.c_str(),"r");
	if(fp==NULL) {
		stringstream ss;
		ss << "Error openning file \"" << filename << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	int err=0;
	std::string keys("marsParam:d,indicatorOfTypeOfLevel:l"); //indexing keys
	char *c_filename = (char *)filename.c_str();
	idx = grib_index_new_from_file(0, c_filename, keys.c_str(), &err);
	if(err!=0) {
		cleanup();
		throw IOException("Failed to index GRIB file \""+filename+"\"", AT);
	}
	indexed=true;
	idx_filename=filename;
}

void GRIBIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	Date UTC_date = date;
	UTC_date.setTimeZone(tz_in);

	const std::string filename = grid2dpath_in+"/"+prefix+UTC_date.toString(Date::NUM).substr(0,10)+"f"+UTC_date.toString(Date::NUM).substr(10,2)+ext;

	read2DGrid(filename, grid_out, parameter, UTC_date);
}

void GRIBIO::read2DGrid(const std::string& filename, Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if(!indexed || idx_filename!=filename) {
		cleanup();
		indexFile(filename);
	}

	//Basic meteo parameters
	if(parameter==MeteoGrids::P) read2DGrid_indexed(1.2, 1, 0, date, grid_out); //PS what is HMO3?
	if(parameter==MeteoGrids::TA) read2DGrid_indexed(11.2, 105, 2, date, grid_out); //T_2M
	if(parameter==MeteoGrids::RH) {
		if(!read2DGrid_indexed(52.2, 105, 2, date, grid_out)) { //RELHUM_2M
			Grid2DObject ta;
			read2DGrid_indexed(11.2, 105, 2, date, ta); //T_2M
			read2DGrid_indexed(17.2, 105, 2, date, grid_out); //TD_2M
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = Atmosphere::DewPointtoRh(grid_out(ii,jj), ta(ii,jj), true);
				}
			}
		}
	}
	if(parameter==MeteoGrids::TSS) read2DGrid_indexed(197.201, 111, 0, date, grid_out); //T_SO
	if(parameter==MeteoGrids::TSG) read2DGrid_indexed(11.2, 1, 0, date, grid_out); //T_G

	//hydrological parameters
	if(parameter==MeteoGrids::HNW) read2DGrid_indexed(61.2, 1, 0, date, grid_out); //tp
	if(parameter==MeteoGrids::ROT) read2DGrid_indexed(90.2, 112, 0, date, grid_out); //RUNOFF
	if(parameter==MeteoGrids::SWE) read2DGrid_indexed(65.2, 1, 0, date, grid_out); //W_SNOW
	if(parameter==MeteoGrids::HS) {
		Grid2DObject snow_density;
		read2DGrid_indexed(133.201, 1, 0, date, snow_density); //RHO_SNOW
		read2DGrid_indexed(65.2, 1, 0, date, grid_out); //W_SNOW
		grid_out.grid2D /= snow_density.grid2D;
	}

	//radiation parameters
	if(parameter==MeteoGrids::ALB) {
		read2DGrid_indexed(84.2, 1, 0, date, grid_out); //ALB_RAD
		grid_out.grid2D /= 100.;
	}
	if(parameter==MeteoGrids::ILWR) read2DGrid_indexed(25.201, 1, 0, date, grid_out); //ALWD_S
	if(parameter==MeteoGrids::ISWR) {
		if(!read2DGrid_indexed(111.250, 1, 0, date, grid_out)) { //GLOB
			Grid2DObject diff;
			read2DGrid_indexed(23.201, 1, 0, date, diff); //diffuse rad, ASWDIFD_S
			read2DGrid_indexed(22.201, 1, 0, date, grid_out); //direct rad, ASWDIR_S
			grid_out.grid2D += diff.grid2D;
		}
	}

	//DEM parameters
	if(parameter==MeteoGrids::DEM) read2DGrid_indexed(8.2, 1, 0, date, grid_out); //HSURF
	if(parameter==MeteoGrids::SLOPE) read2DGrid_indexed(98.202, 1, 0, date, grid_out); //SLO_ANG
	if(parameter==MeteoGrids::AZI) read2DGrid_indexed(99.202, 1, 0, date, grid_out); //SLO_ASP

	//Wind parameters
	if(parameter==MeteoGrids::U) read2DGrid_indexed(33.2, 105, 10, date, grid_out); //U_10M, also in 110, 10 as U
	if(parameter==MeteoGrids::V) read2DGrid_indexed(34.2, 105, 10, date, grid_out); //V_10M, also in 110, 10 as V
	if(parameter==MeteoGrids::W) read2DGrid_indexed(40.2, 109, 10, date, grid_out); //W, 10m
	if(parameter==MeteoGrids::VW_MAX) read2DGrid_indexed(187.201, 105, 10, date, grid_out); //VMAX_10M 10m
	if(parameter==MeteoGrids::DW) {
		if(!read2DGrid_indexed(31.2, 105, 10, date, grid_out)) { //DD_10M
			Grid2DObject V;
			read2DGrid_indexed(34.2, 105, 10, date, V); //V_10M
			read2DGrid_indexed(33.2, 105, 10, date, grid_out); //U_10M
			const double to_deg = 180. / Cst::PI;
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = fmod( atan2( grid_out(ii,jj), V(ii,jj) ) * to_deg + 360., 360.); // turn into degrees [0;360)
				}
			}
		}
	}
	if(parameter==MeteoGrids::VW) {
		if(!read2DGrid_indexed(32.2, 105, 10, date, grid_out)) { //FF_10M
			Grid2DObject V;
			read2DGrid_indexed(34.2, 105, 10, date, V); //V_10M
			read2DGrid_indexed(33.2, 105, 10, date, grid_out); //U_10M
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = sqrt( IOUtils::pow2(grid_out(ii,jj)) + IOUtils::pow2(V(ii,jj)) );
				}
			}
		}
	}

	if(grid_out.isEmpty()) {
		stringstream ss;
		ss << "No suitable data found for parameter " << MeteoGrids::getParameterName(parameter) << " ";
		ss << "at time step " << date.toString(Date::ISO) << " in file \"" << filename << "\"";
		throw NoAvailableDataException(ss.str(), AT);
	}

	//correcting wind speeds
	/*if(parameter==MeteoGrids::U || parameter==MeteoGrids::V || parameter==MeteoGrids::W || parameter==MeteoGrids::VW || parameter==MeteoGrids::VW_MAX) {
		//we need to compute the wind at 7.5m
		Grid2DObject Z0;
		if(read2DGrid_indexed(83.2, 1, 0, date, Z0)) { //Z0
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = Atmosphere::windLogProfile(grid_out(ii,jj), 10., 7.5, Z0(ii,jj));
				}
			}
		} else {
			const double wind_factor = Atmosphere::windLogProfile(1., 10., 7.5, 0.03);
			grid_out.grid2D *= wind_factor;
		}
	}*/
}

void GRIBIO::readDEM(DEMObject& dem_out)
{
	const Date d; //ie: undef. This will be caught when reading the GRIB file
	std::string filename;

	cfg.getValue("DEMFILE", "Input", filename);
	read2DGrid(filename, dem_out, MeteoGrids::DEM, d);
	if(update_dem) {
		dem_out.update();
	} else {
		const int dem_ppt=dem_out.getUpdatePpt();
		if(dem_ppt&DEMObject::SLOPE) {
			Grid2DObject slope;
			read2DGrid(filename, slope, MeteoGrids::SLOPE, d);
			dem_out.slope=slope.grid2D;
			Grid2DObject azi;
			read2DGrid(filename, azi, MeteoGrids::AZI, d);
			dem_out.azi=azi.grid2D;
		}
		if(dem_ppt&DEMObject::NORMAL || dem_ppt&DEMObject::CURVATURE) {
			//we will only update the normals and/or curvatures, then revert update properties
			if(dem_ppt&DEMObject::NORMAL && dem_ppt&DEMObject::CURVATURE) dem_out.setUpdatePpt((DEMObject::update_type)(DEMObject::NORMAL|DEMObject::CURVATURE));
			else if(dem_ppt&DEMObject::NORMAL) dem_out.setUpdatePpt(DEMObject::NORMAL);
			else if(dem_ppt&DEMObject::CURVATURE) dem_out.setUpdatePpt(DEMObject::CURVATURE);

			dem_out.update();
			dem_out.setUpdatePpt((DEMObject::update_type)dem_ppt);
		}

		dem_out.updateAllMinMax();
	}
}

void GRIBIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::scanMeteoPath()
{
	std::list<std::string> dirlist;
	IOUtils::readDirectory(meteopath_in, dirlist, ext);
	dirlist.sort();

	//Check date in every filename and cache it
	std::list<std::string>::iterator it = dirlist.begin();
	while ((it != dirlist.end())) {
		const std::string& filename = *it;
		std::string::size_type spos = filename.find_first_of("0123456789");
		Date date;
		IOUtils::convertString(date, filename.substr(spos,10)+filename.substr(spos+11,2), tz_in); //HACK MeteoSwiss-only file naming?
		std::pair<Date,std::string> tmp(date, filename);

		cache_meteo_files.push_back(tmp);
		it++;
	}
	meteopath_scanned=true;
}

void GRIBIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo,
                             const size_t&)
{
	if(!meteopath_scanned) scanMeteoPath();

	vecMeteo.clear();
	vecMeteo.insert(vecMeteo.begin(), vecPts.size(), std::vector<MeteoData>()); //allocation for the vectors

	double *lats = (double*)malloc(vecPts.size()*sizeof(double));
	double *lons = (double*)malloc(vecPts.size()*sizeof(double));
	double latitudeOfSouthernPole, longitudeOfSouthernPole;
	std::vector<StationData> meta; //metadata for meteo time series
	bool meta_ok=false; //set to true once the metadata have been read

	//find index of first time step
	unsigned int idx_start;
	bool start_found=false;
	for(idx_start=0; idx_start<cache_meteo_files.size(); idx_start++) {
		if(dateStart<cache_meteo_files[idx_start].first) {
			start_found=true;
			break;
		}
	}
	if(start_found==false) return;
	if(idx_start>0) idx_start--; //start with first element before dateStart (useful for resampling)

	try {
		for(unsigned int ii=idx_start; ii<cache_meteo_files.size(); ii++) {
			const Date& date = cache_meteo_files[ii].first;
			if(date>dateEnd) break;
			const std::string filename = meteopath_in+"/"+cache_meteo_files[ii].second;

			if(!indexed || idx_filename!=filename) {
				cleanup();
				indexFile(filename);
			}
			if(!meta_ok) {
				readMeteoMeta(vecPts, meta, latitudeOfSouthernPole, longitudeOfSouthernPole, lats, lons);
				//HACK: remove potential duplicates in the returned points.
				//This can be done by checking on station ID, since it contains the grid point index
				meta_ok=true;
			}

			std::vector<MeteoData> Meteo;
			readMeteoStep(meta, lats, lons, date, Meteo);
			for(unsigned int jj=0; jj<vecPts.size(); jj++)
				vecMeteo[jj].push_back(Meteo[jj]);
		}
	} catch(...) {
		free(lats); free(lons);
		cleanup();
		throw;
	}

	free(lats); free(lons);
}

void GRIBIO::readMeteoMeta(const std::vector<Coords>& vecPts, std::vector<StationData> &stations, double &latitudeOfSouthernPole, double &longitudeOfSouthernPole, double *lats, double *lons)
{
	GRIB_CHECK(grib_index_select_double(idx,"marsParam",8.2),0); //This is the DEM
	GRIB_CHECK(grib_index_select_long(idx,"indicatorOfTypeOfLevel", 1),0);

	grib_handle* h=NULL;
	int err=0;
	if((h = grib_handle_new_from_index(idx,&err)) != NULL) {
		if(!h) {
			cleanup();
			throw IOException("Unable to create grib handle from index", AT);
		}
		const long npoints = vecPts.size();

		GRIB_CHECK(grib_get_double(h,"latitudeOfSouthernPoleInDegrees",&latitudeOfSouthernPole),0);
		GRIB_CHECK(grib_get_double(h,"longitudeOfSouthernPoleInDegrees",&longitudeOfSouthernPole),0);

		long Ni;
		GRIB_CHECK(grib_get_long(h,"Ni",&Ni),0);

		//build GRIB local coordinates for the points
		for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
			lats[ii] = vecPts[ii].getLat() - (90.+latitudeOfSouthernPole);
			lons[ii] = vecPts[ii].getLon() - longitudeOfSouthernPole;
		}

		//retrieve nearest points
		double *outlats = (double*)malloc(npoints*sizeof(double));
		double *outlons = (double*)malloc(npoints*sizeof(double));
		double *values = (double*)malloc(npoints*sizeof(double));
		double *distances = (double*)malloc(npoints*sizeof(double));
		int *indexes = (int *)malloc(npoints*sizeof(int));
		if(grib_nearest_find_multiple(h, 0, lats, lons, npoints, outlats, outlons, values, distances, indexes)!=0) {
			grib_handle_delete(h);
			cleanup();
			throw IOException("Errro when searching for nearest points in DEM", AT);
		}

		//fill metadata
		for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
			StationData sd;
			sd.position.setProj(coordin, coordinparam);
			sd.position.setLatLon(outlats[ii]+(90.+latitudeOfSouthernPole), outlons[ii]+longitudeOfSouthernPole, values[ii]);
			stringstream ss;
			ss << "Point_" << indexes[ii];
			sd.stationID=ss.str();
			stringstream ss2;
			ss2 << "GRIB point (" << indexes[ii] % Ni << "," << indexes[ii] / Ni << ")";
			sd.stationName=ss2.str();
			stations.push_back(sd);
		}

		free(outlats); free(outlons); free(values); free(distances); free(indexes);
		grib_handle_delete(h);
	} else {
		throw NoAvailableDataException("Can not find DEM grid in GRIB file!", AT);
	}
}

bool GRIBIO::readMeteoValues(const double& marsParam, const long& levelType, const long& i_level, const Date& i_date, const long& npoints, double *lats, double *lons, double *values)
{
	GRIB_CHECK(grib_index_select_double(idx,"marsParam",marsParam),0);
	GRIB_CHECK(grib_index_select_long(idx,"indicatorOfTypeOfLevel", levelType),0);

	grib_handle* h=NULL;
	int err=0;
	while((h = grib_handle_new_from_index(idx,&err)) != NULL) {
		if(!h) {
			cleanup();
			throw IOException("Unable to create grib handle from index", AT);
		}

		double timeRange=0.;
		GRIB_CHECK(grib_get_double(h,"timeRangeIndicator", &timeRange),0);
		const Date validity_start=getDate(h);
		const Date validity_end=validity_start+timeRange/24.;
		long level=0;
		if(i_level!=0) GRIB_CHECK(grib_get_long(h,"level", &level),0);
		if(level==i_level && (i_date.isUndef() || (i_date>=validity_start && i_date<=validity_end)) ) {
			double *outlats = (double*)malloc(npoints*sizeof(double));
			double *outlons = (double*)malloc(npoints*sizeof(double));
			double *distances = (double*)malloc(npoints*sizeof(double));
			int *indexes = (int *)malloc(npoints*sizeof(int));
			if(grib_nearest_find_multiple(h, 0, lats, lons, npoints, outlats, outlons, values, distances, indexes)!=0) {
				grib_handle_delete(h);
				cleanup();
				throw IOException("Errro when searching for nearest points in DEM", AT);
			}

			free(outlats); free(outlons); free(distances); free(indexes);
			grib_handle_delete(h);
			return true;
		}
		grib_handle_delete(h);
	}
	return false;
}

void GRIBIO::fillMeteo(double *values, const MeteoData::Parameters& param, const long& npoints, std::vector<MeteoData> &Meteo) {
	for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
		Meteo[ii](param) = values[ii];
	}
}

void GRIBIO::readMeteoStep(std::vector<StationData> &stations, double *lats, double *lons, const Date i_date, std::vector<MeteoData> &Meteo)
{
	const long npoints = stations.size();

	for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
		MeteoData md;
		md.meta = stations[ii];
		md.date = i_date;
		Meteo.push_back(md);
	}

	double *values = (double*)malloc(npoints*sizeof(double));
	double *values2 = (double*)malloc(npoints*sizeof(double)); //for extra parameters

	//basic meteorological parameters
	if(readMeteoValues(1.2, 1, 0, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::P, npoints, Meteo); //PS
	if(readMeteoValues(11.2, 105, 2, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::TA, npoints, Meteo); //T_2M
	if(readMeteoValues(197.201, 111, 0, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::TSS, npoints, Meteo); //T_SO
	if(readMeteoValues(11.2, 1, 0, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::TSG, npoints, Meteo); //T_G
	if(readMeteoValues(52.2, 105, 2, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::RH, npoints, Meteo); //RELHUM_2M
	else if(readMeteoValues(17.2, 105, 2, i_date, npoints, lats, lons, values)) { //TD_2M
		for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
			if(Meteo[ii](MeteoData::TA)!=IOUtils::nodata)
				Meteo[ii](MeteoData::RH) = Atmosphere::DewPointtoRh(values[ii], Meteo[ii](MeteoData::TA), true);
		}
	}

	//hydrological parameters
	if(readMeteoValues(61.2, 1, 0, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::HNW, npoints, Meteo); //tp
	if(readMeteoValues(133.201, 1, 0, i_date, npoints, lats, lons, values)  //RHO_SNOW
	   && readMeteoValues(65.2, 1, 0, i_date, npoints, lats, lons, values2)) { //W_SNOW
		for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
			Meteo[ii](MeteoData::HS) = values2[ii] / values[ii];
		}
	}

	//radiation parameters
	if(readMeteoValues(25.201, 1, 0, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::ILWR, npoints, Meteo); //ALWD_S
	if(readMeteoValues(111.250, 1, 0, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::ISWR, npoints, Meteo); //GLOB
	else {
		if(readMeteoValues(23.201, 1, 0, i_date, npoints, lats, lons, values) //ASWDIFD_S
		   && readMeteoValues(22.201, 1, 0, i_date, npoints, lats, lons, values2)) { //ASWDIR_S
			for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
				Meteo[ii](MeteoData::ISWR) = values[ii] + values2[ii];
			}
		}
	}
	if(readMeteoValues(84.2, 1, 0, i_date, npoints, lats, lons, values)) { //ALB_RAD
		for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
			if(Meteo[ii](MeteoData::ISWR)!=IOUtils::nodata) Meteo[ii](MeteoData::RSWR) = Meteo[ii](MeteoData::ISWR) * values[ii]/100.;
		}
	}

	//Wind parameters
	if(readMeteoValues(187.201, 105, 10, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::VW_MAX, npoints, Meteo); //VMAX_10M
	if(readMeteoValues(31.2, 105, 10, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::DW, npoints, Meteo); //DD_10M
	else {
		if(readMeteoValues(34.2, 105, 10, i_date, npoints, lats, lons, values) //V_10M
		   && readMeteoValues(33.2, 105, 10, i_date, npoints, lats, lons, values2)) { //U_10M
			const double to_deg = 180. / Cst::PI;
			for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
				Meteo[ii](MeteoData::DW) = fmod( atan2( values2[ii], values[ii] ) * to_deg + 360., 360.); // turn into degrees [0;360)
			}
		}
	}
	if(readMeteoValues(32.2, 105, 10, i_date, npoints, lats, lons, values)) fillMeteo(values, MeteoData::VW, npoints, Meteo); //FF_10M
	else {
		if(readMeteoValues(34.2, 105, 10, i_date, npoints, lats, lons, values) //V_10M
		   && readMeteoValues(33.2, 105, 10, i_date, npoints, lats, lons, values2)) { //U_10M
			for(unsigned int ii=0; ii<(unsigned)npoints; ii++) {
				Meteo[ii](MeteoData::VW) =  sqrt( IOUtils::pow2(values[ii]) + IOUtils::pow2(values2[ii]) );
			}
		}
	}

	free(values); free(values2);
}

void GRIBIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                              const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::write2DGrid(const Grid2DObject& /*grid_in*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::cleanup() throw()
{
	if(fp!=NULL) fclose(fp); fp=NULL;
	if(idx!=NULL) grib_index_delete(idx); idx=NULL;
	idx_filename="";
}

#ifndef _METEOIO_JNI
extern "C"
{
#define COMPILE_PLUGIN
#include "exports.h"

	METEOIO_EXPORT void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	METEOIO_EXPORT void* loadObject(const string& classname, const Config& cfg) {
		if(classname == "GRIBIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GRIBIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
