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
 *
 * @section gribio_units Units
 * As specified by WMO.
 *
 * @section gribio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords)
 * - COORDPARAM: extra coordinates parameters (see Coords)
 * - TIME_ZONE: time zone
 * - GRID2DPATH: path where to find the grids
 * - GRIB_DEM_UPDATE: recompute slope/azimuth from the elevations when reading a DEM (default=false,
 * that is we use the slope and azimuth included in the GRIB file)
 *
 */

const double GRIBIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

GRIBIO::GRIBIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	setOptions();
	indexed = false;
	idx=NULL;
	fp = NULL;
}

GRIBIO::GRIBIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	setOptions();
	indexed = false;
	idx=NULL;
	fp = NULL;
}

GRIBIO::GRIBIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	setOptions();
	indexed = false;
	idx=NULL;
	fp = NULL;
}

GRIBIO::~GRIBIO() throw()
{
	cleanup();
}

void GRIBIO::setOptions()
{
	std::string coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	cfg.getValue("TIME_ZONE", "Input", tz_in);
	cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
	update_dem = false;
	cfg.getValue("GRIB_DEM_UPDATE", "Input", update_dem, Config::nothrow);
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
	double angleOfRotationInDegrees;
	GRIB_CHECK(grib_get_double(h,"angleOfRotationInDegrees",&angleOfRotationInDegrees),0);
	if(angleOfRotationInDegrees!=0.) {
		throw InvalidArgumentException("Rotated grids not supported!", AT);
	}

	double ll_latitude;
	double ll_longitude;
	GRIB_CHECK(grib_get_double(h,"latitudeOfFirstGridPointInDegrees",&ll_latitude),0);
	GRIB_CHECK(grib_get_double(h,"longitudeOfFirstGridPointInDegrees",&ll_longitude),0);

	double latitudeOfSouthernPole;
	double longitudeOfSouthernPole;
	GRIB_CHECK(grib_get_double(h,"latitudeOfSouthernPoleInDegrees",&latitudeOfSouthernPole),0);
	GRIB_CHECK(grib_get_double(h,"longitudeOfSouthernPoleInDegrees",&longitudeOfSouthernPole),0);

	Coords llcorner(coordin, coordinparam);
	llcorner.setLatLon( ll_latitude+(90.+latitudeOfSouthernPole), ll_longitude+longitudeOfSouthernPole, 0.);

	double ur_latitude;
	double ur_longitude;
	GRIB_CHECK(grib_get_double(h,"latitudeOfLastGridPointInDegrees",&ur_latitude),0);
	GRIB_CHECK(grib_get_double(h,"longitudeOfLastGridPointInDegrees",&ur_longitude),0);

	const double cntr_latitude = .5*(ll_latitude+ur_latitude)+(90.+latitudeOfSouthernPole);

	double d_i, d_j;
	GRIB_CHECK(grib_get_double(h,"jDirectionIncrementInDegrees",&d_j),0);
	GRIB_CHECK(grib_get_double(h,"iDirectionIncrementInDegrees",&d_i),0);

	cellsize_x = Coords::lon_degree_lenght(cntr_latitude)*d_i;
	cellsize_y = Coords::lat_degree_lenght(cntr_latitude)*d_j;

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

bool GRIBIO::read2DGrid_indexed(grib_index *idx, const double& in_marsParam, const long& i_levelType, const long& i_level, const Date i_date, Grid2DObject& grid_out)
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
	}
	return false;
}

void GRIBIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*i_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
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
	std::string keys="marsParam:d,indicatorOfTypeOfLevel:l"; //indexing keys
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
	const std::string prefix="laf";
	const std::string ext=".grb";
	const std::string filename = grid2dpath_in+"/"+prefix+date.toString(Date::NUM).substr(0,10)+"f"+date.toString(Date::NUM).substr(10,2)+ext;

	read2DGrid(filename, grid_out, parameter, date);
}

void GRIBIO::read2DGrid(const std::string& filename, Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if(!indexed) {
		//the file has not yet been indexed
		indexFile(filename);
	} else if(idx_filename!=filename) {
		 //the file name changed, we have to re-index it
		cleanup();
		indexFile(filename);
	}

	//Basic meteo parameters
	if(parameter==MeteoGrids::P) {
		if(!read2DGrid_indexed(idx, 1.2, 1, 0, date, grid_out)) //PS
			read2DGrid_indexed(idx, 64.202, 1, 0, date, grid_out); //HMO3
	}
	if(parameter==MeteoGrids::TA) read2DGrid_indexed(idx, 11.2, 105, 2, date, grid_out); //T_2M
	if(parameter==MeteoGrids::RH) {
		if(!read2DGrid_indexed(idx, 52.2, 105, 2, date, grid_out)) { //RELHUM_2M
			Grid2DObject ta;
			read2DGrid_indexed(idx, 11.2, 105, 2, date, ta); //T_2M
			read2DGrid_indexed(idx, 17.2, 105, 2, date, grid_out); //TD_2M
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = Atmosphere::DewPointtoRh(grid_out(ii,jj), ta(ii,jj), true);
				}
			}
		}
	}
	if(parameter==MeteoGrids::TSS) read2DGrid_indexed(idx, 197.201, 111, 0, date, grid_out); //T_SO
	if(parameter==MeteoGrids::TSG) read2DGrid_indexed(idx, 11.2, 1, 0, date, grid_out); //T_G

	//hydrological parameters
	if(parameter==MeteoGrids::HNW) read2DGrid_indexed(idx, 61.2, 1, 0, date, grid_out); //tp
	if(parameter==MeteoGrids::ROT) read2DGrid_indexed(idx, 90.2, 112, 0, date, grid_out); //RUNOFF
	if(parameter==MeteoGrids::SWE) read2DGrid_indexed(idx, 65.2, 1, 0, date, grid_out); //W_SNOW

	//radiation parameters
	if(parameter==MeteoGrids::ALB) read2DGrid_indexed(idx, 84.2, 1, 0, date, grid_out); //ALB_RAD
	if(parameter==MeteoGrids::ILWR) read2DGrid_indexed(idx, 25.201, 1, 0, date, grid_out); //ALWD_S
	if(parameter==MeteoGrids::ISWR) {
		if(!read2DGrid_indexed(idx, 111.250, 1, 0, date, grid_out)) { //GLOB
			Grid2DObject diff;
			read2DGrid_indexed(idx, 23.201, 1, 0, date, diff); //diffuse rad, ASWDIFD_S
			read2DGrid_indexed(idx, 22.201, 1, 0, date, grid_out); //direct rad, ASWDIR_S
			grid_out.grid2D += diff.grid2D;
		}
	}

	//DEM parameters
	if(parameter==MeteoGrids::DEM) read2DGrid_indexed(idx, 8.2, 1, 0, date, grid_out); //HSURF
	if(parameter==MeteoGrids::SLOPE) read2DGrid_indexed(idx, 98.202, 1, 0, date, grid_out); //SLO_ANG
	if(parameter==MeteoGrids::AZI) read2DGrid_indexed(idx, 99.202, 1, 0, date, grid_out); //SLO_ASP

	//Wind parameters
	if(parameter==MeteoGrids::U) read2DGrid_indexed(idx, 33.2, 105, 10, date, grid_out); //U_10M, also in 110, 10 as U
	if(parameter==MeteoGrids::V) read2DGrid_indexed(idx, 34.2, 105, 10, date, grid_out); //V_10M, also in 110, 10 as V
	if(parameter==MeteoGrids::W) read2DGrid_indexed(idx, 40.2, 109, 10, date, grid_out); //W, 10m
	if(parameter==MeteoGrids::VW_MAX) read2DGrid_indexed(idx, 187.201, 105, 10, date, grid_out); //VMAX_10M 10m
	if(parameter==MeteoGrids::DW) {
		if(!read2DGrid_indexed(idx, 31.2, 105, 10, date, grid_out)) { //DD_10M
			Grid2DObject V;
			read2DGrid_indexed(idx, 34.2, 105, 10, date, V); //V_10M
			read2DGrid_indexed(idx, 33.2, 105, 10, date, grid_out); //U_10M
			const double to_deg = 180. / Cst::PI;
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = fmod( atan2( grid_out(ii,jj), V(ii,jj) ) * to_deg + 360., 360.); // turn into degrees [0;360)
				}
			}
		}
	}
	if(parameter==MeteoGrids::VW) {
		if(!read2DGrid_indexed(idx, 32.2, 105, 10, date, grid_out)) { //FF_10M
			Grid2DObject V;
			read2DGrid_indexed(idx, 34.2, 105, 10, date, V); //V_10M
			read2DGrid_indexed(idx, 33.2, 105, 10, date, grid_out); //U_10M
			for(unsigned int jj=0; jj<grid_out.nrows; jj++) {
				for(unsigned int ii=0; ii<grid_out.ncols; ii++) {
					grid_out(ii,jj) = sqrt( IOUtils::pow2(grid_out(ii,jj)) + IOUtils::pow2(V(ii,jj)) );
				}
			}
		}
	}
	/*if(parameter==MeteoGrids::U || parameter==MeteoGrids::V || parameter==MeteoGrids::W || parameter==MeteoGrids::VW || parameter==MeteoGrids::VW_MAX) {
		//we need to compute the wind at 7.5m
		Grid2DObject Z0;
		if(read2DGrid_indexed(idx, 83.2, 1, 0, date, Z0)) { //Z0
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
	const Date d;
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

void GRIBIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                             const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
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
	if(idx!=NULL) free(idx); idx=NULL;
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
