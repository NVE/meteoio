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

#include <meteoio/ResamplingAlgorithms2D.h>

#include <errno.h>
#include <grib_api.h>

using namespace std;

namespace mio {
/**
 * @page gribio GRIBIO
 * @section gribio_format Format
 * *Put here the informations about the standard format that is implemented*
 *
 * @section gribio_units Units
 *
 *
 * @section gribio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

const double GRIBIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const unsigned int GRIBIO::MAX_VAL_LEN = 1024; //max value string lengthin GRIB

GRIBIO::GRIBIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	setOptions();
}

GRIBIO::GRIBIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	setOptions();
}

GRIBIO::GRIBIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	setOptions();
}

GRIBIO::~GRIBIO() throw()
{

}

void GRIBIO::setOptions()
{
	cfg.getValue("TIME_ZONE", "Input", tz_in);
	fp = NULL;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	grid2dpath_in = "/local";
}

void GRIBIO::listKeys(grib_handle** h, const std::string& filename)
{
	unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;
	//char* name_space=(char *)"ls"; //name_space=NULL to get all the keys, char* name_space=0;
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

void GRIBIO::listContent_old(const std::string& filename)
{
	grib_handle* h=NULL;
	int err=0, grib_count=0;

	//grib_multi_support_on(0);

	//loop over the messages
	while((h = grib_handle_new_from_file(0,fp,&err)) != NULL) {
		grib_count++;
		//std::cout << "GRIB N. " << grib_count << "\n";
		if(!h) {
			cleanup();
			throw IOException("Unable to create grib handle for \""+filename+"\"", AT);
		}

		long table_id;
		GRIB_CHECK(grib_get_long(h,"table2Version",&table_id),0);
		long param_id;
		GRIB_CHECK(grib_get_long(h,"indicatorOfParameter",&param_id),0);
		long levelType;
		GRIB_CHECK(grib_get_long(h,"indicatorOfTypeOfLevel", &levelType),0); //sfc (surface), pl (pressure level), ml (model level)

		long dataDate, dataTime;
		GRIB_CHECK(grib_get_long(h,"dataDate",&dataDate),0);
		GRIB_CHECK(grib_get_long(h,"dataTime",&dataTime),0);
		const int year=dataDate/10000, month=dataDate/100-year*100, day=dataDate-month*100-year*10000;
		const int hour=dataTime/100, minutes=dataTime-hour*100;
		Date date(year, month, day, hour, minutes, tz_in);

		if(table_id==201 && param_id==133) {
			std::cout << table_id << ":" << param_id << " -> Date=" << date.toString(Date::ISO) << " ";
			long Ni, Nj;
			GRIB_CHECK(grib_get_long(h,"Nx",&Ni),0);
			GRIB_CHECK(grib_get_long(h,"numberOfPointsAlongAMeridian",&Nj),0);
			std::cout << Ni << "x" << Nj << "\n";
			long nb_pts;
			GRIB_CHECK(grib_get_long(h,"numberOfDataPoints",&nb_pts),0);
			std::cout << "nb_pts=" << nb_pts << " ";
			std::cout << "\n";

			/*double *values;
			size_t values_len= 0;
			GRIB_CHECK(grib_get_size(h,"values",&values_len),0);
			values = (double*)malloc(values_len*sizeof(double));

			GRIB_CHECK(grib_get_double_array(h,"values",values,&values_len),0);
			for(size_t i = 0; i < values_len; i++)
				printf("%d %g\n",i+1,values[i]);
			free(values);*/

			//listKeys(&h, filename);
		} else {
			//std::cout << table_id << ":" << param_id << " -> Date=" << date.toString(Date::ISO) << "\n";
		}
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
{
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
	//Grid2DObject tmp_grid(static_cast<unsigned int>(Ni), static_cast<unsigned int>(Nj), cellsize_x, llcorner);
	grid_out.set(Ni, Nj, cellsize_x, llcorner); //HACK
	int i=0;
	for(unsigned int jj=0; jj<(unsigned)Nj; jj++) {
		for(unsigned int ii=0; ii<(unsigned)Ni; ii++)
			//tmp_grid(ii,jj) = values[i++];
			grid_out(ii,jj) = values[i++];
	}
	free(values);
	//grid_out = ResamplingAlgorithms2D::BilinearResampling(tmp_grid, 1., cellsize_x/cellsize_y);
}

void GRIBIO::read2DGrid_intern(const std::string& filename, const double& in_marsParam, const long& i_levelType, const long& i_level, const Date i_date, Grid2DObject& grid_out)
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

		if(marsParam==in_marsParam) {
			const Date date = getDate(h);

			long levelType;
			GRIB_CHECK(grib_get_long(h,"indicatorOfTypeOfLevel", &levelType),0); //sfc (surface), pl (pressure level), ml (model level)
			long level=0;
			if(levelType==105) GRIB_CHECK(grib_get_long(h,"level", &level),0);

			if(levelType==i_levelType && level==i_level && date==i_date) {
				read2Dlevel(h, grid_out);
				return;
			}
		}
	}
}

void GRIBIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*i_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GRIBIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);

	//in 115.260 and 116.260, swiss coordinates for each cell
	const std::string prefix="laf";
	const std::string ext=".grb";
	const std::string filename = grid2dpath_in+"/"+prefix+date.toString(Date::NUM).substr(0,10)+"f"+date.toString(Date::NUM).substr(10,2)+ext;

	fp = fopen(filename.c_str(),"r");
	if(fp==NULL) {
		stringstream ss;
		ss << "Error openning file \"" << filename << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	if(parameter==MeteoGrids::DW) read2DGrid_intern(filename, 31.2, 105, 10, date, grid_out); //10m wind direction, level 10, type 105
	if(parameter==MeteoGrids::VW) read2DGrid_intern(filename, 32.2, 105, 10, date, grid_out); //10m wind speed, level 10, type 105
	if(parameter==MeteoGrids::TA) read2DGrid_intern(filename, 17.2, 105, 2, date, grid_out); //2m TA, type 105, level 2
	if(parameter==MeteoGrids::RH) read2DGrid_intern(filename, 52.2, 105, 2, date, grid_out); //type 105, level 2
	if(parameter==MeteoGrids::TSS) read2DGrid_intern(filename, 11.2, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::HNW) read2DGrid_intern(filename, 61.2, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::ILWR) read2DGrid_intern(filename, 25.201, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::ISWR) read2DGrid_intern(filename, 111.250, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::P) read2DGrid_intern(filename, 1.2, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::HS) read2DGrid_intern(filename, 66.2, 1, 0, date, grid_out); //type 1

	if(parameter==MeteoGrids::DEM) read2DGrid_intern(filename, 8.2, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::SLOPE) read2DGrid_intern(filename, 98.202, 1, 0, date, grid_out); //type 1
	if(parameter==MeteoGrids::AZI) read2DGrid_intern(filename, 99.202, 1, 0, date, grid_out); //type 1

	cleanup();
}

void GRIBIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
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
