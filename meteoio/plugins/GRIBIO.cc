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

#include <errno.h>
#include <grib_api.h>

//#include "GRIBIO_names.cc" //HACK

using namespace std;

namespace mio {
/**
 * @page template GRIBIO
 * @section template_format Format
 * *Put here the informations about the standard format that is implemented*
 *
 * @section template_units Units
 *
 *
 * @section template_keywords Keywords
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

void GRIBIO::listContent(const std::string& filename)
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
			/*int Ni_missing=1;
			grib_is_missing(h, "Ni", &Ni_missing);
			if(Ni_missing==1) std::cout << "Key Ni missing\n";
			else {*/
				//grib_get_native_type(h,options->print_keys[i].name,&(options->print_keys[i].type));
				/*char Ni[100]="";
				size_t len=100;
				grib_get_string(h, "Ni", Ni, &len);*/

				/*long Ni=-1;
				grib_get_long(h,"numberOfPointsAlongAParallel",&Ni);
				std::cout << Ni << " Ni ";*/
			//}
			/*long Ni, Nj;
			GRIB_CHECK(grib_get_long(h,"Nx",&Ni),0);
			GRIB_CHECK(grib_get_long(h,"numberOfPointsAlongAMeridian",&Nj),0);
			std::cout << Ni << "x" << Nj << "\n";*/
			/*long nb_pts;
			GRIB_CHECK(grib_get_long(h,"numberOfDataPoints",&nb_pts),0);
			std::cout << "nb_pts=" << nb_pts << " ";*/
			std::cout << "\n";

			/*double *values;
			size_t values_len= 0;
			GRIB_CHECK(grib_get_size(h,"values",&values_len),0);
			values = (double*)malloc(values_len*sizeof(double));

			GRIB_CHECK(grib_get_double_array(h,"values",values,&values_len),0);
			for(size_t i = 0; i < values_len; i++)
				printf("%d %g\n",i+1,values[i]);
			free(values);*/

			listKeys(&h, filename);
		} else {
			//std::cout << table_id << ":" << param_id << " -> Date=" << date.toString(Date::ISO) << "\n";
		}

		//printf("table2Version=%ld parameter=%ld\n",table_id, param_id);
		/*std::cout << table_id << ":" << param_id << " -> ";
		if(table_id==2)
			 std::cout << parm_table_dwd_002[param_id].name << " (" << parm_table_dwd_002[param_id].comment << ")";
		else if(table_id==201)
			 std::cout << parm_table_dwd_201[param_id].name << " (" << parm_table_dwd_201[param_id].comment << ")";
		else if(table_id==202)
			 std::cout << parm_table_dwd_202[param_id].name << " (" << parm_table_dwd_202[param_id].comment << ")";
		else if(table_id==203)
			 std::cout << parm_table_dwd_203[param_id].name << " (" << parm_table_dwd_203[param_id].comment << ")";
		else if(table_id==204)
			 std::cout << parm_table_dwd_204[param_id].name << " (" << parm_table_dwd_204[param_id].comment << ")";
		else if(table_id==205)
			 std::cout << parm_table_dwd_205[param_id].name << " (" << parm_table_dwd_205[param_id].comment << ")";
		else
			std::cout << "unknown: " << param_id << "\n";
		std::cout << "\n";*/

		//GRIB_CHECK(grib_get_size(h,"values",&values_len),0); //gridded data in "values" key
	}
}

void GRIBIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& i_name)
{
	const std::string filename = grid2dpath_in+"/"+i_name;
	fp = fopen(filename.c_str(),"r");
	if(fp==NULL) {
		stringstream ss;
		ss << "Error openning file \"" << filename << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	listContent(filename);
	cleanup();
}

void GRIBIO::read2DGrid(Grid2DObject& /*grid_out*/, const MeteoGrids::Parameters& /*parameter*/, const Date& /*date*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
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
