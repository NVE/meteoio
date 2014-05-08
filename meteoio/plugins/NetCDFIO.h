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
#ifndef __NetCDFIO_H__
#define __NetCDFIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/meteostats/libinterpol1D.h>

#include <netcdf.h>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace mio {

/**
 * @class NetCDFIO
 * @brief This plug-in allows reading and writing of NetCDF files formatted according to CNRM standard.
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2014-03-13
 */
class NetCDFIO : public IOInterface {
	public:
		enum TimeUnit { seconds, hours, days };

		NetCDFIO(const std::string& configfile);
		NetCDFIO(const NetCDFIO&);
		NetCDFIO(const Config& cfgreader);
		~NetCDFIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readPOI(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void parseInputOutputSection();
		void create_parameters(const int& ncid, const int& did_time, const int& did_points, const size_t& number_of_records,
		                       const size_t& number_of_stations, std::map<size_t, std::string>& map_param_name,
		                       std::map<std::string, double*>& map_data_2D, std::map<std::string, int>& varid);
		void create_meta_data(const int& ncid, const int& did, std::map<std::string, double*>& map_data_1D, std::map<std::string, int>& varid);
		void get_parameters(const std::vector< std::vector<MeteoData> >& vecMeteo, std::map<size_t, std::string>& map_param_name,
		                    std::map<std::string, double*>& map_data_1D, double*& dates);
		void get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data);
		size_t get_dates(const std::vector< std::vector<MeteoData> >& vecMeteo, double*& dates);
		void copy_data(const size_t& number_of_stations, const size_t& number_of_records, const std::vector< std::vector<MeteoData> >& vecMeteo,
                         const std::map<size_t, std::string>& map_param_name, std::map<std::string, double*>& map_data_2D);
		void copy_data(const int& ncid, const std::map<std::string, size_t>& map_parameters, const std::map<std::string, double*> map_data,
		               const size_t& number_of_stations, const size_t& number_of_records, std::vector< std::vector<MeteoData> >& vecMeteo);
		void readData(const int& ncid, const size_t& index_start, const std::vector<Date>& vec_date, const std::map<std::string, size_t>& map_parameters,
		              const MeteoData& meteo_data, std::vector< std::vector<MeteoData> >& vecMeteo);
		void readMetaData(const int& ncid, std::vector<StationData>& vecStation);
		void get_meta_data_ids(const int& ncid, std::map<std::string, int>& map_vid);
		std::string get_varname(const MeteoGrids::Parameters& parameter);
		void get_indices(const int& ncid, const Date& dateStart, const Date& dateEnd, size_t& indexStart, size_t& indexEnd, std::vector<Date>& vecDate);
		void calculate_offset(const std::string& units, NetCDFIO::TimeUnit& time_unit, Date& offset);
		void check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
		                       int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name, const std::string& varname, const Date& date=Date());
		void write2DGrid_internal(const Grid2DObject& grid_in, const std::string& filename, const std::string& varname, const Date& date=Date());
		void fill_data(const Grid2DObject& grid, double*& data);
		void copy_grid(const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
		               const double * const grid, Grid2DObject& grid_out);
		double calculate_cellsize(const size_t& latlen, const size_t& lonlen, const double * const lat, const double * const lon,
		                          double& factor_x, double& factor_y);
		void calculate_dimensions(const Grid2DObject& grid, double*& lat_array, double*& lon_array);
		void add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname);
		void create_latlon_dimensions(const int& ncid, const Grid2DObject& grid, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void create_time_dimension(const int& ncid, int& did_time, int& vid_time);


		/* libnetcdf wrappers */

		//Opening, creating, closing dataset
		void open_file(const std::string& filename, const int& omode, int& ncid);
		void create_file(const std::string& filename, const int& cmode, int& ncid);
		void start_definitions(const std::string& filename, const int& ncid);
		void end_definitions(const std::string& filename, const int& ncid);
		void close_file(const std::string& filename, const int& ncid);

		//Adding variables
		void add_0D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, int& varid);
		void add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid);
		void add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid);
		void add_3D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record,
		                     const int& dimid1, const int& dimid2, int& varid);

		//Adding attributes
		void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
		void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value);
		void get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, std::string& attr_value);

		//Adding dimensions
		void add_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid);


		//Reading data from NetCDF file
		void read_data(const int& ncid, const std::string& varname, const int& varid,
		               const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data);
		void read_data_2D(const int& ncid, const std::string& varname, const int& varid,
		                  const size_t& record, const size_t& count, const size_t& length, double*& data);
		void read_value(const int& ncid, const std::string& varname, const int& varid, double& data);
		void read_value(const int& ncid, const std::string& varname, const int& varid, const size_t& pos, double& data);
		void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data);

		//Writing data to NetCDF file
		void write_data(const int& ncid, const std::string& varname, const int& varid, const double * const data);
		void write_data(const int& ncid, const std::string& varname, const int& varid, const Grid2DObject& grid,
		                const size_t& pos_start, const double * const data);

		//Dealing with variables that have dimension NC_UNLIMITED
		size_t find_record(const int& ncid, const std::string& varname, const int& varid, const double& data);
		size_t add_record(const int& ncid, const std::string& varname, const int& varid, const double& data);
		void write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& pos,
		                  const size_t& length, const double * const data);

		//Dealing with variables and dimensions
		bool check_dim_var(const int& ncid, const std::string& dimname);
		bool check_variable(const int& ncid, const std::string& varname);
		size_t get_1D_var_len(const int& ncid, const std::string& varname);
		void get_variable(const int& ncid, const std::string& varname, int& varid);
		void check_dimensions(const int& ncid, const std::string& varname, const int& varid, const std::vector<std::string>& names);
		void get_dimension(const int& ncid, const std::string& dimname, int& dimid);
		void get_dimension(const int& ncid, const std::string& dimname, int& dimid, size_t& dimlen);
		void get_dimension(const int& ncid, const std::string& varname, const int& varid,
		                   std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen);


		// Private variables
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const std::string cf_time, cf_units, cf_days, cf_seconds, cf_latitude, cf_longitude, cf_altitude, cf_ta, cf_rh, cf_p;
		static const std::string cnrm_points, cnrm_latitude, cnrm_longitude, cnrm_altitude, cnrm_aspect, cnrm_slope, cnrm_ta, cnrm_rh, cnrm_vw, cnrm_dw, cnrm_qair;
		static const std::string cnrm_co2air, cnrm_theorsw, cnrm_neb, cnrm_hnw, cnrm_snowf, cnrm_swr_direct, cnrm_swr_diffuse, cnrm_p, cnrm_ilwr, cnrm_timestep;

		static std::map<std::string, size_t> paramname; ///<Associate a name with meteo parameters in Parameters
		static std::map<std::string, std::string> map_name; ///Associate MeteoIO parameter names with CNRM parameter names
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map

		const Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		double in_dflt_TZ, out_dflt_TZ;     //default time zones
		std::vector<StationData> vecMetaData;
};

} //namespace
#endif
