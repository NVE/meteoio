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

#include <netcdf.h>
#include <string>

namespace mio {

/**
 * @class NetCDFIO
 * @brief This (empty) class is to be used as a template for developing new plugins
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
		void get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data);
		void copy_data(const int& ncid, const std::map<std::string, size_t>& map_parameters, const std::map<std::string, double*> map_data, 
		               const size_t& number_of_stations, const size_t& number_of_records, std::vector< std::vector<MeteoData> >& vecMeteo);
		void readData(const int& ncid, const size_t& index_start, const std::vector<Date>& vec_date, const std::map<std::string, size_t>& map_parameters, const MeteoData& meteo_data, std::vector< std::vector<MeteoData> >& vecMeteo);
		void readMetaData(const int& ncid, std::vector<StationData>& vecStation);
		void copy_grid(const size_t& latlen, const size_t& lonlen, double*& lat, double*& lon, double*& grid, Grid2DObject& grid_out);
		std::string get_varname(const MeteoGrids::Parameters& parameter);
		void get_indices(const int& ncid, const Date& dateStart, const Date& dateEnd, size_t& indexStart, size_t& indexEnd, std::vector<Date>& vecDate);
		void calculate_offset(const std::string& units, NetCDFIO::TimeUnit& time_unit, Date& offset);
		void check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
		                       int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void open_file(const std::string& filename, const int& omode, int& ncid);
		void create_file(const std::string& filename, const int& cmode, int& ncid);
		void create_latlon_dimensions(const int& ncid, const Grid2DObject& grid, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon);
		void create_time_dimension(const int& ncid, int& did_time, int& vid_time);
		void get_variable(const int& ncid, const std::string& varname, int& varid);
		bool check_variable(const int& ncid, const std::string& varname);
		bool check_dim_var(const int& ncid, const std::string& dimname);
		size_t get_1D_var_len(const int& ncid, const std::string& varname);
		void get_dimension(const int& ncid, const std::string& dimname, int& dimid);
		void get_dimension(const int& ncid, const std::string& dimname, int& dimid, size_t& dimlen);
		void get_dimension(const int& ncid, const std::string& varname, const int& varid, 
		                   std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen);
		void get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, std::string& attr_value);
		void read_data(const int& ncid, const std::string& varname, const int& varid,
		               const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data);
		void read_data_2D(const int& ncid, const std::string& varname, const int& varid,
		                  const size_t& record, const size_t& count, const size_t& length, double*& data);
		void read_value(const int& ncid, const std::string& varname, const int& varid, double& data);
		void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data);
		void write_data(const int& ncid, const std::string& varname, const int& varid, double*& data);
		void write_data(const int& ncid, const std::string& varname, const int& varid, const Grid2DObject& grid, const size_t& pos_start, double*& data);
		size_t find_record(const int& ncid, const std::string& varname, const int& varid, const double& data);
		size_t append_record(const int& ncid, const std::string& varname, const int& varid, const double& data);
		void define_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid);
		void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
		void add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid);
		void add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid);
		void add_2D_record(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record, const int& dimid1, const int& dimid2, int& varid);
		void add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname);
		void start_definitions(const std::string& filename, const int& ncid);
		void end_definitions(const std::string& filename, const int& ncid);
		void close_file(const std::string& filename, const int& ncid);
		void read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name, const std::string& varname);
		double calculate_cellsize(const size_t& latlen, const size_t& lonlen, 
                                    double* const& lat, double* const& lon, double& factor);
		void calculate_dimensions(const Grid2DObject& grid, double*& lat_array, double*& lon_array);
		void fill_data(const Grid2DObject& grid, double*& data);

		const Config cfg;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const std::string lat_str, lon_str, z_str, ta_str, rh_str;
		static const std::string cf_time, cf_units, cf_days, cf_seconds;
		static const std::string cnrm_altitude, cnrm_aspect, cnrm_slope, cnrm_ta, cnrm_rh, cnrm_vw, cnrm_dw, cnrm_qair;
		static const std::string cnrm_co2air, cnrm_theorsw, cnrm_neb, cnrm_hnw, cnrm_snowf, cnrm_swr_direct, cnrm_swr_diffuse, cnrm_p, cnrm_ilwr;

		static std::map<std::string, size_t> paramname; ///<Associate a name with meteo parameters in Parameters
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map

		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		double in_dflt_TZ, out_dflt_TZ;     //default time zones
		std::vector<StationData> vecMetaData;
};

} //namespace
#endif