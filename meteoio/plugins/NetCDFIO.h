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
		void open_file(const std::string& filename, const int& omode, int& ncid);
		void create_file(const std::string& filename, const int& cmode, int& ncid);
		void get_variable(const int& ncid, const std::string& varname, int& varid);
		bool check_variable(const int& ncid, const std::string& varname);
		size_t get_1D_var_len(const int& ncid, const std::string& varname);
		void get_dimension(const int& ncid, const std::string& varname, const int& varid, 
		                   std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen);
		void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data);
		void write_data(const int& ncid, const std::string& varname, const int& varid, double*& data);
		void define_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid);
		void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
		void add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid);
		void add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid);
		void add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname);
		void start_definitions(const std::string& filename, const int& ncid);
		void end_definitions(const std::string& filename, const int& ncid);
		void close_file(const std::string& filename, const int& ncid);
		void read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name, const std::string& varname);
		double calculate_cellsize(const size_t& latlen, const size_t& lonlen, 
                                    double* const& lat, double* const& lon, double& factor);
		void calculate_dimensions(const Grid2DObject& grid, double*& lat_array, double*& lon_array);
		void fill_data(const Grid2DObject& grid, double*& data);
		void cleanup() throw();

		const Config cfg;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //namespace
#endif
