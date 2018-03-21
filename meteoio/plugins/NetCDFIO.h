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
#ifndef NetCDFIO_H
#define NetCDFIO_H

#include <meteoio/IOInterface.h>

#include <string>

namespace mio {

class ncParameters {
	public:
		enum Dimensions {NONE, TIME, LATITUDE, LONGITUDE, ALTITUDE, NORTHING, EASTING};
		
		ncParameters(const std::string& filename, const Config& cfg, const std::string& schema, const double& tz_in, const bool& i_debug=false);
		
		std::pair<Date, Date> getDateRange() const;
		std::set<size_t> getParams() const;
		std::vector<Date> getTimestamps() const {return vecTime;}
		Grid2DObject read2DGrid(const size_t& param, const Date& date) const;
		Grid2DObject read2DGrid(const std::string& varname) const;
		Grid2DObject readDEM() const;
		
		void write2DGrid(Grid2DObject grid_in, const Date& date, const std::string& schema);
		
	private:
		typedef struct VAR_ATTR {
			VAR_ATTR() : name(), standard_name(), long_name(), units(), height(IOUtils::nodata), param(IOUtils::npos) {};
			VAR_ATTR(const size_t& prm, const std::string& str1, const std::string& str2, const std::string& str3, const std::string& str4, const double& hgt)
			                     : name(str1), standard_name(str2), long_name(str3), units(str4), height(hgt), param(prm) {};
			std::string toString() const {std::ostringstream os; os << "["  << MeteoGrids::getParameterName(param) << " - " << name << " / " << standard_name << " / " << long_name << " , in " << units << " @ " << height << "]"; return os.str();};
  
			std::string name;
			std::string standard_name;
			std::string long_name;
			std::string units;
			double height;
			size_t param; //mapping to our MeteoGrids::Parameters
		} var_attr;

		typedef struct NC_VARIABLE {
			NC_VARIABLE() : attributes(), units(), scale(1.), offset(0.), nodata(IOUtils::nodata), varid(-1) {};
			NC_VARIABLE(const var_attr& attr, const std::string& i_units, const double& i_scale, const double& i_offset, const double& i_nodata, const int& i_varid)
			                   : attributes(attr),units(i_units), scale(i_scale), offset(i_offset), nodata(i_nodata), varid(i_varid) {};
			std::string toString() const {std::ostringstream os; os << "[" << varid << " - " << "\"" << attributes.name << "\", \"" << units << "\" - packing( *" << scale << ", +" << offset << "), nodata=" << nodata << "]"; return os.str();};
			
			var_attr attributes;
			std::string units;
			double scale, offset, nodata;
			int varid;
		} nc_variable;
		
		typedef struct DIM_ATTRIBUTES {
			DIM_ATTRIBUTES() : name(), long_name(), units(), type(NONE) {};
			DIM_ATTRIBUTES(const Dimensions& i_type, const std::string& i_name, const std::string& i_long_name, const std::string& i_units)
			                     : name(i_name), long_name(i_long_name), units(i_units), type(i_type) {};
			std::string toString() const {std::ostringstream os; os << name << " / " << long_name << " , in " << units; return os.str();};
  
			std::string name;
			std::string long_name;
			std::string units;
			Dimensions type;
		} dim_attributes;
		
		typedef struct NC_DIMENSION {
			NC_DIMENSION() : attributes(), length(0), dimid(-1), varid(-1), isUnlimited(false) {};
			NC_DIMENSION(const dim_attributes& attr) : attributes(attr), length(0), dimid(-1), varid(-1), isUnlimited(false) {};
			NC_DIMENSION(const dim_attributes& attr, const size_t& len, const int& i_dimid, const int& i_varid, const bool& unlimited)
			                     : attributes(attr), length(len), dimid(i_dimid), varid(i_varid), isUnlimited(unlimited) {};
			std::string toString() const {std::ostringstream os; os << "[ " << dimid << "/" << varid << " - " << attributes.toString() << ", length " << length; if (isUnlimited) os << ", unlimited"; os << "]"; return os.str();};
			
			dim_attributes attributes;
			size_t length;
			int dimid, varid;
			bool isUnlimited;
		} nc_dimension;
		
		static std::map< std::string, std::vector<ncParameters::var_attr> > initSchemasVars();
		static std::map< std::string, std::vector<ncParameters::dim_attributes> > initSchemasDims();
		static std::vector<ncParameters::var_attr>  initUserSchemas(const Config& i_cfg);
		static std::string getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name);
		static void getAttribute(const int& ncid, const int& value_id, const std::string& value_name, const std::string& attr_name, double& attr_value);
		static std::vector<double> readDimension(const int& ncid, const nc_dimension& dim);
		static void getTimeTransform(const std::string& time_units, const double& i_TZ, double &o_time_offset, double &o_time_multiplier);
		
		void initVariables(const int& ncid, const std::string& schema_name);
		void initDimensions(const int& ncid);
		
		Grid2DObject read2DGrid(const nc_variable& var, const size_t& time_pos, const bool& m2mm=false, const bool& reZero=false) const;
		std::vector<Date> readTimeDimension(const int& ncid, const nc_dimension& dim) const;
		const ncParameters::var_attr getSchemaAttributes(const std::string& var, const std::string& schema_name) const;
		double calculate_cellsize(double& factor_x, double& factor_y) const;
		void fill2DGrid(Grid2DObject& grid, const double data[], const double& nodata) const;
		
		void create_dimension(const int& ncid, const Dimensions& dimType, const std::string& schema);
		
		static std::map< std::string, std::vector<ncParameters::var_attr> > schemas_vars; ///< all the variables' attributes for all schemas
		static std::map< std::string, std::vector<ncParameters::dim_attributes> > schemas_dims; ///< all the dimensions' attributes for all schemas
		
		std::vector<ncParameters::var_attr> user_schemas; ///< all the variables' attributes for the user defined schema
		std::map<size_t, nc_variable> vars; ///< all the recognized variables for the selected schema_name and current file
		std::map<std::string, nc_variable> unknown_vars; ///< all the unrecognized variables for the current file, as map< name, nc_variable>
		std::vector<Date> vecTime;
		std::vector<double> vecLat, vecLon;
		std::map<ncParameters::Dimensions, nc_dimension> dimensions_map; ///< all the dimensions for the current schema, as found in the current file
		std::string file_and_path;
		std::string coordin, coordinparam;
		double TZ;
		int max_dimension;
		bool wrf_hacks, debug;
};

/**
 * @class NetCDFIO
 * @brief This plug-in allows reading and writing of NetCDF files for gridded data. IT IS NOT YET USABLE!.
 *
 * @ingroup plugins
 */
class NetCDFIO : public IOInterface {
	public:
		NetCDFIO(const std::string& configfile);
		NetCDFIO(const NetCDFIO&);
		NetCDFIO(const Config& cfgreader);

		virtual bool list2DGrids(const Date& start, const Date& end, std::map<Date, std::set<size_t> >& list);
		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void parseInputOutputSection();
		void scanMeteoPath(const std::string& in_path, const std::string& nc_ext, std::vector< std::pair<std::pair<Date,Date>, ncParameters> > &meteo_files);
		void cleanMeteoCache(std::vector< std::pair<std::pair<Date,Date>, ncParameters> > &meteo_files);
		
		const Config cfg;
		std::vector< std::pair<std::pair<Date,Date>, ncParameters> > cache_grid_files; //cache of grid files in GRID2DPATH
		std::vector<MeteoGrids::Parameters> available_params;
		std::string in_schema, out_schema, in_grid2d_path, in_nc_ext;
		double in_dflt_TZ, out_dflt_TZ;
		bool dem_altimeter, debug;
};

} //namespace
#endif
