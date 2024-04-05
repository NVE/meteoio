// SPDX-License-Identifier: LGPL-3.0-or-later
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
#ifndef GRIBIO_H
#define GRIBIO_H

#include <meteoio/IOInterface.h>
#include <meteoio/plugins/libcodes.h>

#include <string>

namespace mio {

using namespace codes;
/**
 * @class GRIBIO
 * @brief This plugin reads GRIB 1 or 2 data files
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2012-01-25
 */
class GRIBIO : public IOInterface {
	public:
		GRIBIO(const std::string& configfile);
		GRIBIO(const Config& cfgreader);

		virtual bool list2DGrids(const Date& /*start*/, const Date& /*end*/, std::map<Date, std::set<size_t> >& /*list*/) {return false;}
		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);

		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);
		
	private:
		void setOptions();

		Coords getGeolocalization(double &cellsize_x, double &cellsize_y, const std::map<std::string,double> &gridParams);
		void read2Dlevel(CodesHandlePtr &h, Grid2DObject& grid_out);
		bool read2DGrid_indexed(const std::string& in_paramId, const long& i_levelType, const long& i_level, const Date i_date, Grid2DObject& grid_out, const long &ensembleNumber=0);
		void read2DGrid(const std::string& filename, Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		void readWind(const std::string& filename, const Date& date);

		void readStations(std::vector<Coords> &vecPoints);

		void scanMeteoPath();

		bool removeDuplicatePoints(std::vector<Coords>& vecPoints, std::vector<double> &lats, std::vector<double> &lons);
		bool readMeteoMeta(std::vector<Coords>& vecPoints, std::vector<StationData> &stations, std::vector<double> &lats, std::vector<double> &lons);
		bool readMeteoValues(const std::string& paramId, const long& levelType, const long& i_level, const Date& i_date, const size_t& npoints, std::vector<double>& lats, std::vector<double>& lons, std::vector<double>& values);
		void fillMeteo(std::vector<double> &values, const MeteoData::Parameters& param, const size_t& npoints, std::vector<MeteoData> &Meteo);
		void readMeteoStep(std::vector<StationData> &stations, std::vector<double> &lats, std::vector<double> &lons, const Date i_date, std::vector<MeteoData> &Meteo);

		const Config cfg;
		std::string grid2dpath_in;
		std::string meteopath_in;
		std::vector<Coords> vecPts; //points to use for virtual stations if METEO=GRIB
		std::vector< std::pair<Date,std::string> > cache_meteo_files; //cache of meteo files in METEOPATH
		std::string meteo_ext; //file extension
		std::string grid2d_ext; //file extension
		std::string grid2d_prefix; //filename prefix, like "laf"
		std::string coordin, coordinparam; //projection parameters
		Grid2DObject VW, DW; //for caching wind fields, since they require quite some calculations
		Date wind_date;
		Coords llcorner;

		CodesIndexPtr file_index; //because it needs to be kept between calls
		std::vector<std::string> paramIdList;
		std::vector<long> ensembleNumbers;
		std::vector<long> levelTypes;
		double latitudeOfNorthernPole, longitudeOfNorthernPole; //for rotated coordinates
		double bearing_offset; //to correct vectors coming from rotated lat/lon, we will add an offset to the bearing
		double cellsize, factor_x, factor_y;

		static const std::string default_ext;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const double tz_in; //GRIB time zone
		bool meteo_initialized; //set to true after we scanned METEOPATH, filed the cache, read the virtual stations from io.ini
		bool llcorner_initialized; //set to true after we properly computed llcorner
		bool update_dem, indexed;
		bool warned_date_in_file;

};

} //namespace
#endif
