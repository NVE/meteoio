/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __PNGIO_H__
#define __PNGIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>
#include <meteoio/Graphics.h>

#include <string>
#include <png.h>

namespace mio {
/**
 * @class PNGIO
 * @brief This (empty) class is to be used as a template for developing new plugins
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2010-06-14
 */
class PNGIO : public IOInterface {
	public:
		PNGIO(void (*delObj)(void*), const Config& i_cfg);

		PNGIO(const std::string& configfile);
		PNGIO(const PNGIO&);
		PNGIO(const Config& cfgreader);
		~PNGIO() throw();

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
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void setOptions();
		void parse_size(const std::string& size_spec, unsigned int& width, unsigned int& height);
		double getScaleFactor(const double& grid_w, const double& grid_h);
		void createMetadata(const Grid2DObject& grid);
		void writeMetadata(png_structp &png_ptr, png_infop &info_ptr);
		Grid2DObject scaleGrid(const Grid2DObject& grid_in);
		void setFile(const std::string& filename, FILE *fp, png_structp& png_ptr, png_infop& info_ptr, const unsigned int &width, const unsigned int &height);
		unsigned int setLegend(const unsigned int &ncols, const unsigned int &nrows, const double &min, const double &max, Array2D<double> &legend_array);
		void writeDataSection(const Grid2DObject &grid, const Array2D<double> &legend_array, const Gradient &gradient, const unsigned int &full_width, png_structp &png_ptr);
		void cleanup(FILE *fp, png_structp png_ptr, png_infop info_ptr);
		std::string decimal_to_dms(const double& decimal);

		Config cfg;
		bool autoscale;
		bool has_legend;
		std::string scaling;
		unsigned int min_w, min_h, max_w, max_h;
		//plus bg and fg colors

		std::vector<std::string> metadata_key, metadata_text;

		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
};

} //namespace
#endif
