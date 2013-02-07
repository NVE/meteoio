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
 * @brief This plugin write 2D grids as PNG images.
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2012-02-26
 */
class PNGIO : public IOInterface {
	public:
		PNGIO(const std::string& configfile);
		PNGIO(const PNGIO&);
		PNGIO(const Config& cfgreader);
		~PNGIO() throw();

		PNGIO& operator=(const PNGIO&); ///<Assignement operator, required because of pointer member

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
		void setFile(const std::string& filename, png_structp& png_ptr, png_infop& info_ptr, const unsigned int &width, const unsigned int &height);
		void writeWorldFile(const Grid2DObject& grid_in, const std::string& filename);
		unsigned int setLegend(const unsigned int &ncols, const unsigned int &nrows, const double &min, const double &max, Array2D<double> &legend_array);
		void writeDataSection(const Grid2DObject &grid, const Array2D<double> &legend_array, const Gradient &gradient, const unsigned int &full_width, const png_structp &png_ptr, png_infop& info_ptr);
		void setPalette(const Gradient &gradient, png_structp& png_ptr, png_infop& info_ptr, png_color *palette);
		void closePNG(png_structp& png_ptr, png_infop& info_ptr, png_color *palette);
		std::string decimal_to_dms(const double& decimal);

		const Config cfg;
		FILE *fp; //since passing fp always fail...
		bool autoscale;
		bool has_legend;
		bool has_world_file; ///< create world file with each file?
		bool optimize_for_speed; ///< optimize for speed instead of compression?
		bool indexed_png; ///< write an indexed png?
		unsigned char nr_levels; ///< number of levels to represent? (less-> smaller file size and faster)
		std::string coordout, coordoutparam; //projection parameters
		std::string grid2dpath;

		std::string scaling;
		unsigned int min_w, min_h, max_w, max_h;

		std::vector<std::string> metadata_key, metadata_text;

		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const unsigned char channel_depth;
		static const unsigned char channel_max_color;
		static const unsigned char transparent_grey;
};

} //namespace
#endif
