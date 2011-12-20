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

#include <string>
#include <png.h>

namespace mio {

class legend {
	public:
		legend(const unsigned int &height, const double &minimum, const double &maximum);
		const Array2D<double> getLegend();

	private:
		Array2D<double> grid;
		void writeLine(const double& val, const unsigned int& px_row);
		void writeChar(const unsigned int i_char[10][6], const double& color, const unsigned int& px_col, const unsigned int& px_row);

		static const unsigned int text_chars_nb; //each label will contain 9 chars
		static const unsigned int char_width; //3 pixels wide + 1 pixel space
		static const unsigned int text_width; //nb chars, 3 pixels wide + 1 pixel space
		static const unsigned int sample_width; //color sample 2 chars wide
		static const unsigned int sample_text_space;
		static const unsigned int legend_plot_space;
		static const unsigned int total_width;

		static const unsigned int char_height;
		static const unsigned int interline;
		static const unsigned int label_height; //1 char + 2 pixels interline
		static const unsigned int nb_labels; //every decile + 0 level
		static const unsigned int total_height;

		static const unsigned int font_0[10][6], font_1[10][6], font_2[10][6], font_3[10][6], font_4[10][6];
		static const unsigned int font_5[10][6], font_6[10][6], font_7[10][6], font_8[10][6], font_9[10][6];
		static const unsigned int font_plus[10][6], font_minus[10][6], font_dot[10][6], font_E[10][6];
};

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

	private:
		void setRGB(double val, const double& min, const double& max, png_byte *ptr);
		void RGBtoHSV(const double r, const double g, const double b, double &h, double &s, double &v);
		void HSVtoRGB(const double h, const double s, const double v, double &r, double &g, double &b);
		void writeMetadata(png_structp &png_ptr, png_infop &info_ptr);
		void cleanup(FILE *fp, png_structp png_ptr, png_infop info_ptr, png_bytep row);
		void cleanup() throw();

		Config cfg;
		static const bool autoscale;
		static const bool has_legend;
		static const double factor;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //namespace
#endif
