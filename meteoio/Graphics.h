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
#ifndef __GRAPHICS_H__
#define __GRAPHICS_H__

#include <meteoio/Array2D.h>

#include <string>
#include <vector>

namespace mio {

/**
 * @class Legend
 * @brief This creates a legend as pixels in a Grid2DObject.
 * This should be used with/by a plugin that would then convert this Grid2DObject into a true
 * graphic file (png, etc).
 *
 * @ingroup graphics
 * @author Mathias Bavay
 * @date   2011-12-23
 */
class legend {
	public:
		legend(const unsigned int &height, const double &minimum, const double &maximum);
		static double getLegendWidth();
		const Array2D<double> getLegend();

		static const int bg_color; ///<marker for solid background
		static const int text_color; ///<marker for solid text

	private:
		Array2D<double> grid;
		void drawLegend(const unsigned int &height, const double &minimum, const double &maximum);
		void writeLine(const double& val, const unsigned int& px_row);
		void writeChar(const unsigned int i_char[10][6], const unsigned int& px_col, const unsigned int& px_row);

		static const unsigned int char_width;
		static const unsigned int char_height;

		static const unsigned int text_chars_nb;
		static const unsigned int char_space;
		static const unsigned int text_width;
		static const unsigned int sample_width;
		static const unsigned int sample_text_space;
		static const unsigned int legend_plot_space;
		static const unsigned int total_width;

		static const unsigned int interline;
		static const unsigned int label_height;
		static const unsigned int nb_labels;
		static const unsigned int total_height;

		static const unsigned int font_0[10][6], font_1[10][6], font_2[10][6], font_3[10][6], font_4[10][6];
		static const unsigned int font_5[10][6], font_6[10][6], font_7[10][6], font_8[10][6], font_9[10][6];
		static const unsigned int font_plus[10][6], font_minus[10][6], font_dot[10][6], font_E[10][6];
};

namespace Color {
	void RGBtoHSV(const double r, const double g, const double b, double &h, double &s, double &v);
	void HSVtoRGB(const double h, const double s, const double v, double &r, double &g, double &b);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Gradient class
/////////////////////////////////////////////////////////////////////////////////////////////////
class Gradient_model {
	public:
		Gradient_model() {setMinMax(0., 0.);}; //do not use this constructor!
		Gradient_model(const double& i_min, const double& i_max) {setMinMax(i_min, i_max);};
		//setBgColor()
		//setFgColor()

		//val must be between 0 and 1 -> check + in doc? TODO
		virtual void getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a) = 0;
	protected:
		double getInterpol(const double& val, const std::vector<double>& X, const std::vector<double>& Y);
		void setMinMax(const double& i_min, const double& i_max);
		void HSV2RGB(const double& h, const double& s, const double& v, unsigned char &r, unsigned char &g, unsigned char &b);

		double max_val, min_val, delta_val;
		std::vector<double> X, v_h,v_s,v_v; ///<control points: vector of X and associated hues, saturations and values. They must be in X ascending order
};

class heat_gradient : public Gradient_model {
	public:
		heat_gradient(const double& i_min, const double& i_max) {setMinMax(i_min, i_max);};
		void getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a);
};

class terrain_gradient : public Gradient_model {
	public:
		terrain_gradient(const double& i_min, const double& i_max);
		void getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a);
};

class slope_gradient : public Gradient_model {
	public:
		slope_gradient(const double& i_min, const double& i_max);
		void getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a);
};

class Gradient {
	public:
		/// This enum provides names for possible color gradients
		typedef enum TYPE {
		            terrain,
		            slope,
		            heat,
		            water
		} Type;

		Gradient() {model=NULL; delta_val=0.;}; //do not use this empty constructor!
		Gradient(const Type& type, const double& min_val, const double &max_val);
		~Gradient() {delete model;};
		//setBgColor()
		//setFgColor()

		//val must be between 0 and 1 -> check + in doc? TODO
		void getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a);

	private:
		double delta_val;
		Gradient_model *model;
};

} //namespace
#endif
