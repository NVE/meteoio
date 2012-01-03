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
#include "Graphics.h"

#include <cmath>

using namespace std;

namespace mio {

/////////////////////////////////////////////////////////////////////////////////////////////////
// Legend class
/////////////////////////////////////////////////////////////////////////////////////////////////

const int legend::bg_color = IOUtils::nodata-1;
const int legend::text_color = IOUtils::nodata-2;

const unsigned int legend::char_width = 6;
const unsigned int legend::char_height = 10;

const unsigned int legend::text_chars_nb = 9; //each label will contain 9 chars
const unsigned int legend::char_space = 1; //horizontal space between two chars
const unsigned int legend::text_width = legend::text_chars_nb*(legend::char_width+legend::char_space); //whole text line
const unsigned int legend::sample_width = legend::char_width*1; //color sample 2 chars wide
const unsigned int legend::sample_text_space = 6; //width of the color sample
const unsigned int legend::legend_plot_space = legend::char_width*1;
const unsigned int legend::total_width = legend::legend_plot_space+legend::sample_width+legend::sample_text_space+legend::text_width;

const unsigned int legend::interline = 5;
const unsigned int legend::label_height = legend::char_height+legend::interline; //1 char + interline
const unsigned int legend::nb_labels = 10+1; //every decile + 0 level
const unsigned int legend::total_height = legend::nb_labels*legend::label_height+legend::interline;

const unsigned int legend::font_0[char_height][char_width] = {{0,0,1,1,0,0}, {0,1,0,0,1,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,1,0,0,1,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_1[char_height][char_width] = {{0,0,0,1,0,0}, {0,0,1,1,0,0}, {0,1,1,1,0,0}, {1,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_2[char_height][char_width] = {{0,1,1,1,1,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,1,1,0}, {0,0,1,1,0,0}, {0,1,1,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,1,1}};
const unsigned int legend::font_3[char_height][char_width] = {{0,1,1,1,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,1,1,1,0}, {0,0,0,0,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,1,1,1,1,0}};
const unsigned int legend::font_4[char_height][char_width] = {{0,0,0,1,1,0}, {0,0,1,1,1,0}, {0,1,1,1,1,0}, {1,1,0,1,1,0}, {1,0,0,1,1,0}, {1,1,1,1,1,1}, {0,0,0,1,1,0}, {0,0,0,1,1,0}, {0,0,0,1,1,0}, {0,0,1,1,1,1}};
const unsigned int legend::font_5[char_height][char_width] = {{1,1,1,1,1,1}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,0,0}, {0,0,0,0,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {1,0,0,0,1,0}, {0,1,1,1,0,0}};
const unsigned int legend::font_6[char_height][char_width] = {{0,1,1,1,1,0}, {1,1,0,0,1,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}};
const unsigned int legend::font_7[char_height][char_width] = {{1,1,1,1,1,1}, {0,0,0,0,1,1}, {0,0,0,1,1,0}, {0,0,0,1,0,0}, {0,0,1,1,0,0}, {0,0,1,0,0,0}, {0,1,1,0,0,0}, {0,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}};
const unsigned int legend::font_8[char_height][char_width] = {{0,0,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}};
const unsigned int legend::font_9[char_height][char_width] = {{0,0,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,1,0,0,1,1}, {0,1,1,1,1,0}};
const unsigned int legend::font_plus[char_height][char_width] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {1,1,1,1,1,1}, {1,1,1,1,1,1}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
const unsigned int legend::font_minus[char_height][char_width] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {1,1,1,1,1,1}, {1,1,1,1,1,1}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
const unsigned int legend::font_dot[char_height][char_width] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_E[char_height][char_width] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,1,1,1,0,0}, {1,1,0,0,1,0}, {1,1,1,1,1,0}, {1,1,0,0,0,0}, {1,1,0,0,1,0}, {0,1,1,1,0,0}};

//create a legend of given height
//if height is insufficient, we don't generate any content, only transparent pixels
legend::legend(const unsigned int &height, const double &minimum, const double &maximum)
{
	drawLegend(height, minimum, maximum);
}

void legend::drawLegend(const unsigned int &height, const double &minimum, const double &maximum)
{
	grid.resize(total_width, height, IOUtils::nodata);
	const double level_inc = (maximum-minimum)/(double)(nb_labels-1); //the infamous interval thing...

	if(height>=total_height) {
		const unsigned int free_space = height-total_height;
		const unsigned int start_legend = free_space/2; //we will start from the bottom

		//fill the background
		for(unsigned int jj=start_legend; jj<(start_legend+total_height); jj++) {
			for(unsigned int ii=0; ii<total_width; ii++)
				grid(ii,jj) = bg_color;
		}

		for(unsigned int l=0; l<nb_labels; l++) {
			const double level_val = level_inc*l+minimum;
			const unsigned int px_row = l*label_height+start_legend;
			writeLine(level_val, px_row);
		}
	}
}

double legend::getLegendWidth() {
	return total_width;
}

void legend::writeLine(const double& val, const unsigned int& px_row)
{
	std::stringstream ss;
	const unsigned int precision = text_chars_nb-6; //full width - (sgn, dot, "e", sgn, two digits exponent)
	ss << setfill (' ') << setw(text_chars_nb) << left << setprecision(precision) << val << endl;

	const unsigned int x_offset = legend_plot_space+sample_width+sample_text_space;

	//write legend colored square
	for(unsigned int j=(px_row+interline); j<(px_row+label_height); j++) {
		for(unsigned int i=legend_plot_space; i<(legend_plot_space+sample_width); i++) {
			grid(i,j) = val;
		}
	}

	for(size_t i=0; i<ss.str().size(); i++) {
		char c=ss.str()[i];
		const unsigned int px_col = i*(char_width+char_space)+x_offset;
		if(c=='0') writeChar(font_0, px_col, px_row);
		if(c=='1') writeChar(font_1, px_col, px_row);
		if(c=='2') writeChar(font_2, px_col, px_row);
		if(c=='3') writeChar(font_3, px_col, px_row);
		if(c=='4') writeChar(font_4, px_col, px_row);
		if(c=='5') writeChar(font_5, px_col, px_row);
		if(c=='6') writeChar(font_6, px_col, px_row);
		if(c=='7') writeChar(font_7, px_col, px_row);
		if(c=='8') writeChar(font_8, px_col, px_row);
		if(c=='9') writeChar(font_9, px_col, px_row);
		if(c=='+') writeChar(font_plus, px_col, px_row);
		if(c=='-') writeChar(font_minus, px_col, px_row);
		if(c=='.') writeChar(font_dot, px_col, px_row);
		if(c=='e') writeChar(font_E, px_col, px_row);
		if(c=='E') writeChar(font_E, px_col, px_row);
	}
}

void legend::writeChar(const unsigned int i_char[char_height][char_width], const unsigned int& px_col, const unsigned int& px_row)
{
	for(unsigned int jj=0; jj<char_height; jj++) {
		for(unsigned int ii=0; ii<char_width; ii++) {
			const unsigned int char_px = i_char[char_height-1-jj][ii]; //we need to swap vertically each char
			if(char_px==0)
				grid(ii+px_col,jj+px_row+interline) = bg_color;
			else
				grid(ii+px_col,jj+px_row+interline) = text_color;
		}
	}
}

const Array2D<double> legend::getLegend()
{
	return grid;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Color class
/////////////////////////////////////////////////////////////////////////////////////////////////

//values between 0 and 1
//see http://www.cs.rit.edu/~ncs/color/t_convert.html or https://secure.wikimedia.org/wikipedia/en/wiki/HSL_and_HSV#Conversion_from_HSL_to_RGB
void Color::RGBtoHSV(const double r, const double g, const double b,
                     double &h, double &s, double &v)
{
	const double minimum = min( min(r,g), b);
	const double maximum = max( max(r,g), b);
	const double delta = maximum - minimum;

	v = maximum;

	if( maximum!=0 )
		s = delta/maximum;
	else { // r = g = b = 0 -> s = 0, v is undefined
		s = 0;
		h = -1;
		return;
	}

	if( r==maximum )
		h = ( g - b ) / delta;		// between yellow & magenta
	else if( g==maximum )
		h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		h = 4 + ( r - g ) / delta;	// between magenta & cyan

	h *= 60;				// degrees
	if( h < 0 )
		h += 360;
}

//values between 0 and 1
//h between 0 and 360
void Color::HSVtoRGB(const double h, const double s, const double v, double &r, double &g, double &b)
{
	if( s==0 ) {
		// achromatic (grey)
		r = g = b = v;
		return;
	}

	const double h_p = h/60;	// sector 0 to 5
	const int i = (int)floor(h_p); //HACK: replace by static_cast<int>
	const double f = h_p - i;	// factorial part of h
	const double p = v * ( 1 - s );
	const double q = v * ( 1 - s * f );
	const double t = v * ( 1 - s * ( 1 - f ) );

	switch( i ) {
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		default:		// case 5:
			r = v;
			g = p;
			b = q;
			break;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Gradient class
/////////////////////////////////////////////////////////////////////////////////////////////////

Gradient::Gradient(const Type& type, const double& i_min, const double& i_max)
{
	delta_val = i_max - i_min;

	if(type==terrain) model = new terrain_gradient(i_min, i_max);
	else if(type==slope) model = new slope_gradient(i_min, i_max);
	else if(type==heat) model = new heat_gradient(i_min, i_max);
	//else if(type==water) model = new water_gradient(i_min, i_max);
}

//val between min_val and max_val
//return values between 0 and 255 per channel
void Gradient::getColor(const double& val, unsigned char& r, unsigned char& g, unsigned char& b, unsigned char& a)
{
	if(val==IOUtils::nodata) {
		r=0; g=0; b=0; a=0;
		return;
	}
	if(val==legend::bg_color) {
		r=255; g=255; b=255; a=255;
		return;
	}
	if(val==legend::text_color) {
		r=0; g=0; b=0; a=255;
		return;
	}

	if(delta_val==0) {
		r=g=b=0; a=255;
		return;
	}

	//get the rgba components by providing val between 0 and 1
	model->getColor(val, r, g, b, a);
}

void Gradient_model::setMinMax(const double& i_min, const double& i_max)
{
	min_val = i_min;
	max_val = i_max;
	delta_val = (max_val-min_val);
}

//we assume that the vectors are sorted by X
double Gradient_model::getInterpol(const double& val, const std::vector<double>& X, const std::vector<double>& Y)
{
	if(X.size()!=Y.size()) {
		std::stringstream ss;
		ss << "For color gradients interpolations, both X and Y vectors must have the same size! ";
		ss << "There are " << X.size() << " abscissa for " << Y.size() << " ordinates.";
		throw IOException(ss.str(), AT);
	}
	size_t i=0;
	while(X[i]<val && i<X.size()) i++;

	if(X[i]==val) return Y[i];
	if(i==0) return Y[0];
	if(i==Y.size()) return Y[ Y.size()-1 ];

	const double y = Y[i-1] + (val-X[i-1])/(X[i]-X[i-1]) * (Y[i]-Y[i-1]);
	return y;
}

void Gradient_model::HSV2RGB(const double& h, const double& s, const double& v, unsigned char &r, unsigned char &g, unsigned char &b)
{
	double r_d, g_d, b_d;
	Color::HSVtoRGB(h, s, v, r_d, g_d, b_d);
	r = static_cast<unsigned char>(r_d*255);
	g = static_cast<unsigned char>(g_d*255);
	b = static_cast<unsigned char>(b_d*255);
}

void heat_gradient::getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a)
{
	const double auto_val = (val-min_val)/delta_val; //autoscale

	const double h = 240. * (1.-auto_val);
	const double v = auto_val*0.75+0.25;
	const double s = 1-auto_val*0.3;

	HSV2RGB(h, s, v, r, g, b);
	a = 255; //no alpha for valid values
}

terrain_gradient::terrain_gradient(const double& i_min, const double& i_max) {
	setMinMax(i_min, i_max);

	//write gradient control points
	X.push_back(-1.); v_h.push_back(198.); v_s.push_back(.50); v_v.push_back(.74); //sea, light blue
	X.push_back(0.); v_h.push_back(198.); v_s.push_back(.50); v_v.push_back(.74); //sea, light blue
	X.push_back(0.); v_h.push_back(144.); v_s.push_back(.58); v_v.push_back(.39); //sea level, dark green
	X.push_back(1200.); v_h.push_back(46.); v_s.push_back(.54); v_v.push_back(.86); //yellow
	X.push_back(2200.); v_h.push_back(4.); v_s.push_back(.71); v_v.push_back(.53); //dark red
	X.push_back(2700.); v_h.push_back(22.); v_s.push_back(.88); v_v.push_back(.41); //maroon
	X.push_back(3000.); v_h.push_back(0.); v_s.push_back(0.); v_v.push_back(.82); //white
}

void terrain_gradient::getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a)
{
	const double h = getInterpol(val, X, v_h);
	const double s = getInterpol(val, X, v_s);
	const double v = getInterpol(val, X, v_v);

	HSV2RGB(h, s, v, r, g, b);
	a = 255; //no alpha for valid values
}

slope_gradient::slope_gradient(const double& i_min, const double& i_max) {
	setMinMax(i_min, i_max);

	//write gradient control points
	X.push_back(0.); v_h.push_back(185.); v_s.push_back(.26); v_v.push_back(.91); //light blue
	X.push_back(25.); v_h.push_back(122.); v_s.push_back(.44); v_v.push_back(.91); //light green
	X.push_back(30.); v_h.push_back(60.); v_s.push_back(.44); v_v.push_back(.91); //light yellow
	X.push_back(35.); v_h.push_back(22.); v_s.push_back(.44); v_v.push_back(.91); //orange
	X.push_back(40.); v_h.push_back(0.); v_s.push_back(.44); v_v.push_back(.91); //red
	X.push_back(45.); v_h.push_back(0.); v_s.push_back(.58); v_v.push_back(.35); //dark red
	X.push_back(50.); v_h.push_back(0.); v_s.push_back(1.); v_v.push_back(0.); //black
}

void slope_gradient::getColor(const double &val, unsigned char &r, unsigned char &g, unsigned char &b, unsigned char &a)
{
	const double h = getInterpol(val, X, v_h);
	const double s = getInterpol(val, X, v_s);
	const double v = getInterpol(val, X, v_v);

	HSV2RGB(h, s, v, r, g, b);
	a = 255; //no alpha for valid values
}

} //namespace
