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

const int legend::bg_color = IOUtils::nodata-1;
const int legend::text_color = IOUtils::nodata-2;

const unsigned int legend::text_chars_nb = 9; //each label will contain 9 chars
const unsigned int legend::char_width = 6+1; //6 pixels wide + 1 pixel space
const unsigned int legend::text_width = legend::text_chars_nb*legend::char_width; //whole text line
const unsigned int legend::sample_width = legend::char_width*1; //color sample 2 chars wide
const unsigned int legend::sample_text_space = 6;
const unsigned int legend::legend_plot_space = legend::char_width*1;
const unsigned int legend::total_width = legend::legend_plot_space+legend::sample_width+legend::sample_text_space+legend::text_width;

const unsigned int legend::char_height = 10;
const unsigned int legend::interline = 5;
const unsigned int legend::label_height = legend::char_height+legend::interline; //1 char + interline
const unsigned int legend::nb_labels = 11; //every decile + 0 level
const unsigned int legend::total_height = legend::nb_labels*legend::label_height+legend::interline;

const unsigned int legend::font_0[10][6] = {{0,0,1,1,0,0}, {0,1,0,0,1,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,1,0,0,1,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_1[10][6] = {{0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,1,1,1,0,0}, {1,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_2[10][6] = {{0,0,1,1,1,0}, {0,1,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,1,1,0}, {0,0,1,1,0,0}, {0,1,1,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,1,1}};
const unsigned int legend::font_3[10][6] = {{0,0,1,1,1,0}, {0,1,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,1,1,1,0}, {0,0,0,0,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,1,1,1,1,0}};
const unsigned int legend::font_4[10][6] = {{0,0,0,1,1,0}, {0,0,1,1,1,0}, {0,1,1,1,1,0}, {1,1,0,1,1,0}, {1,0,0,1,1,0}, {1,1,1,1,1,1}, {0,0,0,1,1,0}, {0,0,0,1,1,0}, {0,0,0,1,1,0}, {0,0,1,1,1,1}};
const unsigned int legend::font_5[10][6] = {{1,1,1,1,1,1}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,0,0}, {0,0,0,0,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {1,0,0,0,1,0}, {0,1,1,1,0,0}};
const unsigned int legend::font_6[10][6] = {{0,1,1,1,1,0}, {1,1,0,0,1,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}};
const unsigned int legend::font_7[10][6] = {{1,1,1,1,1,1}, {0,0,0,0,1,1}, {0,0,0,1,1,0}, {0,0,0,1,0,0}, {0,0,1,1,0,0}, {0,0,1,0,0,0}, {0,1,1,0,0,0}, {0,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}};
const unsigned int legend::font_8[10][6] = {{0,0,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}};
const unsigned int legend::font_9[10][6] = {{0,0,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,1,0,0,1,1}, {0,1,1,1,1,0}};
const unsigned int legend::font_plus[10][6] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {1,1,1,1,1,1}, {1,1,1,1,1,1}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
const unsigned int legend::font_minus[10][6] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {1,1,1,1,1,1}, {1,1,1,1,1,1}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
const unsigned int legend::font_dot[10][6] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_E[10][6] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,1,1,1,0,0}, {1,1,0,0,1,0}, {1,1,1,1,1,0}, {1,1,0,0,0,0}, {1,1,0,0,1,0}, {0,1,1,1,0,0}};

//create a legend of given height
//if hight is insufficient, we don't generate any content, only transparent pixels
legend::legend(const unsigned int &height, const double &minimum, const double &maximum)
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

void legend::writeLine(const double& val, const unsigned int& px_row)
{
	std::stringstream ss;
	ss << setfill (' ') << setw (9) << left << setprecision(3) << val << endl;

	const unsigned int x_offset = legend_plot_space+sample_width+sample_text_space;

	//write legend colored square
	for(unsigned int j=(px_row+interline); j<(px_row+label_height); j++) {
		for(unsigned int i=legend_plot_space; i<(legend_plot_space+sample_width); i++) {
			grid(i,j) = val;
		}
	}

	for(size_t i=0; i<ss.str().size(); i++) {
		char c=ss.str()[i];
		const unsigned int px_col = i*char_width+x_offset;
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

void legend::writeChar(const unsigned int i_char[10][6], const unsigned int& px_col, const unsigned int& px_row)
{
	for(unsigned int jj=0; jj<10; jj++) {
		for(unsigned int ii=0; ii<6; ii++) {
			const unsigned int char_px = i_char[9-jj][ii]; //we need to swap vertically each char
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


Gradient::Gradient(const Type& type, const double& i_min_val, const double& i_max_val)
{
	min_val = i_min_val;
	max_val = i_max_val;
	delta_val = (max_val-min_val);

	//if(type==Type::terrain) {}
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
	getHeat((val-min_val)/delta_val, r, g, b, a);
	//getTerrain((val-min_val)/delta_val, r, g, b, a);
}

//val must be between 0 and 1
void Gradient::getHeat(const double& val, unsigned char& r, unsigned char& g, unsigned char& b, unsigned char& a)
{
	const double h = 240. * (1.-val);
	const double v = val*0.75+0.25;
	const double s = 1-val*0.3;

	double r_d, g_d, b_d;
	Color::HSVtoRGB(h, s, v, r_d, g_d, b_d);
	r = static_cast<unsigned char>(r_d*255);
	g = static_cast<unsigned char>(g_d*255);
	b = static_cast<unsigned char>(b_d*255);
	a = 255; //no alpha for valid values
}

void Gradient::getTerrain(const double& val, unsigned char& r, unsigned char& g, unsigned char& b, unsigned char& a)
{
	std::vector<double> values;
	std::vector<double> v_h,v_s,v_v;

	values.push_back(0.); v_h.push_back(240.); v_s.push_back(1.); v_v.push_back(0.25);
	values.push_back(0.5); v_h.push_back(230.); v_s.push_back(0.9); v_v.push_back(0.5);
	values.push_back(1.); v_h.push_back(0.); v_s.push_back(0.7); v_v.push_back(1.);
	const double h = getInterpol(val, values, v_h);
	const double s = getInterpol(val, values, v_s);
	const double v = getInterpol(val, values, v_v);

	double r_d, g_d, b_d;
	Color::HSVtoRGB(h, s, v, r_d, g_d, b_d);
	r = static_cast<unsigned char>(r_d*255);
	g = static_cast<unsigned char>(g_d*255);
	b = static_cast<unsigned char>(b_d*255);
	a = 255; //no alpha for valid values
}

//we assume that the vectors are sorted by X
double Gradient::getInterpol(const double& val, const std::vector<double>& X, const std::vector<double>& Y)
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
	//if(i==Y.size()) return Y[ Y.size() ]; // not necessary, treated by the formula
	
	const double y = Y[i-1] + (val-X[i-1])/(X[i]-X[i-1]) * (Y[i]-Y[i-1]);
	return y;
}

} //namespace
