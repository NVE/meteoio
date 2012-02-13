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
#include "Graphics.h"

#include <cmath>

using namespace std;

namespace mio {

/////////////////////////////////////////////////////////////////////////////////////////////////
// Legend class
/////////////////////////////////////////////////////////////////////////////////////////////////

const double legend::bg_color = IOUtils::nodata-1;
const double legend::text_color = IOUtils::nodata-2;

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
const unsigned int legend::nb_labels = 12+1; //every decile + 0 level
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

//uniform division of the given range
void legend::simpleLegend(const unsigned int &height, const double &minimum, const double &maximum)
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

//this tries to properly round to legend min, max and steps so that it gets easier to follow
void legend::smartLegend(const unsigned int &height, const double &minimum, const double &maximum)
{
	grid.resize(total_width, height, IOUtils::nodata);

	const double range = maximum-minimum;
	double min_norm=minimum, max_norm=maximum, decade_mult;
	unsigned int step_norm, nb_labels_norm;

	if(range>0.) {
		const double decade = floor(log10(range));
		decade_mult = pow(10.,decade);
		min_norm = floor(minimum/decade_mult*10.)/10., max_norm = ceil(maximum/decade_mult*10.)/10.; //between 0 & 10
		const double range_norm = max_norm-min_norm; //between 0 & 10
		const double step = range_norm/nb_labels; //between 0 & 10 / number of labels -> between 0 & 1

		//for precision issues, it is necessary to keep step_norm as an interger -> we will always use step_norm/10
		if(step<=0.1) step_norm=1;
		else if(step<=0.2) step_norm=2;
		else if(step<=0.5) step_norm=5;
		else step_norm=10;

		min_norm = round(min_norm/step_norm*10)*step_norm/10.;
		max_norm = round(max_norm/step_norm*10)*step_norm/10.;
		nb_labels_norm = (unsigned int)ceil((max_norm-min_norm)/step_norm*10)+1; //because min_norm might have been tweaked
	} else {
		step_norm = 0;
		decade_mult=1.;
		nb_labels_norm=1;
	}

	const unsigned int smart_height = nb_labels_norm*legend::label_height+legend::interline;
	if(height>=smart_height) {
		const unsigned int free_space = height-smart_height;
		const unsigned int start_legend = free_space/2; //we will start from the bottom

		//fill the background
		for(unsigned int jj=start_legend; jj<(start_legend+smart_height); jj++) {
			for(unsigned int ii=0; ii<total_width; ii++)
				grid(ii,jj) = bg_color;
		}

		for(unsigned int l=0; l<nb_labels_norm; l++) {
			const double level_val = (step_norm*l/10.+min_norm)*decade_mult;
			const unsigned int px_row = l*label_height+start_legend;
			writeLine(level_val, px_row);
		}
	}
}

void legend::drawLegend(const unsigned int &height, const double &minimum, const double &maximum){
	smartLegend(height, minimum, maximum);
	//simpleLegend(height, minimum, maximum);
}

unsigned int legend::getLegendWidth() {
	return total_width;
}

void legend::writeLine(const double& val, const unsigned int& px_row)
{
	std::stringstream ss;
	//const unsigned int precision = text_chars_nb-6; //full width - (sgn, dot, "e", sgn, two digits exponent)
	//ss << std::setfill (' ') << std::setw(text_chars_nb) << std::left << std::setprecision(precision) << val << std::endl; //improve this format...
	ss << std::setfill (' ') << std::setw(text_chars_nb) << std::left << val << std::endl;

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

const Array2D<double> legend::getLegend() const {
	return grid;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Color class
/////////////////////////////////////////////////////////////////////////////////////////////////

//values between 0 and 1
//see http://www.cs.rit.edu/~ncs/color/t_convert.html or https://secure.wikimedia.org/wikipedia/en/wiki/HSL_and_HSV#Conversion_from_HSL_to_RGB
void Color::RGBtoHSV(const double& r, const double& g, const double& b,
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

void Color::HSVtoRGB(const double& h, const double& s, const double& v, double &r, double &g, double &b)
{
	if( s==0 ) { //achromatic (grey)
		r = g = b = v;
		return;
	}

	const double h_p = h/60.; //sector 0 to 5
	const int i = static_cast<int>(h_p);
	const double f = h_p - i; //factorial part of h
	const double p = v * ( 1. - s );
	const double q = v * ( 1. - s * f );
	const double t = v * ( 1. - s * ( 1. - f ) );

	switch( i ) {
		case 0:
			r = v; g = t; b = p;
			break;
		case 1:
			r = q; g = v; b = p;
			break;
		case 2:
			r = p; g = v; b = t;
			break;
		case 3:
			r = p; g = q; b = v;
			break;
		case 4:
			r = t; g = p; b = v;
			break;
		default: // case 5:
			r = v; g = p; b = q;
			break;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Gradient class
/////////////////////////////////////////////////////////////////////////////////////////////////
const unsigned char Gradient::channel_max_color = 255;
const unsigned char Gradient::reserved_idx = 5;
const unsigned char Gradient::reserved_cols = 2;

Gradient::Gradient(const Type& type, const double& i_min, const double& i_max, const bool& i_autoscale)
{
	set(type, i_min, i_max, i_autoscale);
	nr_unique_cols = 0;
}

void Gradient::set(const Type& type, const double& i_min, const double& i_max, const bool& i_autoscale)
{
	min = i_min;
	max = i_max;
	delta = i_max - i_min;
	autoscale = i_autoscale;

	if(type==terrain) model = new gr_terrain(i_min, i_max, i_autoscale);
	else if(type==slope) model = new gr_slope(i_min, i_max, i_autoscale);
	else if(type==azi) model = new gr_azi(i_min, i_max, i_autoscale);
	else if(type==heat) model = new gr_heat(i_min, i_max, i_autoscale);
	else if(type==freeze) model = new gr_freeze(i_min, i_max, i_autoscale);
	else if(type==blue) model = new gr_blue(i_min, i_max, i_autoscale);
	else if(type==blue_pink) model = new gr_blue_pink(i_min, i_max, i_autoscale);
	else if(type==pastel) model = new gr_pastel(i_min, i_max, i_autoscale);
	else if(type==bg_isomorphic) model = new gr_bg_isomorphic(i_min, i_max, i_autoscale);
}

void Gradient::setNrOfLevels(const unsigned char& i_nr_unique_levels) {
	if(i_nr_unique_levels<=reserved_idx) {
		stringstream ss;
		ss << "Insufficient number of colors requested for gradient: ask for more than ";
		ss << reserved_idx << " colors!";
		throw InvalidArgumentException(ss.str(), AT);
	}
	nr_unique_cols = i_nr_unique_levels - reserved_idx;
}

//val between min_val and max_val
//return values between 0 and 255 per channel
//getColor: take value between min & max, as defined in the constructor (use min/max for rescaling gradient control points). Use autoscale bool only for specific adjustments (ie: remove sea level blue color in autoscale, etc) ie: autoscaling is done purely by the caller, who specifies the min/max for the gradient (and that should be enough)

void Gradient::getColor(const double& val, unsigned char& r, unsigned char& g, unsigned char& b, bool& a) const
{
	if(model==NULL) {
		throw UnknownValueException("Please set the color gradient before using it!", AT);
	}

	if(val==IOUtils::nodata) {
		r=0; g=0; b=0; a=true;
		return;
	}
	if(val==legend::bg_color) {
		r=channel_max_color-1; g=channel_max_color-1; b=channel_max_color-1; a=false;
		return;
	}
	if(val==legend::text_color) {
		r=0; g=0; b=0; a=false;
		return;
	}
	if(autoscale && delta==0) { //constant data throughout the grid & autoscale are no friends...
		r=g=b=channel_max_color/2; a=false;
		return;
	}

	a=false;
	double r_d,g_d,b_d;
	double val_norm;
	if(nr_unique_cols==0) {
		if(autoscale && val<min) val_norm=0.;
		else if(autoscale && val>max) val_norm=1.;
		else val_norm = (val-min)/delta;
	} else {
		if(autoscale && val<min) val_norm=0.;
		else if(autoscale && val>max) val_norm=1.;
		else val_norm = (static_cast<unsigned int>( (val-min)/delta*(double)nr_unique_cols )) / (double)nr_unique_cols;
	}
	model->getColor(val_norm, r_d, g_d, b_d);
	r = static_cast<unsigned char>(r_d*channel_max_color);
	g = static_cast<unsigned char>(g_d*channel_max_color);
	b = static_cast<unsigned char>(b_d*channel_max_color);
}

void Gradient::getColor(const double& val, unsigned char& index) const
{
	if(model==NULL) {
		throw UnknownValueException("Please set the color gradient before using it!", AT);
	}

	if(val==IOUtils::nodata) {
		index=0;
		return;
	}
	if(val==legend::bg_color) {
		index=1;
		return;
	}
	if(val==legend::text_color) {
		index=2;
		return;
	}
	if(delta==0) { //otherwise constant data throughout the grid makes a division by zero...
		index=nr_unique_cols/2 + reserved_idx;
		return;
	}

	if(nr_unique_cols==0) {
		throw UnknownValueException("Please define the number of colors for indexed gradients!", AT);
	}

	//watch out!! the palette contains some reserved values at the begining
	if(val<min) index=reserved_idx-2;
	else if(val>max) index=reserved_idx-1;
	else index = static_cast<unsigned char>( (val-min)/delta*(double)nr_unique_cols ) + reserved_idx;
}

void Gradient::getPalette(std::vector<unsigned char> &palette, size_t &nr_colors) const
{
	if(model==NULL) {
		throw UnknownValueException("Please set the color gradient before using it!", AT);
	}
	if(nr_unique_cols==0) {
		throw UnknownValueException("Please define the number of colors for indexed gradients!", AT);
	}

	palette.clear();

	//transparent color
	palette.push_back(channel_max_color); palette.push_back(channel_max_color); palette.push_back(channel_max_color);
	//legend background color
	palette.push_back(channel_max_color-1); palette.push_back(channel_max_color-1); palette.push_back(channel_max_color-1);
	//legend text color
	palette.push_back(0); palette.push_back(0); palette.push_back(0);

	double r_d, g_d, b_d;

	//underflow data color
	model->getColor(-0.1, r_d, g_d, b_d);
	palette.push_back( static_cast<unsigned char>(r_d*channel_max_color) );
	palette.push_back( static_cast<unsigned char>(g_d*channel_max_color) );
	palette.push_back( static_cast<unsigned char>(b_d*channel_max_color) );

	//overflow data color
	model->getColor(1.1, r_d, g_d, b_d);
	palette.push_back( static_cast<unsigned char>(r_d*channel_max_color) );
	palette.push_back( static_cast<unsigned char>(g_d*channel_max_color) );
	palette.push_back( static_cast<unsigned char>(b_d*channel_max_color) );

	//all normal colors
	for(unsigned char ii=0; ii<=nr_unique_cols; ii++) {
		const double val_norm = (double)ii/(double)nr_unique_cols;
		model->getColor(val_norm, r_d, g_d, b_d);
		palette.push_back( static_cast<unsigned char>(r_d*channel_max_color) );
		palette.push_back( static_cast<unsigned char>(g_d*channel_max_color) );
		palette.push_back( static_cast<unsigned char>(b_d*channel_max_color) );
	}

	const size_t nr_entries = palette.size();
	if((nr_entries%3) != 0) {
		std::stringstream ss;
		ss << "Error when creating color palette: " << nr_entries << " data points ";
		ss << "for 3 channels palette is impossible! (colors are interlaced)";
		throw IOException(ss.str(), AT);
	}
	nr_colors = nr_entries/3;
}

//we assume that the vectors are sorted by X
double Gradient_model::getInterpol(const double& val, const std::vector<double>& X, const std::vector<double>& Y) const
{
	const size_t nr = X.size();
#ifndef NOSAFECHECKS
	if(Y.size()!=nr) {
		std::stringstream ss;
		ss << "For color gradients interpolations, both X and Y vectors must have the same size! ";
		ss << "There are " << X.size() << " abscissa for " << Y.size() << " ordinates.";
		throw IOException(ss.str(), AT);
	}
	if(nr==0) {
		throw IOException("Empty vector of control points for color gradient interpolation", AT);
	}
#endif

	size_t i=0;
	while(i<nr && X[i]<val) i++; //find index of first element greater than val

	if(i==0) return Y[0];
	if(i>=nr) return Y[ nr-1 ];
	if(X[i]==val) return Y[i];

	const double y = Y[i-1] + (val-X[i-1])/(X[i]-X[i-1]) * (Y[i]-Y[i-1]);
	return y;
}

void Gradient_model::HSV2RGB(const double& h, const double& s, const double& v, unsigned char &r, unsigned char &g, unsigned char &b) const
{
	double r_d, g_d, b_d;
	Color::HSVtoRGB(h, s, v, r_d, g_d, b_d);
	r = static_cast<unsigned char>(r_d*255);
	g = static_cast<unsigned char>(g_d*255);
	b = static_cast<unsigned char>(b_d*255);
}

void Gradient_model::getColor(const double &val, double &r, double &g, double &b) const
{
	const double h = getInterpol(val, X, v_h);
	const double s = getInterpol(val, X, v_s);
	const double v = getInterpol(val, X, v_v);

	Color::HSVtoRGB(h, s, v, r, g, b);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// Various Gradients
/////////////////////////////////////////////////////////////////////////////////////////////////

void gr_heat::getColor(const double &i_val, double &r, double &g, double &b) const
{
	//normalize the input for under and over range
	double val=i_val;
	if(val<0.) val=0.;
	if(val>1.) val=1.;

	const double h = 240. * (1.-val);
	const double v = val*0.75+0.25;
	const double s = 1.-val*0.3;

	Color::HSVtoRGB(h, s, v, r, g, b);
}

gr_blue_pink::gr_blue_pink(const double& /*i_min*/, const double& /*i_max*/, const bool& /*i_autoscale*/) {
	X.push_back(0.); v_h.push_back(0.); v_s.push_back(0.); v_v.push_back(.95); //almost white
	X.push_back(0.2); v_h.push_back(172.); v_s.push_back(.4); v_v.push_back(.95); //light blue
	X.push_back(.4); v_h.push_back(213.); v_s.push_back(.4); v_v.push_back(.95); //violet
	X.push_back(.6); v_h.push_back(255.); v_s.push_back(.4); v_v.push_back(.95); //violet/blue
	X.push_back(.8); v_h.push_back(278.); v_s.push_back(.4); v_v.push_back(.95); //violet
	X.push_back(1.); v_h.push_back(359.); v_s.push_back(.3); v_v.push_back(.95); //red
}

gr_freeze::gr_freeze(const double& i_min, const double& i_max, const bool& /*i_autoscale*/) {
	//we want the yellow/green step to always be at zero celsius
	double begin=0., middle=0.5, end=1.;
	if(i_min<0. && i_max>0.) {
		const double range = i_max-i_min;
		middle = -i_min / range;
	} else if(i_min<=0. && i_max<=0.) {
		middle = 1.;
		end = 1000.;
	} else if(i_min>=0. && i_max>=0.) {
		begin = -1000.;
		middle = 0.;
	}
	X.push_back(begin); v_r.push_back(0.); v_g.push_back(0.); v_b.push_back(1.); //blue
	X.push_back(middle); v_r.push_back(1.); v_g.push_back(1.); v_b.push_back(0.); //yellow
	X.push_back(middle); v_r.push_back(0.); v_g.push_back(1.); v_b.push_back(0.); //green
	X.push_back(end); v_r.push_back(1.); v_g.push_back(0.); v_b.push_back(0.); //red
}

void gr_freeze::getColor(const double &val, double &r, double &g, double &b) const {
	//interpolation on RGB values
	r = getInterpol(val, X, v_r);
	g = getInterpol(val, X, v_g);
	b = getInterpol(val, X, v_b);
}

gr_blue::gr_blue(const double& /*i_min*/, const double& /*i_max*/, const bool& /*i_autoscale*/) {
	X.push_back(0.); v_h.push_back(0.); v_s.push_back(0.); v_v.push_back(.99); //5
	X.push_back(.16667); v_h.push_back(180.); v_s.push_back(.2); v_v.push_back(.99); //10
	X.push_back(.33334); v_h.push_back(193.); v_s.push_back(.32); v_v.push_back(.97); //20
	X.push_back(.5); v_h.push_back(205.); v_s.push_back(.43); v_v.push_back(.94); //50
	X.push_back(.66667); v_h.push_back(219.); v_s.push_back(.55); v_v.push_back(.91); //80
	X.push_back(.83335); v_h.push_back(231.); v_s.push_back(.66); v_v.push_back(.88); //120
	X.push_back(1.); v_h.push_back(244.); v_s.push_back(.78); v_v.push_back(.85); //200
	X.push_back(1.); v_h.push_back(270.); v_s.push_back(1.); v_v.push_back(.8); //200
}

gr_bg_isomorphic::gr_bg_isomorphic(const double& /*i_min*/, const double& /*i_max*/, const bool& /*i_autoscale*/) {
	X.push_back(0.); v_h.push_back(47.); v_s.push_back(.92); v_v.push_back(0.); //black
	X.push_back(.5); v_h.push_back(178.); v_s.push_back(.58); v_v.push_back(.67); //light blue
	X.push_back(1.); v_h.push_back(84.); v_s.push_back(.67); v_v.push_back(1.); //light green
}

gr_pastel::gr_pastel(const double& /*i_min*/, const double& /*i_max*/, const bool& /*i_autoscale*/) {
	X.push_back(0.); v_h.push_back(0.); v_s.push_back(1.); v_v.push_back(0.); //black
	X.push_back(1./7.); v_h.push_back(185.); v_s.push_back(.26); v_v.push_back(.56); //light blue
	X.push_back(2./7.); v_h.push_back(122.); v_s.push_back(.44); v_v.push_back(.91); //light green
	X.push_back(4./7.); v_h.push_back(60.); v_s.push_back(.44); v_v.push_back(.91); //light yellow
	X.push_back(5./7.); v_h.push_back(22.); v_s.push_back(.44); v_v.push_back(.91); //orange
	X.push_back(6./7.); v_h.push_back(0.); v_s.push_back(.44); v_v.push_back(.91); //red
	X.push_back(7./7.); v_h.push_back(0.); v_s.push_back(.4); v_v.push_back(.7); //dark red
}

gr_terrain::gr_terrain(const double& /*i_min*/, const double& i_max, const bool& i_autoscale) {
	if(i_autoscale) {
		X.push_back(0.); v_h.push_back(144.); v_s.push_back(.50); v_v.push_back(.39); //sea level, dark green
		X.push_back(.25); v_h.push_back(46.); v_s.push_back(.54); v_v.push_back(.86); //yellow
		X.push_back(.5); v_h.push_back(4.); v_s.push_back(.71); v_v.push_back(.53); //dark red
		X.push_back(.75); v_h.push_back(22.); v_s.push_back(.88); v_v.push_back(.41); //maroon
		X.push_back(1.); v_h.push_back(22.); v_s.push_back(.2); v_v.push_back(.5); //light maroon
	} else {
		const double snow_line = i_max - 500.;
		X.push_back(0.); v_h.push_back(198.); v_s.push_back(.50); v_v.push_back(.74); //sea, light blue
		X.push_back(0.); v_h.push_back(144.); v_s.push_back(.50); v_v.push_back(.39); //sea level, dark green
		X.push_back(.4); v_h.push_back(46.); v_s.push_back(.54); v_v.push_back(.86); //yellow, 1200m
		X.push_back(.73333); v_h.push_back(4.); v_s.push_back(.71); v_v.push_back(.53); //dark red, 2200m
		X.push_back(.9); v_h.push_back(22.); v_s.push_back(.88); v_v.push_back(.41); //maroon, 2700m
		X.push_back(.98333); v_h.push_back(22.); v_s.push_back(.36); v_v.push_back(.79); //light maroon, 2950m
		X.push_back(1.); v_h.push_back(0.); v_s.push_back(0.); v_v.push_back(.7); //light gray == permanent snow line, 3000m
		for(size_t i=0; i<X.size(); i++) X[i] = X[i]*snow_line/i_max; //snow line is the reference

		X.push_back(1.); v_h.push_back(0.); v_s.push_back(0.); v_v.push_back(.95); //almost white == fully glaciated line
	}
}

gr_slope::gr_slope(const double& /*i_min*/, const double& /*i_max*/, const bool& /*i_autoscale*/) {
	//usually, between 0 and 50
	X.push_back(0.); v_h.push_back(185.); v_s.push_back(.26); v_v.push_back(.56); //light blue
	X.push_back(.5); v_h.push_back(122.); v_s.push_back(.44); v_v.push_back(.91); //light green, 25
	X.push_back(.6); v_h.push_back(60.); v_s.push_back(.44); v_v.push_back(.91); //light yellow, 30
	X.push_back(.7); v_h.push_back(22.); v_s.push_back(.44); v_v.push_back(.91); //orange, 35
	X.push_back(.8); v_h.push_back(0.); v_s.push_back(.44); v_v.push_back(.91); //red, 40
	X.push_back(.9); v_h.push_back(0.); v_s.push_back(.58); v_v.push_back(.35); //dark red, 45
	X.push_back(1.); v_h.push_back(0.); v_s.push_back(1.); v_v.push_back(0.); //black, 50
}

gr_azi::gr_azi(const double& /*i_min*/, const double& /*i_max*/, const bool& i_autoscale) {
	if(i_autoscale) {
		X.push_back(0.); v_h.push_back(113.); v_s.push_back(.66); v_v.push_back(.91); //light green
		X.push_back(1.); v_h.push_back(360.); v_s.push_back(.66); v_v.push_back(.91); //light red
	} else {
		X.push_back(0.); v_h.push_back(240.); v_s.push_back(.78); v_v.push_back(1.); //blue
		X.push_back(.25); v_h.push_back(310.); v_s.push_back(.78); v_v.push_back(1.); //magenta
		X.push_back(.5); v_h.push_back(360.); v_s.push_back(1.); v_v.push_back(1.); //red, increasing hue
		X.push_back(.5); v_h.push_back(0.); v_s.push_back(1.); v_v.push_back(1.); //red, back to hue=0
		X.push_back(.75); v_h.push_back(28.); v_s.push_back(.78); v_v.push_back(1.); //orange
		X.push_back(1.); v_h.push_back(240.); v_s.push_back(.78); v_v.push_back(1.); //back to blue
	}
}

} //namespace
