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
#include "PNGIO.h"
#include <meteoio/ResamplingAlgorithms2D.h>

#include <algorithm>

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

const unsigned int legend::font_0[10][6] = {{0,0,1,1,0,0}, {0,1,0,0,1,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,0,1,1,0,1}, {1,0,0,0,0,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,1,0,0,1,0}, {0,0,1,1,0,0}};
const unsigned int legend::font_1[10][6] = {{0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,1,1,1,0,0}, {1,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,0,1,1,0,0}, {0,1,1,1,1,0}};
const unsigned int legend::font_2[10][6] = {{0,0,1,1,1,0}, {0,1,0,0,0,1}, {0,0,0,0,0,1}, {0,0,0,0,0,1}, {0,0,0,0,1,1}, {0,0,0,1,1,0}, {0,0,1,1,0,0}, {0,1,1,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,1,1}};
const unsigned int legend::font_3[10][6] = {{0,0,1,1,1,0}, {0,1,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,1,1,1,0}, {0,0,0,0,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {0,1,1,1,1,0}};
const unsigned int legend::font_4[10][6] = {{0,0,0,1,1,0}, {0,0,1,1,1,0}, {0,1,1,1,1,0}, {1,1,0,1,1,0}, {1,0,0,1,1,0}, {1,1,1,1,1,1}, {0,0,0,1,1,0}, {0,0,0,1,1,0}, {0,0,0,1,1,0}, {0,0,1,1,1,1}};
const unsigned int legend::font_5[10][6] = {{1,1,1,1,1,1}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,0,0}, {0,0,0,0,1,0}, {0,0,0,0,1,1}, {0,0,0,0,1,1}, {1,0,0,0,1,0}, {0,1,1,1,0,0}};
const unsigned int legend::font_6[10][6] = {{0,1,1,1,1,0}, {1,1,0,0,1,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,1,1,0,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,0,1,1,0,0}};
const unsigned int legend::font_7[10][6] = {{1,1,1,1,1,1}, {0,0,0,0,1,1}, {0,0,0,1,1,0}, {0,0,0,1,0,0}, {0,0,1,1,0,0}, {0,0,1,0,0,0}, {0,1,1,0,0,0}, {0,1,0,0,0,0}, {1,1,0,0,0,0}, {1,1,0,0,0,0}};
const unsigned int legend::font_8[10][6] = {{0,0,1,1,0,0}, {0,1,0,0,1,0}, {0,1,0,0,1,0}, {0,1,0,0,1,0}, {0,0,1,1,0,0}, {0,1,0,0,1,0}, {1,1,0,0,1,1}, {1,1,0,0,1,1}, {0,1,0,0,1,0}, {0,0,1,1,0,0}};
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

/**
 * @page template PNGIO
 * @section template_format Format
 * *Put here the informations about the standard format that is implemented*
 *
 * @section template_units Units
 *
 *
 * @section template_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - etc
 */

const double PNGIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const double PNGIO::factor = 2.0; //image scale factor
const bool PNGIO::autoscale = true;
const bool PNGIO::has_legend = true;

PNGIO::PNGIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PNGIO::PNGIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PNGIO::PNGIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PNGIO::~PNGIO() throw()
{

}

void PNGIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                             std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                             const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                              const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::write2DGrid(const Grid2DObject& grid_in, const std::string& filename)
{
	FILE *fp;
	png_structp png_ptr=NULL;
	png_infop info_ptr=NULL;
	png_bytep row=NULL;

	//scale input image
	Grid2DObject grid = ResamplingAlgorithms2D::BilinearResampling(grid_in, factor);
	const double ncols = grid.ncols, nrows = grid.nrows;
	const double min = grid.grid2D.getMin();
	const double max = grid.grid2D.getMax();

	// Open file for writing (binary mode)
	if (!IOUtils::validFileName(filename)) {
		throw InvalidFileNameException(filename, AT);
	}
	fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		throw FileAccessException(filename, AT);
	}

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fclose(fp);
		throw IOException("Could not allocate write structure", AT);
	}

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fclose(fp);
		png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		throw IOException("Could not allocate info structure", AT);
	}

	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		cleanup(fp, png_ptr, info_ptr, row);
		throw IOException("Error during png creation", AT);
	}

	png_init_io(png_ptr, fp);

	unsigned int full_width=ncols;
	Array2D<double> legend_array;
	if(has_legend) {
		legend leg(nrows, min, max);
		legend_array = leg.getLegend();
		unsigned int nx, ny;
		legend_array.size(nx,ny);
		full_width += nx;
	}
	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, full_width, nrows,
	             8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	writeMetadata(png_ptr, info_ptr);
	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (4 bytes per pixel - RGBA)
	row = (png_bytep) malloc(4 * full_width * sizeof(png_byte));

	// Write image data
	for(int y=nrows-1 ; y>=0 ; y--) {
		for(unsigned int x=0 ; x<ncols ; x++) {
			setRGB( grid(x,y), min, max, &(row[x*4]) );
		}
		for(unsigned int x=ncols; x<full_width; x++) {
			setRGB( legend_array(x-ncols,y), min, max, &(row[x*4]) );
		}
		png_write_row(png_ptr, row);
	}

	png_write_end(png_ptr, NULL);
	cleanup(fp, png_ptr, info_ptr, row);
}

void PNGIO::cleanup() throw()
{

}

//val between 0 and 1
//return values between 0 and 255 per channel
void PNGIO::setRGB(double val, const double& min, const double& max, png_byte *ptr)
{
	if(val==IOUtils::nodata) {
		ptr[0]=0; ptr[1]=0; ptr[2]=0; ptr[3]=0;
		return;
	}
	if(val==legend::bg_color) {
		ptr[0]=255; ptr[1]=255; ptr[2]=255; ptr[3]=255;
		return;
	}
	if(val==legend::text_color) {
		ptr[0]=0; ptr[1]=0; ptr[2]=0; ptr[3]=255;
		return;
	}

	ptr[3]=255; //no alpha for valid values

	if(autoscale) val = (val-min)/(max-min);

	const double h = 240. * (1.-val);
	const double v = val*0.75+0.25;
	double r,g,b;
	HSVtoRGB(h, 1., v, r, g, b);
	ptr[0] = static_cast<int>(r*255);
	ptr[1] = static_cast<int>(g*255);
	ptr[2] = static_cast<int>(b*255);
}

//values between 0 and 1
//see http://www.cs.rit.edu/~ncs/color/t_convert.html or https://secure.wikimedia.org/wikipedia/en/wiki/HSL_and_HSV#Conversion_from_HSL_to_RGB
void PNGIO::RGBtoHSV(const double r, const double g, const double b,
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
void PNGIO::HSVtoRGB(const double h, const double s, const double v, double &r, double &g, double &b)
{
	if( s==0 ) {
		// achromatic (grey)
		r = g = b = v;
		return;
	}

	const double h_p = h/60;	// sector 0 to 5
	const int i = (int)floor(h_p);
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

void PNGIO::writeMetadata(png_structp &png_ptr, png_infop &info_ptr)
{
	const std::string version = "MeteoIO "+getLibVersion();
	char version_c[79]=""; strncpy(version_c, version.c_str(), 79);
	const std::string logname = IOUtils::getLogName();
	char logname_c[79]=""; strncpy(logname_c, logname.c_str(), 79);
	png_text info_text[3];
	info_text[0].key = (char*)"Title";
	info_text[0].text = (char*)"Gridded data";
	info_text[0].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[1].key = (char*)"Author";
	info_text[1].text = logname_c;
	info_text[1].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[2].key = (char*)"Software";
	info_text[2].text = version_c;
	info_text[2].compression = PNG_TEXT_COMPRESSION_NONE;
	//TODO: add geolocalization tags: latitude, longitude, altitude, cellsize + indicator of position: llcorner
	//add data set timestamp
	png_set_text(png_ptr, info_ptr, info_text, 3);
}

void PNGIO::cleanup(FILE *fp, png_structp png_ptr, png_infop info_ptr, png_bytep row)
{
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) free(row);
}

#ifndef _METEOIO_JNI
extern "C"
{
#define COMPILE_PLUGIN
#include "exports.h"

	METEOIO_EXPORT void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	METEOIO_EXPORT void* loadObject(const string& classname, const Config& cfg) {
		if(classname == "PNGIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new PNGIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
