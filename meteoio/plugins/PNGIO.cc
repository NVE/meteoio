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

	if(autoscale) {
		const double min = grid.grid2D.getMin();
		const double max = grid.grid2D.getMax();
		grid.grid2D -= min;
		grid.grid2D /= (max-min);
	}

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

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, ncols, nrows,
	             8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	writeMetadata(png_ptr, info_ptr);
	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (4 bytes per pixel - RGBA)
	row = (png_bytep) malloc(4 * ncols * sizeof(png_byte));

	// Write image data
	for (int y=nrows-1 ; y>=0 ; y--) {
		for (unsigned int x=0 ; x<ncols ; x++) {
			setRGB2( grid(x,y), &(row[x*4]) );
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
void PNGIO::setRGB(const double val, png_byte *ptr)
{
	if(val==IOUtils::nodata) {
		ptr[0]=0; ptr[1]=0; ptr[2]=0; ptr[3]=0;
		return;
	}

	ptr[3]=255; //no alpha for valid values
	int v = (int)(val * 768);
	if (v < 0) v = 0;
	if (v > 768) v = 768;
	int offset = v % 256;

	if (v<256) {
		ptr[0]=0; ptr[1]=0; ptr[2]=offset;
	}
	else if (v<512) {
		ptr[0]=0; ptr[1]=offset; ptr[2]=255-offset;
	}
	else {
		ptr[0]=offset; ptr[1]=255-offset; ptr[2]=0;
	}
}

//val between 0 and 1
//return values between 0 and 255 per channel
void PNGIO::setRGB2(const double val, png_byte *ptr)
{
	if(val==IOUtils::nodata) {
		ptr[0]=0; ptr[1]=0; ptr[2]=0; ptr[3]=0;
		return;
	}

	ptr[3]=255; //no alpha for valid values

	const double h = 240. * (1.-val);
	const double v = val*0.75+0.25;
	double r,g,b;
	HSVtoRGB(h, 1., v, r, g, b);
	ptr[0] = static_cast<int>(r*255);
	ptr[1] = static_cast<int>(g*255);
	ptr[2] = static_cast<int>(b*255);
}

//values between 0 and 1
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
	//TODO: add geolocalization
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
