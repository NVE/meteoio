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
#include <meteoio/Graphics.h>

#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page template PNGIO
 * @section template_format Format
 * *Put here the informations about the standard format that is implemented*
 * Finally, the naming scheme for meteo grids should be: YYYYMMDDHHmm_{MeteoGrids::Parameters}.png
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

void PNGIO::read2DGrid(Grid2DObject&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PNGIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters& , const Date&)
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
	writeMetadata(grid, png_ptr, info_ptr);
	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (4 bytes per pixel - RGBA)
	row = (png_bytep) malloc(4 * full_width * sizeof(png_byte));

	// Write image data
	Gradient gradient(Gradient::heat, min, max);
	for(int y=nrows-1 ; y>=0 ; y--) {
		for(unsigned int x=0 ; x<ncols ; x++) {
			const unsigned int i=x*4;
			unsigned char r,g,b,a;
			gradient.getColor(grid(x,y), r,g,b,a);
			row[i]=r; row[i+1]=g; row[i+2]=b; row[i+3]=a;
		}
		for(unsigned int x=ncols; x<full_width; x++) {
			//setRGB( legend_array(x-ncols,y), min, max, &(row[x*4]) );
			const unsigned int i=x*4;
			unsigned char r,g,b,a;
			gradient.getColor(legend_array(x-ncols,y), r,g,b,a);
			row[i]=r; row[i+1]=g; row[i+2]=b; row[i+3]=a;
		}
		png_write_row(png_ptr, row);
	}

	png_write_end(png_ptr, NULL);
	cleanup(fp, png_ptr, info_ptr, row);
}

void PNGIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	std::stringstream ss;
	ss << date.toString(Date::NUM) << "_" << MeteoGrids::getParameterName(parameter) << ".png";
	write2DGrid(grid_in, ss.str());
}

void PNGIO::cleanup() throw()
{

}

void PNGIO::writeMetadata(const Grid2DObject& grid, png_structp &png_ptr, png_infop &info_ptr)
{
	const std::string version = "MeteoIO "+getLibVersion();
	char version_c[79]=""; strncpy(version_c, version.c_str(), 79);
	const std::string logname = IOUtils::getLogName();
	char logname_c[79]=""; strncpy(logname_c, logname.c_str(), 79);
	char latitude[12]=""; snprintf(latitude, 12, "%.6f", grid.llcorner.getLat());
	char longitude[12]=""; snprintf(longitude, 12, "%.6f", grid.llcorner.getLon());
	char cellsize[12]=""; snprintf(cellsize, 12, "%.2f", (grid.cellsize*factor)); //HACK: not working...
	
	png_text info_text[7];
	info_text[0].key = (char*)"Title";
	info_text[0].text = (char*)"Gridded data"; //HACK: write meteogrid parameter
	info_text[0].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[1].key = (char*)"Author";
	info_text[1].text = logname_c;
	info_text[1].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[2].key = (char*)"Software";
	info_text[2].text = version_c;
	info_text[2].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[3].key = (char*)"Position"; //HACK: see exif key
	info_text[3].text = (char*)"llcorner";
	info_text[3].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[4].key = (char*)"Latitude"; //HACK: see exif key
	info_text[4].text = latitude;
	info_text[4].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[5].key = (char*)"Longitude";
	info_text[5].text = longitude;
	info_text[5].compression = PNG_TEXT_COMPRESSION_NONE;
	info_text[6].key = (char*)"Cellsize";
	info_text[6].text = cellsize;
	info_text[6].compression = PNG_TEXT_COMPRESSION_NONE;
	//TODO: add geolocalization tags: latitude, longitude, altitude, cellsize + indicator of position: llcorner
	//add data set timestamp
	png_set_text(png_ptr, info_ptr, info_text, 7);
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
