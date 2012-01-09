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
 * @page pngio PNGIO
 * @section template_format Format
 * This plugin write data to the Portable Network Graphics format (see https://secure.wikimedia.org/wikipedia/en/wiki/Portable_Network_Graphics).
 * No data read has been implemented, because reading an existing file would require the exact knowlege of the color gradient that has been used
 * to create it. When writing grids, various color gradients will be used depending on the parameter that the data represents. Nodata values
 * are represented by transparent pixels (transparency is acheived through a transparent color instead of a true alpha channel for size and performance).
 * Finally, the naming scheme for meteo grids should be: YYYYMMDDHHmm_{MeteoGrids::Parameters}.png
 *
 * @section template_units Units
 * All units are MKSA except temperatures that are expressed in celcius.
 *
 * @section template_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Output] section
 * - GRID2DPATH: meteo grids directory where to read the grids; [Output] section
 * - PNG_LEGEND: plot legend on the side of the graph? (default: true)
 * - PNG_MIN_SIZE: guarantee that a 2D plot will have at least the given size
 * - PNG_MAX_SIZE: guarantee that a 2D plot will have at most the given size
 * - PNG_SCALING: scaling algorithm, either nearest or bilinear (default=bilinear)
 * - PNG_AUTOSCALE: autoscale for the color gradient? (default=true)
 * - PNG_WORLD_FILE: create world file with each file? (default=false)
 *
 * The size are specified as width followed by height, with the separator being either a space, 'x' or '*'. If a minimum and a maximum size are given, the average of the smallest and largest permissible sizes will be used.
 * The world file is used for geolocalization and goes alongside the graphics output. By convention,
 * the file has the same name as the image file, with the third letter of the extension jammed with a w: tif->tfw, jpg->jqw.
 * The format is the following:
 * @code
 *    5.000000000000 (size of pixel in x direction)
 *    0.000000000000 (rotation term for row)
 *    0.000000000000 (rotation term for column)
 *    -5.000000000000 (size of pixel in y direction)
 *    492169.690845528910 (x coordinate of centre of upper left pixel in map units)
 *    5426523.318065105000 (y coordinate of centre of upper left pixel in map units)
 * @endcode
 *
 * @section example Example use
 * @code
 * GRID2D = PNG
 * png_legend = false
 * png_min_size = 400x400
 * png_max_size = 1366*768
 * @endcode
 */

const double PNGIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)
const unsigned char PNGIO::channel_depth = 8;
const unsigned char PNGIO::channel_max_color = 255;
const unsigned char PNGIO::transparent_grey = channel_max_color;

PNGIO::PNGIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	setOptions();
}

PNGIO::PNGIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	setOptions();
}

PNGIO::PNGIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	setOptions();
}

void PNGIO::setOptions()
{
	cfg.getValue("COORDSYS", "Output", coordout);
	cfg.getValue("COORDPARAM", "Output", coordoutparam, Config::nothrow);
	cfg.getValue("GRID2DPATH", "Output", grid2dpath);
	//cfg.getValue("TIME_ZONE", "Output", tz_out, Config::nothrow);

	//get size specifications
	std::string min_size, max_size;
	min_w = min_h = max_w = max_h = IOUtils::unodata;
	cfg.getValue("PNG_MIN_SIZE", "Output", min_size, Config::nothrow);
	if(min_size!="") parse_size(min_size, min_w, min_h);
	cfg.getValue("PNG_MAX_SIZE", "Output", max_size, Config::nothrow);
	if(max_size!="") parse_size(max_size, max_w, max_h);

	autoscale = true;
	cfg.getValue("PNG_AUTOSCALE", "Output", autoscale, Config::nothrow);
	has_legend = true;
	cfg.getValue("PNG_LEGEND", "Output", has_legend, Config::nothrow);
	scaling = "bilinear";
	cfg.getValue("PNG_SCALING", "Output", scaling, Config::nothrow);
	has_world_file=false;
	cfg.getValue("PNG_WORLD_FILE", "Output", has_world_file, Config::nothrow);

	if(has_legend) { //we need to save room for the legend
		if(min_w!=IOUtils::unodata) min_w -= legend::getLegendWidth();
		if(max_w!=IOUtils::unodata) max_w -= legend::getLegendWidth();
	}
}

void PNGIO::parse_size(const std::string& size_spec, unsigned int& width, unsigned int& height)
{
	char rest[32] = "";
	if(sscanf(size_spec.c_str(), "%u %u%31s", &width, &height, rest) < 2)
	if(sscanf(size_spec.c_str(), "%u*%u%31s", &width, &height, rest) < 2)
	if(sscanf(size_spec.c_str(), "%ux%u%31s", &width, &height, rest) < 2) {
		std::stringstream ss;
		ss << "Can not parse PNGIO size specification \"" << size_spec << "\"";
		throw InvalidFormatException(ss.str(), AT);
	}
	std::string tmp(rest);
	IOUtils::trim(tmp);
	if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
		std::stringstream ss;
		ss << "Invalid PNGIO size specification \"" << size_spec << "\"";
		throw InvalidFormatException(ss.str(), AT);
	}
}

double PNGIO::getScaleFactor(const double& grid_w, const double& grid_h)
{
	if(grid_w==0 || grid_h==0) {
		return 1.;
	}

	double min_factor = IOUtils::nodata;
	if(min_w!=IOUtils::unodata) { //min_w & min_w are read together
		const double min_w_factor = (double)min_w / (double)grid_w;
		const double min_h_factor = (double)min_h / (double)grid_h;
		min_factor = std::max(min_w_factor, min_h_factor);
	}

	double max_factor = IOUtils::nodata;
	if(max_w!=IOUtils::unodata) { //max_w & max_h are read together
		const double max_w_factor = (double)max_w / (double)grid_w;
		const double max_h_factor = (double)max_h / (double)grid_h;
		max_factor = std::min(max_w_factor, max_h_factor);
	}

	if(min_factor==IOUtils::nodata && max_factor==IOUtils::nodata)
		return 1.; //no user given specification
	if(min_factor!=IOUtils::nodata && max_factor!=IOUtils::nodata)
		return (min_factor+max_factor)/2.; //both min & max -> average

	//only one size specification provided -> return its matching factor
	if(min_factor!=IOUtils::nodata)
		return min_factor;
	else
		return max_factor;
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

Grid2DObject PNGIO::scaleGrid(const Grid2DObject& grid_in)
{ //scale input image
	const double factor = getScaleFactor(grid_in.ncols, grid_in.nrows);
	if(scaling=="nearest")
		return ResamplingAlgorithms2D::NearestNeighbour(grid_in, factor);
	else if(scaling=="bilinear")
		return ResamplingAlgorithms2D::BilinearResampling(grid_in, factor);
	else {
		stringstream ss;
		ss << "Grid scaling algorithm \"" << scaling << "\" unknown";
		throw UnknownValueException(ss.str(), AT);
	}
}

void PNGIO::setFile(const std::string& filename, FILE *fp, png_structp& png_ptr, png_infop& info_ptr, const unsigned int &width, const unsigned int &height)
{
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
		cleanup(fp, png_ptr, info_ptr);
		throw IOException("Error during png creation", AT);
	}

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth). Alpha channel with PNG_COLOR_TYPE_RGB_ALPHA
	png_set_IHDR(png_ptr, info_ptr, width, height,
	             channel_depth, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	//set transparent color (ie: cheap transparency: leads to smaller files and shorter run times)
	png_color_16 trans_rgb_value = {transparent_grey, transparent_grey, transparent_grey, transparent_grey, transparent_grey};
	png_set_tRNS(png_ptr, info_ptr, 0, 0, &trans_rgb_value);
}

unsigned int PNGIO::setLegend(const unsigned int &ncols, const unsigned int &nrows, const double &min, const double &max, Array2D<double> &legend_array)
{
	if(has_legend) {
		legend leg(nrows, min, max);
		legend_array = leg.getLegend();
		unsigned int nx, ny;
		legend_array.size(nx,ny);
		return (ncols+nx);
	} else {
		return ncols;
	}
}

void PNGIO::writeDataSection(const Grid2DObject &grid, const Array2D<double> &legend_array, const Gradient &gradient, const unsigned int &full_width, png_structp &png_ptr)
{
	const double ncols = grid.ncols, nrows = grid.nrows;

	// Allocate memory for one row (3 bytes per pixel - RGB)
	const unsigned char channels = 3;
	png_bytep row=NULL;
	row = (png_bytep) malloc(channels * full_width * sizeof(png_byte));

	// Write image data
	for(int y=nrows-1 ; y>=0 ; y--) {
		unsigned int x=0;
		for(; x<ncols ; x++) {
			const unsigned int i=x*channels;
			unsigned char r,g,b;
			bool a;
			gradient.getColor(grid(x,y), r,g,b,a);
			if(a==true) {
				row[i]=transparent_grey; row[i+1]=transparent_grey; row[i+2]=transparent_grey;
			} else {
				row[i]=r; row[i+1]=g; row[i+2]=b;
			}
		}
		for(; x<full_width; x++) {
			const unsigned int i=x*channels;
			unsigned char r,g,b;
			bool a;
			gradient.getColor(legend_array(x-ncols,y), r,g,b,a);
			if(a==true) {
				row[i]=transparent_grey; row[i+1]=transparent_grey; row[i+2]=transparent_grey;
			} else {
				row[i]=r; row[i+1]=g; row[i+2]=b;
			}
		}
		png_write_row(png_ptr, row);
	}

	free(row);
}

void PNGIO::writeWorldFile(const Grid2DObject& grid_in, const std::string& filename)
{
	const string world_file = IOUtils::removeExtension(filename)+".pnw";
	const double cellsize = grid_in.cellsize;
	Coords world_ref=grid_in.llcorner;
	world_ref.setProj(coordout, coordoutparam);
	world_ref.moveByXY(0., grid_in.nrows*cellsize);

	std::ofstream fout;
	fout.open(world_file.c_str());
	if (fout.fail()) {
		throw FileAccessException(world_file, AT);
	}

	try {
		fout << std::setprecision(12) << cellsize << "\n";
		fout << "0.000000000000\n";
		fout << "0.000000000000\n";
		fout << std::setprecision(12) << -cellsize << "\n";
		fout << std::setprecision(12) << world_ref.getEasting() << "\n";
		fout << std::setprecision(12) << world_ref.getNorthing() << "\n";
	} catch(...) {
		fout.close();
		throw FileAccessException("Failed when writing to PNG world file \""+world_file+"\"", AT);
	}

	fout.close();
}

void PNGIO::write2DGrid(const Grid2DObject& grid_in, const std::string& filename)
{
	string full_name = grid2dpath+"/"+filename;
	FILE *fp=NULL;
	png_structp png_ptr=NULL;
	png_infop info_ptr=NULL;

	//scale input image
	const Grid2DObject grid = scaleGrid(grid_in);
	const double ncols = grid.ncols, nrows = grid.nrows;
	const double min = grid.grid2D.getMin();
	const double max = grid.grid2D.getMax();

	Gradient gradient(Gradient::heat, min, max, autoscale);
	gradient.setNrOfLevels(50);

	Array2D<double> legend_array; //it will remain empty if there is no legend
	const unsigned int full_width = setLegend(ncols, nrows, min, max, legend_array);

	setFile(full_name, fp, png_ptr, info_ptr, full_width, nrows);
	if(has_world_file) writeWorldFile(grid, full_name);

	createMetadata(grid);
	metadata_key.push_back("Title"); //adding generic title
	metadata_text.push_back("Unknown Gridded data");
	writeMetadata(png_ptr, info_ptr);

	writeDataSection(grid, legend_array, gradient, full_width, png_ptr);

	png_write_end(png_ptr, NULL);
	cleanup(fp, png_ptr, info_ptr);
}

void PNGIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	std::string filename;
	if(parameter==MeteoGrids::DEM || parameter==MeteoGrids::SLOPE || parameter==MeteoGrids::AZI)
		filename = grid2dpath + "/" + MeteoGrids::getParameterName(parameter) + ".png";
	else
		filename = grid2dpath + "/" + date.toString(Date::NUM) + "_" + MeteoGrids::getParameterName(parameter) + ".png";

	FILE *fp=NULL;
	png_structp png_ptr=NULL;
	png_infop info_ptr=NULL;

	//scale input image
	Grid2DObject grid = scaleGrid(grid_in);
	const double ncols = grid.ncols, nrows = grid.nrows;
	double min = grid.grid2D.getMin();
	double max = grid.grid2D.getMax();

	Gradient gradient;
	if(parameter==MeteoGrids::DEM) {
		if(!autoscale) {
			min = 0.; //we want a 3000 snow line with a full scale legend
			max = 3500.;
			gradient.set(Gradient::terrain, min, max, autoscale); //max used as snow line reference
		} else
			gradient.set(Gradient::terrain, min, max, autoscale);
	} else if(parameter==MeteoGrids::SLOPE) {
		gradient.set(Gradient::slope, min, max, autoscale);
	} else if(parameter==MeteoGrids::AZI) {
		if(!autoscale) {
			min = 0.;
			max = 360.;
		}
		gradient.set(Gradient::azi, min, max, autoscale);
	} else if(parameter==MeteoGrids::HS) {
		if(!autoscale) {
			min = 0.; max = 3.5;
		}
		gradient.set(Gradient::blue, min, max, autoscale);
	} else if(parameter==MeteoGrids::TA) {
		grid.grid2D -= Cst::t_water_freezing_pt; //convert to celsius
		if(!autoscale) {
			min = -10.; max = 10.;
		} else {
			min -= Cst::t_water_freezing_pt;
			max -= Cst::t_water_freezing_pt;
		}
		gradient.set(Gradient::heat, min, max, autoscale);
	} else if(parameter==MeteoGrids::TSS) {
		grid.grid2D -= Cst::t_water_freezing_pt; //convert to celsius
		if(!autoscale) {
			min = -10.; max = 10.;
		} else {
			min -= Cst::t_water_freezing_pt;
			max -= Cst::t_water_freezing_pt;
		}
		gradient.set(Gradient::freeze, min, max, autoscale);
	} else if(parameter==MeteoGrids::RH) {
		if(!autoscale) {
			min = 0.; max = 1.;
		}
		gradient.set(Gradient::bg_isomorphic, min, max, autoscale);
	} else if(parameter==MeteoGrids::ALB) {
		if(!autoscale) {
			min = 0.; max = 1.;
		}
		gradient.set(Gradient::bg_isomorphic, min, max, autoscale);
	} else if(parameter==MeteoGrids::SWE) {
		if(!autoscale) {
			min = 0.; max = 2000.;
		}
		gradient.set(Gradient::blue, min, max, autoscale);
	} else {
		gradient.set(Gradient::heat, min, max, autoscale);
	}
	gradient.setNrOfLevels(30);

	Array2D<double> legend_array; //it will remain empty if there is no legend
	const unsigned int full_width = setLegend(ncols, nrows, min, max, legend_array);

	setFile(filename, fp, png_ptr, info_ptr, full_width, nrows);
	if(has_world_file) writeWorldFile(grid, filename);

	createMetadata(grid);
	metadata_key.push_back("Title"); //adding title
	metadata_text.push_back( MeteoGrids::getParameterName(parameter)+" on "+date.toString(Date::ISO) );
	metadata_key.push_back("Simulation Date");
	metadata_text.push_back( date.toString(Date::ISO) );
	metadata_key.push_back("Simulation Parameter");
	metadata_text.push_back( MeteoGrids::getParameterName(parameter) );
	writeMetadata(png_ptr, info_ptr);

	writeDataSection(grid, legend_array, gradient, full_width, png_ptr);

	png_write_end(png_ptr, NULL);
	cleanup(fp, png_ptr, info_ptr);
}

void PNGIO::createMetadata(const Grid2DObject& grid)
{
	const double lat = grid.llcorner.getLat();
	const double lon = grid.llcorner.getLon();
	stringstream ss;

	metadata_key.push_back("Creation Time");
	Date cr_date;
	cr_date.setFromSys();
	metadata_text.push_back( cr_date.toString(Date::ISO) );
	metadata_key.push_back("Author");
	metadata_text.push_back(IOUtils::getLogName());
	metadata_key.push_back("Software");
	metadata_text.push_back("MeteoIO "+getLibVersion());
	metadata_key.push_back("Position");
	metadata_text.push_back("llcorner");
	metadata_key.push_back("Cellsize");
	ss.str(""); ss << fixed << setprecision(2) << grid.cellsize;
	metadata_text.push_back(ss.str());
	metadata_key.push_back("Latitude");
	ss.str(""); ss << fixed << setprecision(6) << lat;
	metadata_text.push_back(ss.str());
	metadata_key.push_back("Longitude");
	ss.str(""); ss << fixed << setprecision(6) << lon;
	metadata_text.push_back(ss.str());

	if(lat<0.) {
		metadata_key.push_back("LatitudeRef");
		metadata_text.push_back("S");
		metadata_key.push_back("GPSLatitude");
		metadata_text.push_back(decimal_to_dms(-lat));
	} else {
		metadata_key.push_back("LatitudeRef");
		metadata_text.push_back("N");
		metadata_key.push_back("GPSLatitude");
		metadata_text.push_back(decimal_to_dms(lat));
	}
	if(lon<0.) {
		metadata_key.push_back("LongitudeRef");
		metadata_text.push_back("W");
		metadata_key.push_back("GPSLongitude");
		metadata_text.push_back(decimal_to_dms(-lon));
	} else {
		metadata_key.push_back("LongitudeRef");
		metadata_text.push_back("E");
		metadata_key.push_back("GPSLongitude");
		metadata_text.push_back(decimal_to_dms(lon));
	}

	//add data set timestamp
}

void PNGIO::writeMetadata(png_structp &png_ptr, png_infop &info_ptr)
{
	const size_t nr = metadata_key.size();
	png_text *info_text;
	info_text = (png_text *)calloc(sizeof(png_text), nr);
	char **key, **text;
	key = (char**)calloc(sizeof(char)*80, nr);
	text = (char**)calloc(sizeof(char)*80, nr);

	for(size_t ii=0; ii<nr; ii++) {
		key[ii] = (char *)calloc(sizeof(char), 80);
		text[ii] = (char *)calloc(sizeof(char), 80);
		strncpy(key[ii], metadata_key[ii].c_str(), 80);
		strncpy(text[ii], metadata_text[ii].c_str(), 80);
		info_text[ii].key = key[ii];
		info_text[ii].text = text[ii];
		info_text[ii].compression = PNG_TEXT_COMPRESSION_NONE;
	}

	png_set_text(png_ptr, info_ptr, info_text, nr);
	free(info_text);
	for(size_t ii=0; ii<nr; ii++) {
		free(key[ii]);
		free(text[ii]);
	}
	free(key);
	free(text);

	png_write_info(png_ptr, info_ptr);
}

std::string PNGIO::decimal_to_dms(const double& decimal) {
	std::stringstream dms;
	const int d = static_cast<int>( floor(decimal) );
	const double m = floor( ((decimal - (double)d)*60.)*100. ) / 100.;
	const double s = 3600.*(decimal - (double)d) - 60.*m;

	dms << d << "/1 " << static_cast<int>(m*100) << "/100 " << fixed << setprecision(6) << s << "/1";
	return dms.str();
}

void PNGIO::cleanup(FILE *fp, png_structp png_ptr, png_infop info_ptr)
{
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
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
