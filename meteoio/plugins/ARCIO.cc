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
#include "ARCIO.h"
#include <errno.h>
#include <string.h>
#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page arc ARC
 * @section arc_format Format
 * This is for reading grid data in the ARC-GIS format, or more properly, ESRI ascii grid format (see http://en.wikipedia.org/wiki/ESRI_grid). We consider the following specification (in the absence of an official specification):
 * - a single space character is used as field spearator
 * - the header data is right aligned to the 23rd column
 * - float header data has 3 digits precision
 * - all grid data is written as float (which might cause some trouble for some softwares)
 *
 * These specifications should reflect commonly accepted practise.
 *
 * Finally, the naming scheme for meteo grids should be: YYYY-MM-DDTHH.mm_{MeteoGrids::Parameters}.asc
 *
 * @section lus_format Land Use Format
 * The landuse codes are coming from PREVAH and have the format 1LLDC where:
 * - LL is the land use code as given in the table given below
 * - D is the soil depth
 * - C is the field capacity
 *
 * <center><table border="0">
 * <tr><td>
 * <table border="1">
 * <tr><th>land use (vegetation)</th><th>Prevah land use classes</th></tr>
 * <tr><td>01</td><td>water</td></tr>
 * <tr><td>02</td><td>settlement</td></tr>
 * <tr><td>03</td><td>coniferous forest</td></tr>
 * <tr><td>04</td><td>decidous forest</td></tr>
 * <tr><td>05</td><td>mixed forest</td></tr>
 * <tr><td>06</td><td>cereals</td></tr>
 * <tr><td>07</td><td>pasture</td></tr>
 * <tr><td>08</td><td>bush</td></tr>
 * <tr><td>09</td><td>undefined</td></tr>
 * <tr><td>10</td><td>undefined</td></tr>
 * <tr><td>11</td><td>road</td></tr>
 * <tr><td>12</td><td>undefined</td></tr>
 * <tr><td>13</td><td>firn</td></tr>
 * <tr><td>14</td><td>bare ice</td></tr>
 * <tr><td>15</td><td>rock</td></tr>
 * </table></td><td><table border="1">
 * <tr><th>land use (vegetation)</th><th>Prevah land use classes</th></tr>
 * <tr><td>16</td><td>undefined</td></tr>
 * <tr><td>17</td><td>undefined</td></tr>
 * <tr><td>18</td><td>fruit</td></tr>
 * <tr><td>19</td><td>vegetables</td></tr>
 * <tr><td>20</td><td>wheat</td></tr>
 * <tr><td>21</td><td>alpine vegetation</td></tr>
 * <tr><td>22</td><td>wetlands</td></tr>
 * <tr><td>23</td><td>rough pasture</td></tr>
 * <tr><td>24</td><td>subalpine meadow</td></tr>
 * <tr><td>25</td><td>alpine meadow</td></tr>
 * <tr><td>26</td><td>bare soil vegetation</td></tr>
 * <tr><td>27</td><td>free</td></tr>
 * <tr><td>28</td><td>corn</td></tr>
 * <tr><td>29</td><td>grapes</td></tr>
 * <tr><td>30-99</td><td>undefined</td></tr>
 * </table></td></tr>
 * </table></center>
 *
 * @section arc_units Units
 * The distances are assumed to be in meters.
 *
 * @section arc_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - GRID2DPATH: meteo grids directory where to read/write the grids; [Input] and [Output] sections
 * - A3D_VIEW: use Alpine3D's grid viewer naming scheme (default=false)? [Input] and [Output] sections.
 * - DEMFILE: for reading the data as a DEMObject
 * - LANDUSE: for interpreting the data as landuse codes
 * - DAPATH: path+prefix of file containing data assimilation grids (named with ISO 8601 basic date and .sca extension, example ./input/dagrids/sdp_200812011530.sca)
 */

ARCIO::ARCIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	a3d_view_in = false;
	cfg.getValue("A3D_VIEW", "Input", a3d_view_in, Config::nothrow);
	a3d_view_out = false;
	cfg.getValue("A3D_VIEW", "Output", a3d_view_out, Config::nothrow);
	getGridPaths();
}

ARCIO::ARCIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	a3d_view_in = false;
	cfg.getValue("A3D_VIEW", "Input", a3d_view_in, Config::nothrow);
	a3d_view_out = false;
	cfg.getValue("A3D_VIEW", "Output", a3d_view_out, Config::nothrow);
	getGridPaths();
}

ARCIO::ARCIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	a3d_view_in = false;
	cfg.getValue("A3D_VIEW", "Input", a3d_view_in, Config::nothrow);
	a3d_view_out = false;
	cfg.getValue("A3D_VIEW", "Output", a3d_view_out, Config::nothrow);
	getGridPaths();
}

void ARCIO::getGridPaths() {
	grid2dpath_in="", grid2dpath_out="";
	string tmp="";
	cfg.getValue("GRID2D", "Input", tmp, Config::nothrow);
	if (tmp == "ARC") //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
	tmp="";
	cfg.getValue("GRID2D", "Output", tmp, Config::nothrow);
	if (tmp == "ARC") //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Output", grid2dpath_out);
}

ARCIO::~ARCIO() throw()
{
	cleanup();
}

void ARCIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void ARCIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name)
{
	int i_ncols, i_nrows;
	unsigned int ncols, nrows;
	double xllcorner, yllcorner, cellsize, plugin_nodata;
	double tmp;
	std::string line="";
	std::map<std::string, std::string> header; // A map to save key value pairs of the file header

	if (!IOUtils::validFileName(full_name)) {
		throw InvalidFileNameException(full_name, AT);
	}
	if (!IOUtils::fileExists(full_name)) {
		throw FileNotFoundException(full_name, AT);
	}

	fin.clear();
	fin.open (full_name.c_str(), ifstream::in);
	if (fin.fail()) {
		stringstream ss;
		ss << "Error openning file \"" << full_name << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	//Go through file, save key value pairs
	try {
		IOUtils::readKeyValueHeader(header, fin, 6, " ");
		IOUtils::getValueForKey(header, "ncols", i_ncols);
		IOUtils::getValueForKey(header, "nrows", i_nrows);
		IOUtils::getValueForKey(header, "xllcorner", xllcorner);
		IOUtils::getValueForKey(header, "yllcorner", yllcorner);
		IOUtils::getValueForKey(header, "cellsize", cellsize);
		IOUtils::getValueForKey(header, "NODATA_value", plugin_nodata);

		//HACK!! would it be possible for getValueForKey() to do this transparently? (with a user flag)
		i_ncols = IOUtils::standardizeNodata(i_ncols, plugin_nodata);
		i_nrows = IOUtils::standardizeNodata(i_nrows, plugin_nodata);
		xllcorner = IOUtils::standardizeNodata(xllcorner, plugin_nodata);
		yllcorner = IOUtils::standardizeNodata(yllcorner, plugin_nodata);
		cellsize = IOUtils::standardizeNodata(cellsize, plugin_nodata);

		if ((i_ncols==0) || (i_nrows==0)) {
			throw IOException("Number of rows or columns in 2D Grid given is zero, in file: " + full_name, AT);
		}
		if((i_ncols<0) || (i_nrows<0)) {
			throw IOException("Number of rows or columns in 2D Grid read as \"nodata\", in file: " + full_name, AT);
		}
		ncols = (unsigned int)i_ncols;
		nrows = (unsigned int)i_nrows;

		//compute/check WGS coordinates (considered as the true reference) according to the projection as defined in cfg
		Coords location(coordin, coordinparam);
		location.setXY(xllcorner, yllcorner, IOUtils::nodata);

		//Initialize the 2D grid
		grid_out.set(ncols, nrows, cellsize, location);

		//Read one line after the other and parse values into Grid2DObject
		for (unsigned int kk=nrows-1; (kk < nrows); kk--) {
			getline(fin, line, eoln);
			std::istringstream iss(line);
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);

			for (unsigned int ll=0; ll < ncols; ll++){
				iss >> std::skipws >> tmp;
				if (iss.fail()) {
					stringstream ss;
					ss << "Can not read column " << ll+1 << " of data line " << nrows-kk << " in file " << full_name << ": ";
					ss << ncols << " columns of doubles expected";
					throw InvalidFormatException(ss.str(), AT);
				}
				grid_out.grid2D(ll, kk) = IOUtils::standardizeNodata(tmp, plugin_nodata);
			}

		}
	} catch(const std::exception& e) {
		cleanup();
		std::stringstream msg;
		msg << "[E] Error when reading ARC grid \"" << full_name << "\" : " << e.what();
		throw InvalidFormatException(msg.str(), AT);
	}
	cleanup();
}

void ARCIO::read2DGrid(Grid2DObject& grid_out, const std::string& filename) {
	read2DGrid_internal(grid_out, grid2dpath_in+"/"+filename);
}

void ARCIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if(a3d_view_in) {
		string ext="";
		if(parameter==MeteoGrids::HS)
			ext="sdp";
		else {
			ext = MeteoGrids::getParameterName(parameter);
			IOUtils::toLower(ext);
		}
		read2DGrid_internal(grid_out, grid2dpath_in + "/" + date.toString(Date::NUM)+"."+ext );
	} else {
		std::string date_str = date.toString(Date::ISO);
		std::replace( date_str.begin(), date_str.end(), ':', '.');
		read2DGrid_internal(grid_out, grid2dpath_in + "/" + date_str + "_" + MeteoGrids::getParameterName(parameter) + ".asc");
	}
}

void ARCIO::readDEM(DEMObject& dem_out)
{
	string filename="";
	cfg.getValue("DEMFILE", "Input", filename);
	read2DGrid_internal(dem_out, filename);
}

void ARCIO::readLanduse(Grid2DObject& landuse_out)
{
	string filename="";
	cfg.getValue("LANDUSEFILE", "Input", filename);
	read2DGrid_internal(landuse_out, filename);
}

void ARCIO::readAssimilationData(const Date& date_in, Grid2DObject& da_out)
{
	int yyyy, MM, dd, hh, mm;
	date_in.getDate(yyyy, MM, dd, hh, mm);
	string filepath="";

	cfg.getValue("DAPATH", "Input", filepath);

	stringstream ss;
	ss.fill('0');
	ss << filepath << "/" << setw(4) << yyyy << setw(2) << MM << setw(2) <<  dd << setw(2) <<  hh <<  setw(2) <<  mm << ".sca";

	read2DGrid_internal(da_out, ss.str());
}

void ARCIO::readStationData(const Date&, std::vector<StationData>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARCIO::readMeteoData(const Date&, const Date&, std::vector< std::vector<MeteoData> >&,
                          const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARCIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARCIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARCIO::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	std::string full_name = grid2dpath_out+"/"+name;
	fout.open(full_name.c_str());
	if (fout.fail()) {
		stringstream ss;
		ss << "Error openning file \"" << full_name << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	try {
		Coords llcorner=grid_in.llcorner;
		//we want to make sure that we are using the provided projection parameters
		//so that we output is done in the same system as the inputs
		llcorner.setProj(coordout, coordoutparam);

		fout << fixed << showpoint << setprecision(6);
		fout << "ncols " << setw(23-6) << grid_in.ncols << "\n";
		fout << "nrows " << setw(23-6) << grid_in.nrows << "\n";
		fout << "xllcorner " << setw(23-10) << setprecision(3) << llcorner.getEasting() << "\n";
		fout << "yllcorner " << setw(23-10) << setprecision(3) << llcorner.getNorthing() << "\n";
		fout << "cellsize " << setw(23-9) << setprecision(3) << grid_in.cellsize << "\n";
		fout << "NODATA_value " << (int)(IOUtils::nodata) << "\n";

		for (unsigned int kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
			for (unsigned int ll=0; ll < grid_in.ncols; ll++){
				fout << grid_in.grid2D(ll, kk) << " ";
			}
			fout << "\n";
		}
	} catch(...) {
		cout << "[E] error when writing ARC grid \"" << full_name << "\" " << AT << ": "<< endl;
		cleanup();
		throw;
	}

	cleanup();
}

void ARCIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	//the path will be added by write2DGrid
	if(a3d_view_out) {
		string ext="";
		if(parameter==MeteoGrids::HS)
			ext="sdp";
		else {
			ext = MeteoGrids::getParameterName(parameter);
			IOUtils::toLower(ext);
		}
		write2DGrid(grid_in, date.toString(Date::NUM)+"."+ext );
	} else {
		std::string date_str = date.toString(Date::ISO);
		std::replace( date_str.begin(), date_str.end(), ':', '.');
		write2DGrid(grid_in, date_str + "_" + MeteoGrids::getParameterName(parameter) + ".asc");
	}
}

#ifndef _METEOIO_JNI
extern "C"
{
#define COMPILE_PLUGIN
#include "exports.h"

	METEOIO_EXPORT void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	METEOIO_EXPORT void* loadObject(const std::string& classname, const Config& cfg) {
		if(classname == "ARCIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new ARCIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
