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
 * - COORDIN: input coordinate system (see Coordinate)
 * - COORDPARAM: extra input coordinates parameters (see Coordinate)
 * - DEMFILE: for reading the data as a DEMObject
 * - LANDUSE: for interpreting the data as landuse codes
 * - DAPATH: path+prefix of file containing data assimilation grids (named with ISO 8601 basic date and .sca extension, example ./input/dagrids/sdp_200812011530.sca)
 */

using namespace std;

void ARCIO::getProjectionParameters() {
	//get projection parameters
	try {
		cfg.getValue("COORDSYS", "INPUT", coordsys);
		cfg.getValue("COORDPARAM", "INPUT", coordparam, ConfigReader::nothrow);
	} catch(std::exception& e){
		//problems while reading values for COORDIN or COORDPARAM
		std::cerr << "[E] " << AT << ": reading configuration file: " << "\t" << e.what() << std::endl;
		throw;
	}
}

ARCIO::ARCIO(void (*delObj)(void*), const std::string& filename) : IOInterface(delObj), cfg(filename)
{
	getProjectionParameters();
}

ARCIO::ARCIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	getProjectionParameters();
}

ARCIO::ARCIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	getProjectionParameters();
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

void ARCIO::read2DGrid(Grid2DObject& grid_out, const std::string& filename)
{
	int i_ncols, i_nrows;
	unsigned int ncols, nrows;
	double xllcorner, yllcorner, cellsize, plugin_nodata;
	double tmp_val;
	std::vector<std::string> tmpvec;
	std::string line="";
	std::map<std::string, std::string> header; // A map to save key value pairs of the file header

	if (!IOUtils::validFileName(filename)) {
		throw InvalidFileNameException(filename, AT);
	}
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}
  
	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
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
			throw IOException("Number of rows or columns in 2D Grid given is zero, in file: " + filename, AT);
		}
		if((i_ncols<0) || (i_nrows<0)) {
			throw IOException("Number of rows or columns in 2D Grid read as \"nodata\", in file: " + filename, AT);
		}
		ncols = (unsigned int)i_ncols;
		nrows = (unsigned int)i_nrows;

		//compute/check WGS coordinates (considered as the true reference) according to the projection as defined in cfg
		Coords location(coordsys, coordparam);
		location.setXY(xllcorner, yllcorner, IOUtils::nodata);
		
		//Initialize the 2D grid
		grid_out.set(ncols, nrows, cellsize, location);
		
		//Read one line after the other and parse values into Grid2DObject
		for (unsigned int kk=nrows-1; (kk < nrows); kk--) {
			getline(fin, line, eoln); //read complete line
			//cout << "k:" << kk << "\n" << line << endl;

			if (IOUtils::readLineToVec(line, tmpvec) != ncols) {
				throw InvalidFormatException("Premature End " + filename, AT);
			}
			
			for (unsigned int ll=0; ll < ncols; ll++){
				if (!IOUtils::convertString(tmp_val, tmpvec[ll], std::dec)) {
					throw ConversionFailedException("For Grid2D value in line: " + line + " in file " + filename, AT);
				}
				
				if(tmp_val<=plugin_nodata) {
					//replace file's nodata by uniform, internal nodata
					grid_out.grid2D(ll, kk) = IOUtils::nodata;
				} else {
					grid_out.grid2D(ll, kk) = tmp_val;
				}
			}
		}
	} catch(std::exception& e) {
		cleanup();
		throw;
	}
	cleanup();
}

void ARCIO::readDEM(DEMObject& dem_out)
{
	string filename="";
	cfg.getValue("DEMFILE", "INPUT", filename);
	read2DGrid(dem_out, filename);
}

void ARCIO::readLanduse(Grid2DObject& landuse_out)
{
	string filename="";
	cfg.getValue("LANDUSEFILE", "INPUT", filename); // cout << tmp << endl;
	read2DGrid(landuse_out, filename);
}

void ARCIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	int yyyy, MM, dd, hh, mm;
	date_in.getDate(yyyy, MM, dd, hh, mm);
	string filepath="";

	cfg.getValue("DAPATH", "INPUT", filepath); // cout << tmp << endl;
  
	stringstream ss;
	ss.fill('0');
	ss << filepath << "/" << setw(4) << yyyy << setw(2) << MM << setw(2) <<  dd << setw(2) <<  hh <<  setw(2) <<  mm << ".sca";

	read2DGrid(da_out, ss.str());
}

void ARCIO::readStationData(const Date_IO&, std::vector<StationData>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARCIO::readMeteoData(const Date_IO&, const Date_IO&, 
					 std::vector< std::vector<MeteoData> >&, 
					 std::vector< std::vector<StationData> >&,
					 const unsigned int&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARCIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&, 
					  const std::vector< std::vector<StationData> >&,
					  const std::string&)
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
	//uncomment if no overwriting should be allowed
	/* 
	 if (IOUtils::fileExists(name)) {
	   throw IOException("File " + name + " already exists!", AT);
	 }
	*/

	Coords llcorner=grid_in.llcorner;
	//we want to make sure that we are using the provided projection parameters
	//so that we output is done in the same system as the inputs
	llcorner.setProj(coordsys, coordparam);

	fout.open(name.c_str());
	if (fout.fail()) {
		throw FileAccessException(name, AT);
	}
	
	fout << fixed << showpoint << setprecision(6);

	try {
		fout << "ncols " << setw(23-6) << grid_in.ncols << endl;
		fout << "nrows " << setw(23-6) << grid_in.nrows << endl;
		fout << "xllcorner " << setw(23-10) << setprecision(3) << llcorner.getEasting() << endl;
		fout << "yllcorner " << setw(23-10) << setprecision(3) << llcorner.getNorthing() << endl;
		fout << "cellsize " << setw(23-9) << setprecision(3) << grid_in.cellsize << endl;
		fout << "NODATA_value " << (int)(IOUtils::nodata) << endl;

		for (unsigned int kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
			for (unsigned int ll=0; ll < grid_in.ncols; ll++){
				fout << grid_in.grid2D(ll, kk) << " ";
			}
			fout << endl;
		}
	} catch(...) {
		cout << "[E] " << AT << ": "<< endl;
		cleanup();
		throw;
	}

	cleanup();
}

#ifndef _METEOIO_JNI
extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
  
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "ARCIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new ARCIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif
