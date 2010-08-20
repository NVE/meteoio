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
#include "ARPSIO.h"

using namespace std;

namespace mio {
/**
 * @page arps ARPSIO
 * @section arps_format Format
 * This is for reading grid data in the ARPS grid format (it transparently supports both true ARPS ascii grids and grids modified by the ARPSGRID utility). Currently, only DEM reading is implemented.
 *
 * @section arps_units Units
 *
 *
 * @section arps_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - DEMFILE: path and file containing the DEM; [Input] section
 * - ARPS_X: x coordinate of the lower left corner of the grids; [Input] section
 * - ARPS_Y: y coordinate of the lower left corner of the grids; [Input] section
 * - ARPSPATH: path to the input directory where to find the arps files to be read as grids; [Input] section //NOT USED YET
 */

const double ARPSIO::plugin_nodata = -999.; //plugin specific nodata value

ARPSIO::ARPSIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	dimx=dimy=dimz=0;
	cellsize=0.;
	fin=NULL;
	is_true_arps=true;
}

ARPSIO::ARPSIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	dimx=dimy=dimz=0;
	cellsize=0.;
	fin=NULL;
	is_true_arps=true;
}

ARPSIO::ARPSIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	dimx=dimy=dimz=0;
	cellsize=0.;
	fin=NULL;
	is_true_arps=true;
}

ARPSIO::~ARPSIO() throw()
{
	cleanup();
}

void ARPSIO::read2DGrid(Grid2DObject& grid_out, const std::string& /*_name*/)
{
	std::string meteopathname;
	std::string parameter, grid_name;
	unsigned int layer;
	//HACK: build parameter and grid_name from _name
	//example: _name="NW08::u" -> grid_name=NW08 and parameter=u
	//idem for layer?
	cfg.getValue("ARPSPATH", "Input", meteopathname);
	const std::string _filename = meteopathname + grid_name + ".asc";

	openGridFile(_filename);
	readGridLayer(parameter, layer, grid_out);

	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::read3DGrid(Grid3DObject& grid_out, const std::string& /*in_name*/)
{
	std::string meteopathname;
	std::string parameter, grid_name;
	//HACK: build parameter and grid_name from _name
	//example: _name="NW08::u" -> grid_name=NW08 and parameter=u
	cfg.getValue("ARPSPATH", "Input", meteopathname);
	const std::string _filename = meteopathname + grid_name + ".asc";

	openGridFile(_filename);

	//resize the grid just in case
	grid_out.grid3D.resize(dimx, dimy, dimz);

	// Read until the parameter is found
	moveToMarker(parameter);

	//read the data we are interested in
	for (unsigned int ix = 0; ix < dimx; ix++) {
		for (unsigned int iy = 0; iy < dimy; iy++) {
			for (unsigned int iz = 0; iz < dimz; iz++) {
				double tmp;
				if(fscanf(fin," %16lf%*[\n]",&tmp)==1) {
					grid_out.grid3D(ix,iy,iz) = tmp;
				} else {
					cleanup();
					throw InvalidFormatException("Failure in reading 3D grid in file "+_filename, AT);
				}
			}
		}
	}

	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readDEM(DEMObject& dem_out)
{
	std::string _filename;
	cfg.getValue("DEMFILE", "Input", _filename);
	openGridFile(_filename);
	if(is_true_arps) {
		readGridLayer(std::string("zp coordinat"), 1 ,dem_out);
	} else {
		readGridLayer(std::string("zp_coordinat"), 1 ,dem_out);
	}
}

void ARPSIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
					 std::vector< std::vector<MeteoData> >& /*vecMeteo*/, 
					 std::vector< std::vector<StationData> >& /*vecStation*/,
					 const unsigned int&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                          const std::vector< std::vector<StationData> >& /*vecStation*/,
                          const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void ARPSIO::initializeGRIDARPS()
{
	double v1, v2;

	//go to read the sizes
	moveToMarker("nnx");
	//finish reading the line and move to the next one
	if(fscanf(fin,"%*[^\n]")!=0) {
		cleanup();
		throw InvalidFormatException("Error in file format of file "+filename, AT);
	}
	if (fscanf(fin," %u %u %u \n",&dimx,&dimy,&dimz)!=3) {
		cleanup();
		throw InvalidFormatException("Can not read dimx, dimy, dimz from file "+filename, AT);
	}
	if (dimx==0 || dimy==0 || dimz==0) {
		cleanup();
		throw IndexOutOfBoundsException("Invalid dimx, dimy, dimz from file "+filename, AT);
	}

	//initializing cell size
	moveToMarker("x_coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two x coordinates from file "+filename, AT);
	}
	const double cellsize_x = v2 - v1;
	moveToMarker("y_coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two y coordinates from file "+filename, AT);
	}
	const double cellsize_y = v2 - v1;
	if(cellsize_x!=cellsize_y) {
		cleanup();
		throw InvalidFormatException("Only square cells currently supported! Non compliance in file "+filename, AT);
	}
	cellsize = cellsize_y;
}

void ARPSIO::initializeTrueARPS(const char curr_line[ARPS_MAX_LINE_LENGTH])
{
	double v1, v2;

	//go to read the sizes
	if (sscanf(curr_line," nx = %u, ny = %u, nz = %u ",&dimx,&dimy,&dimz)!=3) {
		cleanup();
		throw InvalidFormatException("Can not read dimx, dimy, dimz from file "+filename, AT);
	}
	if (dimx==0 || dimy==0 || dimz==0) {
		cleanup();
		throw IndexOutOfBoundsException("Invalid dimx, dimy, dimz from file "+filename, AT);
	}

	//initializing cell size
	moveToMarker("x coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two x coordinates from file "+filename, AT);
	}
	const double cellsize_x = v2 - v1;
	moveToMarker("y coordinate");
	if (fscanf(fin,"%lg %lg",&v1,&v2)!=2) {
		cleanup();
		throw InvalidFormatException("Can not read first two y coordinates from file "+filename, AT);
	}
	const double cellsize_y = v2 - v1;
	if(cellsize_x!=cellsize_y) {
		cleanup();
		throw InvalidFormatException("Only square cells currently supported! Non compliance in file "+filename, AT);
	}
	cellsize = cellsize_y;
}

void ARPSIO::openGridFile(const std::string& in_filename)
{
	unsigned int v1;
	filename = in_filename;

	if((fin=fopen(filename.c_str(),"r")) == NULL) {
		cleanup();
		throw FileAccessException("Can not open file "+filename, AT);
	}

	//identify if the file is an original arps file or a file modified by ARPSGRID
	char dummy[ARPS_MAX_LINE_LENGTH];
	for (int j=0; j<5; j++) {
		//the first easy difference in the structure happens at line 5
		if(fgets(dummy,ARPS_MAX_STRING_LENGTH,fin)==NULL) {
			cleanup();
			throw InvalidFormatException("Fail to read header lines of file "+filename, AT);
		}
	}
	if (sscanf(dummy," nx = %u, ny = ", &v1)<1) {
		//this is an ASCII file modified by ARPSGRID
		is_true_arps=false;
		initializeGRIDARPS();
	} else {
		//this is a true ARPS file
		initializeTrueARPS(dummy);
	}

	//get llcorner
	cfg.getValue("ARPS_X", "Input", xcoord, Config::dothrow);
	cfg.getValue("ARPS_Y", "Input", ycoord, Config::dothrow);

	//come back to the begining of the file
	rewind(fin);
}

void ARPSIO::cleanup() throw()
{
	if (fin!=NULL) {//close fin if open
		fclose(fin);
		fin=NULL;
	}
	/*if (fin.is_open()) {//close fin if open
		fin.close();
	}*/
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

/** @brief Read a specific layer for a given parameter from the ARPS file
 * @param parameter The parameter to extract. This could be any of the following:
 *        - x_coordinate for getting the X coordinates of the mesh
 *        - y_coordinate for getting the Y coordinates of the mesh
 *        - zp_coordinat for getting the Z coordinates of the mesh
 *        - u for getting the u component of the wind field
 *        - v for getting the v component of the wind field
 *        - w for getting the w component of the wind field
 * @param layer     Index of the layer to extract (1 to dimz)
 * @param grid      [out] grid containing the values. The grid will be resized if necessary.
*/
void ARPSIO::readGridLayer(const std::string& parameter, const unsigned int& layer, Grid2DObject& grid)
{
	if(layer<1 || layer>dimz) {
		cleanup();
		stringstream tmp;
		tmp << "Layer " << layer << " does not exist in ARPS file " << filename << " (nr layers=" << dimz << " ";
		throw IndexOutOfBoundsException(tmp.str(), AT);
	}

	//resize the grid just in case
	Coords llcorner(coordin, coordinparam);
	llcorner.setXY(xcoord, ycoord, IOUtils::nodata);
	grid.set(dimx, dimy, cellsize, llcorner);

	// Read until the parameter is found
	moveToMarker(parameter);

	// move to the begining of the layer of interest
	if(layer>1) {
		double tmp;
		const unsigned int jmax=dimx*dimy*(layer-1);
		for (unsigned int j = 0; j < jmax; j++)
			if(fscanf(fin," %16lf%*[\n]",&tmp)==EOF) {
				cleanup();
				throw InvalidFormatException("Fail to skip data layers in file "+filename, AT);
			}
	}

//HACK TODO swap ix/iy (for speed) and number cells from llcorner (cf ARCIO)
	//read the data we are interested in
	for (unsigned int ix = 0; ix < dimx; ix++) {
		for (unsigned int iy = 0; iy < dimy; iy++) {
			double tmp;
			if(fscanf(fin," %16lf%*[\n]",&tmp)==1) {
				grid.grid2D(ix,iy) = tmp;
			} else {
				cleanup();
				throw InvalidFormatException("Fail to read data layer in file "+filename, AT);
			}
		}
	}
}

void ARPSIO::moveToMarker(const std::string& marker)
{
	char dummy[ARPS_MAX_LINE_LENGTH];
	do {
		fscanf(fin," %[^\t\n] ",dummy);
	} while (!feof(fin) && strcmp(dummy,marker.c_str()) != 0);
	if(feof(fin)) {
		cleanup();
		const std::string message = "End of file "+filename+" should NOT have been reached when looking for "+marker;
		throw InvalidFormatException(message, AT);
	}
}

#ifndef _METEOIO_JNI
extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const std::string& classname, const Config& cfg) {
		if(classname == "ARPSIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new ARPSIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
