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

#define ARPS_MAX_LINE_LENGTH 6000
#define ARPS_MAX_STRING_LENGTH 256

using namespace std;

namespace mio {
/**
 * @page arps ARPSIO
 * @section arps_format Format
 * This is for reading grid data in the ARPS grid format
 *
 * @section arps_units Units
 *
 *
 * @section arps_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - METEOPATH: path to the output directory; [Output] section
 * - METEOFILE#: input meteo data file, e.g. METEOFILE1, METEOFILE2; [Input] section
 * - STATION#: station name as listed in the METAFILE, e.g. STATION1, STATION2; [Input] section
 * - METAFILE: filename of the meta data file; [Input] section
 * - NROFSTATIONS: integer, the number of stations for which meteo files are provided; [Input] section
 */

const double ARPSIO::plugin_nodata = -999.; //plugin specific nodata value

ARPSIO::ARPSIO(void (*delObj)(void*), const std::string& filename) : IOInterface(delObj), cfg(filename)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	dimx=dimy=dimz=0;
	fin=NULL;
}

ARPSIO::ARPSIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	dimx=dimy=dimz=0;
	fin=NULL;
}

ARPSIO::ARPSIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	dimx=dimy=dimz=0;
	fin=NULL;
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

void ARPSIO::read3DGrid(Grid3DObject& grid_out, const std::string& /*_name*/)
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
				fscanf(fin," %16lf%*[\n]",&tmp);
				grid_out.grid3D(ix,iy,iz) = tmp;
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
	readGridLayer(std::string("zp_coordinat"), 1 ,dem_out);
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

void ARPSIO::openGridFile(const std::string& _filename)
{
	filename = _filename;

	if((fin=fopen(filename.c_str(),"r")) == 0) {
		cleanup();
		throw FileAccessException("Can not open file "+filename, AT);
	}

	//now, initialize the sizes
	char dummy[ARPS_MAX_LINE_LENGTH];
	for (int j=0; j<12; j++) {
		fgets(dummy,ARPS_MAX_STRING_LENGTH,fin);
	}
	if (fscanf(fin,"%ud %ud %ud\n",&dimx,&dimy,&dimz)!=3) {
		cleanup();
		throw InvalidFormatException("Can not read dimx, dimy, dimz from file "+filename, AT);
	}
	if (dimx==0 || dimy==0 || dimz==0) {
		cleanup();
		throw IndexOutOfBoundsException("Invalid dimx, dimy, dimz from file "+filename, AT);
	}
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
	grid.grid2D.resize(dimx, dimy);

	// Read until the parameter is found
	moveToMarker(parameter);

	// move to the begining of the layer of interest
	if(layer>1) {
		double tmp;
		const unsigned int jmax=dimx*dimy*(layer-1);
		for (unsigned int j = 0; j < jmax; j++)
			fscanf(fin," %16lf%*[\n]",&tmp);
	}

	//read the data we are interested in
	for (unsigned int ix = 0; ix < dimx; ix++) {
		for (unsigned int iy = 0; iy < dimy; iy++) {
			double tmp;
			fscanf(fin," %16lf%*[\n]",&tmp);
			grid.grid2D(ix,iy) = tmp;
		}
	}
}

void ARPSIO::moveToMarker(const std::string& marker)
{
	char dummy[ARPS_MAX_LINE_LENGTH];
	do {
		fscanf(fin," %s ",dummy);
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

	void* loadObject(const string& classname, const string& filename) {
		if(classname == "ARPSIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new ARPSIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
