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
#include "PGMIO.h"

using namespace std;

namespace mio {
/**
 * @page pgmio PGMIO
 * @section pgmio_format Format
 * This reads a grid file in PGM format (see http://www.fileformat.info/format/pbm/egff.htm). This is a graphic format that is supported by a wide range of graphics programs (Gimp, Irfanview, Paint Shop Pro, gqview, etc). This allows to write a grid as an image (one pixel equals one cell), read an image as a grid (useful for creating synthetic DEMs). Since there is no geolocalization information in this format, such data is either encoded as a comment (when writing a file) a read from io.ini (for reading).
 *
 * Please keep in mind that only a finite number of greyscales are used, making a discretization of the data. Moreover, we consider that a color of "0" is NODATA.
 *
 * @section pgmio_units Units
 * Cellsize is in meters, x/y coords in the section's coordinate system.
 *
 * @section pgmio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - PGM_XCOORD: lower left x coordinate; [Input] section
 * - PGM_YCOORD: lower left y coordinate; [Input] section
 * - PGM_CELLSIZE: cellsize in meters; [Input] section
 * - PGM_MIN: minimum value in real world coordinates to match with the minimum value read out of the PGM file (such minimum being greater than 0 because 0 is NODATA)
 * - PGM_MAX: maximum value in real world coordinates to match with the maximum value read out of the PGM file
 */

const double PGMIO::plugin_nodata = -999.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

PGMIO::PGMIO(void (*delObj)(void*), const Config& i_cfg) : IOInterface(delObj), cfg(i_cfg)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PGMIO::PGMIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PGMIO::PGMIO(const Config& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

PGMIO::~PGMIO() throw()
{

}

size_t PGMIO::getNextHeader(std::vector<std::string>& vecString, const std::string& filename) {
	std::string line="";

	while(!fin.eof()) {
		getline(fin, line);
		IOUtils::trim(line);
		if(line.size()>0 && line.at(0)!='#') {
			return IOUtils::readLineToVec(line, vecString);
		}
	}
	throw IOException("Can not read necessary header lines in " + filename, AT);
}

void PGMIO::read2DGrid(Grid2DObject& grid_out, const std::string& filename)
{
	unsigned int ncols, nrows, nr_colors;
	double xllcorner, yllcorner, cellsize;
	const double plugin_nodata=0.;
	double tmp_val, val_min, val_max;
	std::vector<std::string> tmpvec;
	std::string line="";

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

	const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	//Go through file, save key value pairs
	try {
		//read header: magic value
		if(getNextHeader(tmpvec, filename)!=1) {
			throw IOException("Can not read necessary header in " + filename, AT);
		}
		//read header: image width and height
		if(getNextHeader(tmpvec, filename)!=2) {
			throw IOException("Can not read necessary header in " + filename, AT);
		}
		IOUtils::convertString(ncols, tmpvec[0]);
		IOUtils::convertString(nrows, tmpvec[1]);
		//read header: number of greys
		if(getNextHeader(tmpvec, filename)!=1) {
			throw IOException("Can not read necessary header in " + filename, AT);
		}
		IOUtils::convertString(nr_colors, tmpvec[0]);

		cfg.getValue("PGM_XCOORD", "Input", xllcorner, Config::dothrow);
		cfg.getValue("PGM_YCOORD", "Input", yllcorner, Config::dothrow);
		cfg.getValue("PGM_CELLSIZE", "Input", cellsize, Config::dothrow);
		cfg.getValue("PGM_MIN", "Input", val_min, Config::dothrow);
		cfg.getValue("PGM_MAX", "Input", val_max, Config::dothrow);

		Coords location(coordin, coordinparam);
		location.setXY(xllcorner, yllcorner, IOUtils::nodata);

		//Initialize the 2D grid
		grid_out.set(ncols, nrows, cellsize, location);

		//initialize scale factor
		const double scale_factor = (val_max-val_min)/(double)(nr_colors-2); //because 256 colors = 0 to 255!! and color0 = nodata

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

				if(tmp_val==plugin_nodata) {
					//replace file's nodata by uniform, internal nodata
					grid_out.grid2D(ll, kk) = IOUtils::nodata;
				} else {
					grid_out.grid2D(ll, kk) = (tmp_val-1)*scale_factor+val_min; //because color0 = nodata
				}
			}
		}
	} catch(std::exception& e) {
		cleanup();
		throw;
	}
	cleanup();
}

void PGMIO::readDEM(DEMObject& dem_out)
{
	string filename="";
	cfg.getValue("DEMFILE", "Input", filename);
	read2DGrid(dem_out, filename);
}

void PGMIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                          std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                          const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                           const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	const unsigned int nr_colors = 256;
	fout.open(name.c_str());
	if (fout.fail()) {
		throw FileAccessException(name, AT);
	}

	Coords llcorner=grid_in.llcorner;
	//we want to make sure that we are using the provided projection parameters
	//so that we output is done in the same system as the inputs
	llcorner.setProj(coordout, coordoutparam);

	fout << fixed << showpoint << setprecision(6);

	try {
		//writing the header
		fout << "P2\n";
		fout << "#Generated by MeteoIO - http://slfsmm.indefero.net/p/meteoio\n";
		fout << "#llcorner latitude = " << setprecision(6) << llcorner.getLat() << "\n";
		fout << "#llcorner longitude = " << setprecision(6) << llcorner.getLon() << "\n";
		fout << "#cellsize = " << setprecision(2) << grid_in.cellsize << " m\n";
		fout << grid_in.ncols << " " << grid_in.nrows << "\n";
		fout << nr_colors << "\n";

		//writing the data
		const double max_value = grid_in.grid2D.getMax();
		const double min_value = grid_in.grid2D.getMin();
		const double scaling = 1./(max_value - min_value) * (double)(nr_colors-1); //so we keep color 0 for nodata

		for (unsigned int kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
			for (unsigned int ll=0; ll < grid_in.ncols; ll++) {
				const double value = grid_in.grid2D(ll, kk);
				if(value!=IOUtils::nodata)
					fout << (int)floor((grid_in.grid2D(ll, kk)-min_value)*scaling)+1 << " ";
				else
					fout << "0" << " ";
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

void PGMIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
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

	METEOIO_EXPORT void* loadObject(const string& classname, const Config& cfg) {
		if(classname == "PGMIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new PGMIO(deleteObject, cfg);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
#endif

} //namespace
