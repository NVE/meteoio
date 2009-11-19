#ifndef __MAINPAGE_H__
#define __MAINPAGE_H__
 /**
 * @mainpage Welcome to MeteoIO
 * @section intro_sec Introduction
 * This library aims at making data access easy and safe for numerical simulations in environmental sciences requiring general meteorological data. It's main design goals are:
 * - providing data format/protocol independent data access
 * - providing safe and robust I/O
 * - making I/O code as unobtrusive and simple as possible for the user
 * - providing ready to use data to the user, which means transparent caching, filtering, resampling, spatial interpolation.
 * - enabling unattended use from an IO point of view
 * - offering high modularity so that individual elements of the library can easily be replaced/expanded/added
 * - by its modularity, help interdisciplinary development, each module being targeted at a specific developer profile
 * 
 * This library is available under LPGL version 3 or above, see <a href="http://www.gnu.org/licenses/lgpl.txt">www.gnu.org</a>.
 *
 * @section table_of_content Table of content
 * -# \subpage quick_overview "Quick overview" of the functionnality provided by MeteoIO
 * -# \subpage plugins "Available plugins" and usage
 * -# \subpage filters "Available filters" and usage
 * -# \subpage dev_plugins "Plugins" developer's guide
 * -# \subpage dev_filters "Filters" developer's guide
 * -# \subpage examples "Examples"
 */

 /**
 * @page quick_overview Quick overview
 * This library contains various classes that have been designed to deal with various sets of problems. This page shows the different sets of problems and what kind of functionnality the library offers to tackle them.
 * 
 *
 * @section iohandler_sec Data reading
 * The class IOHandler provides the meteorological data from the sources selected by the user in its configuration file. This class inherits from IOInterface and is implemented through plugins that are responsible for implementing a given data access (see \ref dev_plugins "Plugins developer's guide" for more information). It therefore proposes a uniform, standardized access to the data that can be meteorological data, gridded data (including Digital Elevation Model (DEM) data or variations like for landuse codes) and tables of coordinates (for special processing at users selected locations). A buffered version of this class exists: BufferedIOHandler that should be prefered. The description of the plugins and their usage can be found in \ref plugins "Available plugins".
 * This class also transparently calls the filtering class, FilterAlgorithms in order to filter the data according to the configuration of the user.
 * 
 *
 * @section meteo Meteorological data
 * The data structures designed for storing meteorological data have been split into two classes: MeteoData and StationData.
 * @subsection meteodata MeteoData
 * The class MeteoData stores the measurement data coming from some idealized station. It contains the widest set of meteorological measurements. It can be compared, assigned and set (either using a constructor or by calling a set method). Its meteorological parameters can be directly accessed or using a param() method that takes a enum in order to be able to cycle through the parameters.
 * @subsection stationdata StationData
 * The class StationData contains the metadata of a weather station, that is mostly its location. It supports the comparison operators.
 * @subsection getmeteo_sec Getting the data
 * The getMeteoData method defined in the IOHandler class provides a vector of MeteoData and StationData for the requested time step. More details are given in \ref iohandler_sec .
 * 
 *
 * @section arrays Arrays related functionnalities
 * @subsection arrays_sec Arrays
 * The classes Array, Array2D and Array3D are designed for the storage and handling of 1D, 2D, 3D arrays in memory. These classes provide access to a given element (read/write), sizing or resizing of an existing array as well as clearing an array. They also provide the minimum and the maximum of the values that are stored in the array. Finally, a subset of an array can be extracted.
 * @subsection grids_sec Grids
 * Built on top of the arrays, defined as classes Grid2DObject and Grid3DObject, the grids add the geolocalization. This means that the coordinates of the lower-left corner of the array are stored as well as the cellsize. They can be built manually, or by providing an array. A subset constructor is available, allowing to extract a subset of the grid. It is also possible to get the lat/long (in WGS84) coordinates matching an (i,j) coordinate in the grid. Finally, It is possible to test for geolocalization equality (ie: do two grids have the same geolocalization).
 * @subsection dem_sec Digital Elevation Models
 * The last layer for gridded data is class DEMObject. Various parameters that are specific to Digital Elevation Models (DEM) are added: for each grid point, the slope, the azimuth, the curvature as well as the normal vector are defined (an optional parameter can be used to select the algorithm to be used). The minimums and maximums (over the grid) for each of these parameters are available. A subset of the DEM can be extracted using the subset constructor.
 * 
 *
 * @section proj_sec Geographic projections
 * The class MapProj is dedicated to geographic projections. It can use both internal algorithms and projections provided by <a href="http://trac.osgeo.org/proj/">libproj4</a>.
 * @subsection coord_conv Coordinate conversion
 * The class MapProj takes one or two arguments describing the coordinate system of the input data and then converts back and forth with lat/long WGS84. It can be used to construct a local coordinate system, that is to say a metric grid whose origin is chosen by the user (through the lat/long parameters provided to the constructor). This is useful when working with multiple gridded coordinate system in order to get a common system that would still allow easy distances calculations.
 * @subsection dist_sec Distances
 * A few method used internally to work with custom, local grids are exposed to the user in order to easily compute distances beetwen points (using their lat/long). The algorithms can optionnaly be chosen (otherwise a default choice is used).
 * 
 *
 * @section interpol_sec Interpolations
 * @subsection interpol2d_sec Spatial interpolations
 * The class Meteo2DInterpolator receives a Digital Elevation Model (DEM) in its constructor as well as two vectors, one of MeteoData the other one of StationData. Then it allows filling 2D grid (as Grid2DObject) with spatially interpolated meteorological parameters.
 * @subsection interpol1d_sec 1D interpolations
 * The Interpol1D is dedicated to 1D interpolations.
 *
 * @section config_sec Configuration files handling
 * In order to offer a consistent interface to the user as well as make it easy to read configuration parameters, the class ConfigReader is exposed. Once constructed with a configuration file name, each key's parameter can be retrieved with a call to the templatized getValue() method.
 * 
 *
 * @section date_sec Dates handling
 * Dates should be constructed as a Date_IO object. Then, it is easy to built a new date from a julian date, from an ISO formatted date string, from a date split in fields or from a UNIX date (number of seconds since epoch). Then the various representation of the date can be retrieved, the date arithmetics can be done (for example, get the date that is 1 year, 3 months, 15 hours after 2008-12-01T11:54:00) as well as comparisons. The date printing can be controlled by keywords.
 * 
 *
 * @section exceptions_sec Exceptions
 * A few customized exceptions have been defined in IOException : these exceptions have to do with I/O, parameter parsing, argument validity, etc and consistently print usefull debuging information when thrown.
 * 
 *
 * @section misc_sec Miscellaneous
 * The IOUtils class is a static class that contains a few helper functions, in order to deal with things as diverse as units conversions, checking for file presence, equality within a given epsilon, string parsing.
 * 
 */

//Plugins overview given in IOHandler.cc

 /**
 * @page filters Filters overview
 *
 */

 /**
 * @page dev_plugins Plugins developer's guide
 *
 */

 /**
 * @page dev_filters Filters developer's guide
 *
 */

/**
 * @page examples Examples
 * Here is a simple exmaple showing how to get some meteorological data into the MeteoData and StationData vectors.
 * \code 
 * #include <iostream>
 * #include "MeteoIO.h"
 * 
 * int main(int argc, char** argv) {
 * 	(void)argc;
 * 	//provide date as ISO formatted, for example 2008-12-01T15:35:00
 * 	Date_IO d1;
 * 	std::vector<MeteoData> vecMeteo;
 * 	std::vector<StationData> vecStation;
 * 
 * 	IOHandler *raw_io = NULL;
 * 	BufferedIOHandler *io = NULL;
 * 
 * 	try {
 * 		ConfigReader cfg("io.ini");
 * 		raw_io = new IOHandler(cfg);
 * 		io = new BufferedIOHandler(*raw_io, cfg);
 * 	} catch (IOException& e){
 * 		std::cout << "Problem with IOHandler creation, cause: " << e.what() << std::endl;
 * 	}
 * 	
 * 	try {
 * 		convertString(d1,argv[1]);
 * 		io->readMeteoData(d1, vecMeteo, vecStation);
 * 	} catch (IOException& e){
 * 		std::cout << "Problem when reading data, cause: " << e.what() << std::endl;
 * 	}
 * 	
 * 	//writing some data out in order to prove that it really worked!
 * 	for (unsigned int ii=0; ii < vecMeteo.size(); ii++) {
 * 		std::cout << "---------- Station: " << (ii+1) << " / " << vecStation.size() << std::endl;
 * 		std::cout << vecStation[ii].toString() << std::endl;
 * 		std::cout << vecMeteo[ii].toString() << std::endl;
 * 	}
 * 
 * 	delete io;
 * 	delete raw_io;
 * 
 * 	return 0;
 * }
 * \endcode
 *
 * Now, we can also read a Digital Elevation Model, extract a sub set as defined by some geographical coordinates and distances and write it back to disk:
 * \code
 * #include "MeteoIO.h"
 * 
 * int main(void) {
 * 	const double lat1=46.1592, lon1=8.12993;
 * 	const double dist_x=700, dist_y=1200;
 * 	DEMObject dem;
 * 	IOHandler *raw_io = NULL;
 * 	int i,j;
 * 	
 * 	try {
 * 		ConfigReader cfg("io.ini");
 * 		raw_io = new IOHandler(cfg);
 * 	} catch (IOException& e){
 * 		std::cout << "Problem with IOHandler creation, cause: " << e.what() << std::endl;
 * 	}
 * 	raw_io->readDEM(dem);
 * 	dem.WGS84_to_grid(lat1, lon1, i,j);
 * 
 * 	const int ncols = (int)ceil(dist_x/dem.cellsize);
 * 	const int nrows = (int)ceil(dist_y/dem.cellsize);
 * 
 * 	DEMObject sub_dem(dem, i, j, ncols, nrows);
 * 	io->write2DGrid(sub_dem,"sub_dem.dem");
 * 
 * 	return 0;
 * }
 * \endcode
 * 
 * The next example shows how to compute and output spatial interpolations using previously read meteorological data and DEM.
 * \code
 * void spatial_interpolations(IOHandler& io, DEMObject& dem, std::vector<MeteoData>& vecMeteo, 
 * 		std::vector<StationData>& vecStation);
 * {
 * 	Grid2DObject    p(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
 * 	Grid2DObject hnw(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
 * 	Grid2DObject   vw(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
 * 	Grid2DObject   rh(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
 * 	Grid2DObject   ta(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
 * 	Meteo2DInterpolator mi(dem, vecMeteo, vecStation);
 * 
 * 	mi.interpolate(hnw, rh, ta, vw, p);
 * 	
 * 	std::cout << "Writing the Grids to *.2d files" << std::endl;
 * 	io->write2DGrid(ta, "output/ta.2d");
 * 	io->write2DGrid(p, "output/p.2d");
 * 	io->write2DGrid(vw, "output/vw.2d");
 * 	io->write2DGrid(nswc, "output/nswc.2d");
 * 	io->write2DGrid(rh, "output/rh.2d");
 * 	
 * }
 * \endcode
 */

#endif