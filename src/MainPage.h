/***********************************************************************************/
/*  Copyright 2009-2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS */
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
#ifndef __MAINPAGE_H__
#define __MAINPAGE_H__

namespace mio {
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
 * -# End User documentation
 *    -# \subpage general "General concepts"
 *    -# \subpage plugins "Available plugins" and usage
 *    -# \subpage coords "Available coordinate systems" and usage
 *    -# \subpage filters "Available filters" and usage
 *    -# \subpage interpol2d "Available spatial interpolations" and usage
 * -# Programing using MeteoIO
 *    -# \subpage quick_overview "Quick overview" of the functionnality provided by MeteoIO
 *    -# \subpage examples "Usage examples"
 * -# Expanding MeteoIO
 *    -# How to \subpage dev_plugins "write a Plugin"
 *    -# How to \subpage dev_filters "write a Filter"
 */

 /**
 * @page general General concepts
 * Since MeteoIO is a library, you, as an end user, will have a limited direct exposure to it: the library is called by the program that you are using, not directly by yourself. You will basically have to set some parameters in a configuration file that defines how MeteoIO has to behave. This configuration file is often named "io.ini" and follows the INI file format standard (see http://en.wikipedia.org/wiki/INI_file). In order to understand how this file is structured, let us first have a look at the general structure of MeteoIO and afterward the structure of this configuration file and where to find the available configuration parameters.
 *
 * @section MeteoIO_structure General MeteoIO structure
 * MeteoIO can be seen as a set of modules that is focused on the handling of input/output operations (including data preparation) for numerical simulations in the realm of earth sciences. On the visible side, it offers the following modules:
 * - a set of plugins for accessing the data (for example, a plugin might be responsible for fetching the raw data from a given database)
 * - a set of filters and processing elements for applying transformations to the data (for example, a filter might remove all data that is out of range)
 * - a set of spatial interpolation algorithms (for example, such an algorithm might perform Inverse DIstance Weighting for filling a grid with spatially interpolated data)
 *
 * Moreover, a few assumptions are made about the data that you are using: each data point has to be associated with a geographic location (defined by some sort of coordinates) and very often you will also need to provide a Digital Elevation Model. Therefore, you will also notice a few extra modules that come to play on the visible side:
 * - a module to deal with Digital Elevation Models. Such module will, for example, interpret a grid of data as a grid of elevations and compute a grid of slopes.
 * - a module to deal with coordinate systems. Such module will require you to define which coordinate system are your data in and transparently handle potential coordinate conversions in the program that you are using.
 * - a module to deal with... configuration files. The program that you are using might be using this module for other configuration files.
 *
 * @section Config Configuration file
 * @subsection Config_syntax Configuration file syntax
 * The configuration inputs/outputs file is divided in sections. Each section name is enclosed in brackets.
 * @code
 * [My_section]
 * bla bla bla
 * @endcode
 *
 * Within each section, you might have comments and/or key/value pairs. The comments start with a "#" or a ";" sign and run until the end of the line. A whole line might be commented out, or only a fraction of it. A key/value pair is a keyword, followed by a "=" sign, followed by the value to associate with this key.
 * @code
 * #This is a commented out line
 * PI = 3.14159   #key/value pair and a comment
 * @endcode
 *
 * A valid value can be an integer, a float, or string, a list of keywords, a mixed list of keywords and numbers...
 * @code
 * TA::algorithms = IDW_LAPSE CST_LAPSE
 * @endcode
 *
 * @subsection Config_structure Configuration file structure
 * MeteoIO imposes a minimum structure to the configuration file: It must contain the [General], [Input] and [Output] sections. If any filter is to be used, a [Filters] section has to be present and if any spatial interpolation is to be used, an [Interpolations2D] section has to be present. A minimal set of keys has to be there, an potentially a number of optional keys. Moreover, the program that you are using might also impose you some specific keys or sections.
 * The keys and their location in the configuration file (ie: to which section they belong) depends on the module that is actually using them. The optional keys depend on the specific options handled by each specific module (or plugin, or algorithm). Therefore, we can draw the following skeleton:
 * @code
 * [General]
 * PLUGINPATH	= ../../lib	#optional, if not provided the plugins will be search in the normal library path
 *
 * [Input]
 * COORDSYS	= CH1903	#mandatory: which coordinate system is used for the geographic coordinates
 * COORDPARAM	= -999		#extra arguments for the chosen coordinate system (often, none)
 *
 * DEM		= ARC		#plugin to use for reading DEM information
 * #this might be followed by any number of arguments that are specific to this plugin
 *
 * METEO	= A3D		#plugin to use for reading meteorological data
 * #this might be followed by any number of arguments that are specific to this plugin
 *
 * [Output]
 * COORDSYS	= CH1903
 *
 * GRID2D	= ARC		#plugin to use for writing 2D grids
 *
 * [Filters]
 * TA::filter1	= min_max	#first filter to use on the parameter TA
 * TA::arg1	= 230 330	#arguments for this first filter
 * TA::filter2	= rate		#second filter to use (in chronological order)
 * TA::arg2	= 0.01		#arguments for this second filter
 * #add any extra filter that you want to use. They will be applied serially
 *
 * [Interpolations2D]
 * TA::algorithms = IDW_LAPSE CST_LAPSE	#list of algorithms to consider for use for spatially interpolating parameter TA
 * TA::cst_lapse = -0.008 		#parameter for a specific interpolation algorithm for parameter TA
 *
 * @endcode
 *
 * @subsection Finding_docs Where to find the proper documentation
 * As can be seen from the previous example, each plugin, each filter or each interpolation algorithm might have its own parameters. Therefore, this is the documentation of each specific plugin/filter/algorithm that has to be used in order to figure out what can be configured when it (see the next sections in the welcome page).
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
 * The class Coords is dedicated to geographic projections. It can use both internal algorithms and projections provided by <a href="http://trac.osgeo.org/proj/">libproj4</a>.
 * @subsection coord_conv Coordinate conversion
 * The class Coords takes one or two arguments describing the coordinate system of the input data and then converts back and forth with lat/long WGS84. It can be used to construct a local coordinate system, that is to say a metric grid whose origin is chosen by the user (through the lat/long parameters provided to the constructor). This is useful when working with multiple gridded coordinate system in order to get a common system that would still allow easy distances calculations. See the supported \ref Coordinate_types "projections".
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
 * Dates should be constructed as a Date object. Then, it is easy to built a new date from a julian date, from an ISO formatted date string, from a date split in fields or from a UNIX date (number of seconds since epoch). Then the various representation of the date can be retrieved, the date arithmetics can be done (for example, get the date that is 1 year, 3 months, 15 hours after 2008-12-01T11:54:00) as well as comparisons. The date printing can be controlled by keywords.
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
//Filters overview given in FilterAlgorithms.cc
//Plugin development given in IOInterface.h

 /**
 * @page dev_filters Filter developer's guide
 *
 * In order to add a new filter to the already existing set of filters the developer only needs to edit
 * the class FilterAlgorithms. It is important to understand that the filter operate on a "per parameter" basis.
 * This means that a filter might be executed for the parameter TA and another one for the parameter HNW, so the filter
 * algorithm only has to deal with a generic filtering method based on double values.
 *
 * To implement a new filter two steps are necessary:
 *
 * -# Registering the filter within the function FilterAlgorithms::initStaticData(), by simply adding a line that
 *    resembles:
 *    @code
 *    filterMap["filtername"]=FilterProperties(true, (unsigned int)1, Date(0.0), &FilterAlgorithms::ImplementedFilter);
 *    @endcode
 *    Where the filtername can be freely chosen (although it has to be unique) and the parameters to FilterProperties
 *    are (in this order) a boolean stating whether this filter is only a check or whether it changes data,
 *    an unsigned integer stating how many points are minimally needed for this filter to be executed, a Date object
 *    marking the minimal period that data supplied to this filter needs to span and finally a function pointer
 *    pointing towards the actual implementation of the filter.
 * -# Implementation of the filter, a static function which always has the following interface and returns a boolean:
 *    @code
 *    bool FilterAlgorithms::FilterName(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
 *			  const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
 *			  const unsigned int& paramindex,
 *			  std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
 *    @endcode
 *    vecM and vecS represent the raw data as it is read through the readMeteoData functions of the MeteoIO plugins.
 *    the unsigned integer pos is the index of either the elements within vecM and vecS that represents the data
 *    data for the given Date object called date, or if there is no exact match possible then pos is the index of
 *    the first tuple with a date greater (that is newer) than the date requested (in this case a resampling will have
 *    to interpolate the data for the date requested). paramindex is the MeteoData parameter that this filter shall be
 *    run upon (see MeteoData for the enumeration of parameters - e.g. TA(=0), HNW(=6), ISWR(=1)). The two vectors
 *    vecFilteredM and vecFilteredS correspond and they hold the filtered data. Typically they hold one or two elements
 *    depending on whether resampling is required or not. In case resampling is not required they only hold one value
 *    each representing the exact date requested. Changes on the vecFilteredM and vecFilteredS vectors will be
 *    propagated back to the BufferedIOHandler and thus to the user. If the filter has to alter data it may do so only
 *    on these two vectors. Typically filters that perform checks only (e.g. check whether values are within certain
 *    bounds), need not look at anything but the data in these vectors, whereas filters that have more complicated
 *    schemes of operation (like accumulation or resampling) might take into account the "raw" data as accessible
 *    through vecM and vecS. Helper functions that parse the arguments to the filters through a ConfigReader obejct
 *    are available.
 *
 *
 * Here an example implementation of the MaximumFilter, which checks whether a value is greater then an argument
 * supplied to the filter and if so changes the value either to IOUtils::nodata (normal operation) or to the
 * maximum value supplied in the argument (soft mode of operation). An example section in the io.ini file supplied
 * to the ConfigReader could look like this:
 * @code
 * [Filters]
 * TA::filter1 = max
 * TA::arg1    = soft 280
 * @endcode
 * Which has the following interpretation: Apply filter max (max-value-filter) to the parameter TA (air temperature)
 * in case that a value is greater than 280 degrees Kelvin change that value to 280.
 *
 * Now for the actual implementation of the filter:
 * @code
 * bool FilterAlgorithms::MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
 *	         const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
 *		    const unsigned int& paramindex,
 *		    std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
 * {
 *      (void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
 *      //parse arguments and check whether they are valid
 *      bool isSoft = false;
 *      std::vector<double> vecArgs;
 *      parseFilterArguments("max", _vecArgs, 1, 1, isSoft, vecArgs);
 *
 *      //Run actual MaxValue filter over all relevant meteo data
 *      for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
 *      double& value = vecFilteredM[ii].param(paramindex);
 *      if (value == IOUtils::nodata) continue;
 *        if (value>vecArgs[0]){
 *        	if (isSoft) value=vecArgs[0];
 *        	else value=IOUtils::nodata;
 *        }
 *      }
 *
 *      return true;
 * }
 * @endcode
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
 * 	Date d1;
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
 * 		std::cout << vecStation[ii] << std::endl;
 * 		std::cout << vecMeteo[ii] << std::endl;
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
 * 	const double dist_x=700, dist_y=1200;
 * 	DEMObject dem;
 * 	IOHandler *raw_io = NULL;
 *	ConfigReader *cfg = NULL;
 * 	int i,j;
 *
 * 	try {
 * 		cfg = new ConfigReader("io.ini");
 * 		raw_io = new IOHandler(cfg);
 * 	} catch (IOException& e){
 * 		std::cout << "Problem with IOHandler creation, cause: " << e.what() << std::endl;
 * 	}
 * 	raw_io->readDEM(dem);
 *	Coords point(*cfg);
 *	point.setLatLon(46.1592, 8.12993);
 *	dem.WGS84_to_grid(point, i,j);
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
 * 	Grid2DObject ta;
 * 	Meteo2DInterpolator mi(dem, vecMeteo, vecStation);
 *
 * 	mi.interpolate(MeteoData::TA, ta);
 *
 * 	std::cout << "Writing the Grids to files" << std::endl;
 * 	io->write2DGrid(ta, "ta_2d.asc");
 *
 * }
 * \endcode
 */

} //end namespace mio
#endif
