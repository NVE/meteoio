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
#ifndef __IOINTERFACE_H__
#define __IOINTERFACE_H__

#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"
#include "DEMObject.h"
#include "Date_IO.h"
#include "DynamicLibrary.h"
#include "Array.h"
#include "Array2D.h"

#include <vector>

/**
 * @class IOInterface
 * @brief An abstract class representing the IO Layer of the software Alpine3D. For each type of IO (File, DB, Webservice, etc)
 * a derived class is to be created that holds the specific implementation of the purely virtual member funtions. 
 * The IOHandler class is a wrapper class that is able to deal with all above implementations of the IOInterface abstract base class.
 *
 * @author Thomas Egger
 * @date   2009-01-08
 */
class IOInterface : public PluginObject {
	public:

		IOInterface(void (*delObj)(void*));
		virtual ~IOInterface();

		/**
		* @brief A virtual copy constructor used for cloning objects
		*    
		* Example Usage:
		* @code
		* IOInterface *io1,*io2;
		* io1 = new IOHandler("io.ini");
		* io2 = io1->clone(); 
		* @endcode
		* @return A pointer to the cloned object.
		*/
		//virtual IOInterface* clone() const = 0;

		/**
		* @brief A generic function for parsing 2D grids into a Grid2DObject. The string parameter shall be used for addressing the 
		* specific 2D grid to be parsed into the Grid2DObject.
		* @param grid_out A Grid2DObject instance 
		* @param parameter A std::string representing some information for the function on what grid to retrieve
		*/ 
		virtual void read2DGrid(Grid2DObject& grid_out, const string& parameter="") = 0;

		/**
		* @brief Parse the DEM (Digital Elevation Model) into the Grid2DObject
		*    
		* Example Usage:
		* @code
		* Grid2DObject dem;
		* IOHandler io1("io.ini");
		* io1.readDEM(dem);
		* @endcode
		* @param dem_out A Grid2DObject that holds the DEM
		*/
		virtual void readDEM(DEMObject& dem_out) = 0;

		/**
		* @brief Parse the landuse model into the Grid2DObject
		*    
		* Example Usage:
		* @code
		* Grid2DObject landuse;
		* IOHandler io1("io.ini");
		* io1.readLanduse(landuse);
		* @endcode
		* @param landuse_out A Grid2DObject that holds the landuse model
		*/
		virtual void readLanduse(Grid2DObject& landuse_out) = 0;


		/**
		* @brief Fill vecMeteo and vecStation with a time series of objects  
		* corresponding to the interval indicated by dateStart and dateEnd.
		*
		* Matching rules:
		* - if dateStart and dateEnd are the same: return exact match for date
		* - if dateStart > dateEnd: return first data set with date > dateStart
		* - read in all data starting with dateStart until dateEnd
		* - if there is no data at all then the vectors will be empty, no exception will be thrown
		*    
		* Example Usage:
		* @code
		* vector< vector<MeteoData> > vecMeteo;      //empty vector
		* vector< vector<StationData> > vecStation;  //empty vector
		* Date_IO d1(2008,06,21,11,00);       //21.6.2008 11:00
		* Date_IO d2(2008,07,21,11,00);       //21.7.2008 11:00
		* IOHandler io1("io.ini");
		* io1.readMeteoData(d1, d2, vecMeteo, vecStation);
		* @endcode
		* @param dateStart   A Date_IO object representing the beginning of an interval (inclusive)
		* @param dateEnd     A Date_IO object representing the end of an interval (inclusive)
		* @param vecMeteo    A vector of vector<MeteoData> objects to be filled with data
		* @param vecStation  A vector of vector<StationData> objects to be filled with data
		* @param stationindex (optional) update only the station given by this index
		*/
		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos) = 0;

		/**
		* @brief Parse the assimilation data into a Grid2DObject for a certain date represented by the Date_IO object
		*    
		* Example Usage:
		* @code
		* Grid2DObject adata;
		* Date_IO d1(2008,06,21,11,00);       //21.6.2008 11:00
		* IOHandler io1("io.ini");
		* io1.readAssimilationData(d1, adata);
		* @endcode
		* @param date_in A Date_IO object representing the date of the assimilation data
		* @param da_out  A Grid2DObject that holds the assimilation data for every grid point
		*/
		virtual void readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out) = 0;

		/**
		* @brief Read a list of points by their grid coordinates
		* This allows for example to get a list of points where to produce more detailed outputs.
		* @param pts (CSpecialPTSArray) A vector of points coordinates
		*/
		virtual void readSpecialPoints(CSpecialPTSArray& pts) = 0;

		/**
		* @brief Write a Grid2DObject
		* @param grid_in (Grid2DObject) The grid to write
		* @param name (string) Identifier usefull for the output plugin (it could become part of a file name, a db table, etc)
		*/
		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name="") = 0;
};

#endif
