#ifndef __IOINTERFACE_H__
#define __IOINTERFACE_H__

#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"
#include "Date_IO.h"
#include "DynamicLibrary.h"
#include "Array.h"
#include "Array2D.h"

#include <vector>

/**
 * @class IOInterface
 * @brief An abstract class representing the IO Layer of the software Alpine3D. For each type of IO (File, DB, Webservice, etc)
 * a derived class is to be created that holds the specific implementation of the purely virtual member funtions. So far 
 * the following children have been implemented (by keyword for the io.ini key/value config file):
 * - A3D for reading original Alpine3D meteo files (no extra requirements)
 * - BOSCHUNG for reading Boshung xml meteo files (requires libxml)
 * - IMIS for reading meteo data out of the IMIS database (requires Oracle's OCCI library)
 * - GEOTOP for reading original GeoTop meteo files (no extra requirements)
 * - GSN for reading meteo data out of the Global Sensor Network web service interface (requires GSoap)
 * - ARC for reading ESRI/ARC DEM files (no extra requirements)
 * - GRASS for reading Grass DEM files (no extra requirements)
 * 
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
		virtual void readDEM(Grid2DObject& dem_out) = 0;

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
