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
 *        a derived class is to be created that holds the specific implementation of the purely virtual member funtions. So far 
 *        the following children have been implemented:
 * - class A3DIO is able to interface with ASCII files
 * - class BoschungIO is able to interface with the specific IO system of the company Boschung
 * - class IOHandler is a wrapper class that is able to deal with all above implementations of the IOInterface abstract base class.
 *   Which implementation is used for a specific member function can be configured in a key/value file.
 *
 * @author Thomas Egger
 * @date   2009-01-08
 */
class IOInterface : public PluginObject {
	public:

		IOInterface(void (*delObj)(void*));

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
		* @brief Get the grid size of a 2D grid in cells. The underlying mechanism may vary 
		* (e.g. A3DIO reads the header of the DEM file to retrieve the information)
		*    
		* Example Usage:
		* @code
		* int nx, ny;
		* IOHandler io1("io.ini");
		* io1.get2DGridSize(nx, ny);
		* @endcode
		* @param nx Number of columns of the 2D grid
		* @param ny Number of rows of the 2D grid
		*/
		virtual void get2DGridSize(int& nx, int& ny) = 0;

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
		* @brief See MeteoData::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation).
		*/
		//virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo) = 0;

		/**
		* @brief Fill vector<MeteoData> and vector<StationData> objects with multiple datasets 
		* corresponding to the time indicated by the Date_IO object.
		* Matching rule: Find first data set for every station which has an event time (measurement time) that is greater 
		* (newer) or equal to the time represented by the Date_IO object parameter. The vector<StationData> object holds multiple 
		* StationData objects representing meta information about the meteo stations that recorded the meteo data.
		*    
		* Example Usage:
		* @code
		* vector<MeteoData> vecMeteo;      //empty vector
		* vector<StationData> vecStation;  //empty vector
		* Date_IO d1(2008,06,21,11,00);       //21.6.2008 11:00
		* IOHandler io1("io.ini");
		* io1.readMeteoData(d1, vecMeteo, vecStation);
		* @endcode
		* @param date_in     A Date_IO object representing the approximate date/time for the sought MeteoData objects
		* @param vecMeteo    A vector of MeteoData objects to be filled with data
		* @param vecStation  A vector of StationData objects to be filled with data
		*/
		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  unsigned int stationindex=IOUtils::npos) = 0;

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

		virtual void readSpecialPoints(CSpecialPTSArray& pts) = 0;

		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name="") = 0;
};

#endif
