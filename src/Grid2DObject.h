#ifndef __GRID2DOBJECT_H__
#define __GRID2DOBJECT_H__

#include "Array2D.h"
#include "IOExceptions.h"

/**
 * @class Grid2DObject
 * @brief A class to represent 2D Grids. Typical application as DEM or Landuse Model.
 *
 * @author Thomas Egger
 * @date   2008-12-20
 */
#ifdef _POPC_
class Grid2DObject : POPBase {
	public:
		/**
		* @brief Serialize method for POPC. Used to marshall data and send it from an object to another.
		*
		* @param buf pointer to the communication buffer
		* @param pack indicates if the data is sent or received
		*/
		virtual void Serialize(POPBuffer &buf, bool pack);
#else
class Grid2DObject{
#endif
	public:
		/**
		* @brief Default constructor.
		* Initializes all variables to 0, except lat/long which are initialized to IOUtils::nodata
		*/
		Grid2DObject();
		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize);

		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const Array2D<double>& grid2D_in);

		/**
		* @brief constructs an object as a subset of another grid object
		* @param _grid2Dobj (const Grid2DObject&) initial grid object
		* @param _nx (const unsigned int) starting column of the subset
		* @param _ny (const unsigned int) starting row of the subset
		* @param _ncols (const unsigned int) number of columns of the subset
		* @param _nrows (const unsigned int) number of rows of the subset
		*/
		Grid2DObject(const Grid2DObject& _grid2Dobj,
				   const unsigned int& _nx, const unsigned int& _ny, //Point in the plane
				   const unsigned int& _ncols, const unsigned int& _nrows); //dimensions of the sub-plane

		/**
		* @brief Set all variables in one go.
		* @param ncols (unsigned int) number of colums in the grid2D
		* @param nrows (unsigned int) number of rows in the grid2D
		* @param xllcorner (double) x-coordinate of lower left corner
		* @param yllcorner (double) y-coordinate of lower left corner
		* @param latitude (double) decimal latitude
		* @param longitude (double) decimal longitude
		* @param cellsize (double) value for cellsize in grid2D
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize);
		/**
		* @brief Set all variables in one go. Notably the member grid2D of type 
		* Array2D\<double\> will be destroyed and recreated to size ncols x nrows.
		* @param ncols (unsigned int) number of colums in the grid2D
		* @param nrows (unsigned int) number of rows in the grid2D
		* @param xllcorner (double) x-coordinate of lower left corner
		* @param yllcorner (double) y-coordinate of lower left corner
		* @param latitude (double) decimal latitude
		* @param longitude (double) decimal longitude
		* @param cellsize (double) value for cellsize in grid2D
		* @param grid2D_in (CArray\<double\>&) grid to be copied by value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const Array2D<double>& grid2D_in); //TODO: const CArray would be better...
		
		/**
		* @brief check if the current Grid2DObject has the same geolocalization attributes
		* as another Grid2DObject. The grid coordinates (xllcorner & yllcorner) are NOT
		* checked as these might be tweaked for convenience (like between input grid and local grid)
		* @param target (Grid2DObject) grid to compare to
		* @return (bool) true if same geolocalization
		*/
		bool isSameGeolocalization(const Grid2DObject& target);

		Array2D<double> grid2D;
		unsigned int ncols, nrows;
		double xllcorner, yllcorner, cellsize;
		double latitude, longitude;

 protected:
		void setValues(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude, const double& cellsize);

		void checkCoordinates();
};
#include "DEMObject.h"

#endif
