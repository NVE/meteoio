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
		void Serialize(paroc_buffer &buf, bool pack);
#else
class Grid2DObject{
#endif  
	public:
		/**
		* Default constructor.
		* Initializes all variables to 0, except nodata, which is initialized to -9999.0
		*/
		Grid2DObject();
		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const double& nodata);

		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const double& nodata,
		       CArray2D<double>& grid2D_in);

		/**
		* @brief Set all variables in one go.
		* @param ncols (unsigned int&) number of colums in the grid2D
		* @param nrows (unsigned int&) number of rows in the grid2D
		* @param xllcorner (double&) x-coordinate of lower left corner
		* @param yllcorner (double&) y-coordinate of lower left corner
		* @param latitude (double&) decimal latitude
		* @param longitude (double&) decimal longitude
		* @param cellsize (double&) value for cellsize in grid2D
		* @param nodata (double&) value representing a NODATA value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const double& nodata);
		/**
		* @brief Set all variables in one go. Notably the member grid2D of type CArray2D<double> 
		* will be destroyed and recreated to size ncols x nrows.
		*
		* @param ncols (unsigned int&) number of colums in the grid2D
		* @param nrows (unsigned int&) number of rows in the grid2D
		* @param xllcorner (double&) x-coordinate of lower left corner
		* @param yllcorner (double&) y-coordinate of lower left corner
		* @param latitude (double&) decimal latitude
		* @param longitude (double&) decimal longitude
		* @param cellsize (double&) value for cellsize in grid2D
		* @param nodata (double&) value representing a NODATA value
		* @param grid2D_in (CArray<double>&) grid to be copied by value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const double& nodata,
			CArray2D<double>& grid2D_in);


		/**
		* @brief Serialize method for POPC. Used to marshall data and send it from an object to another.
		*
		* @param buf pointer to the communication buffer
		* @param pack indicates if the data is sent or received
		*/  
		CArray2D<double> grid2D;
		unsigned int ncols, nrows;
		double xllcorner, yllcorner, cellsize, nodata;
		double latitude, longitude;
};

#endif
