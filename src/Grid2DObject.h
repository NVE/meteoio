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
#ifndef __GRID2DOBJECT_H__
#define __GRID2DOBJECT_H__

#include "Array2D.h"
#include "IOExceptions.h"
#include "MapProj.h"

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
			const double& cellsize/*, const MapProj& proj=Grid2DObject::NULL_proj*/);

		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const Array2D<double>& grid2D_in/*,
			const MapProj& proj=Grid2DObject::NULL_proj*/);

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
		* @brief Converts WGS84 coordinates into grid coordinates (i,j)
		* If any coordinate is outside of the grid, the matching coordinate is set to (unsigned)-1
		* @note This computation is only precise to ~1 meter (in order to be faster). Please use 
		* IOUtils::VincentyDistance for up to .5 mm precision.
		* @param _latitude point latitude
		* @param _longitude point longitude
		* @param i matching X coordinate in the grid
		* @param j matching Y coordinate in the grid
		* @return EXIT_SUCESS or EXIT_FAILURE if the given point was outside the grid (sets (i,j) to closest values within the grid)
		*/
		int WGS84_to_grid(const double& _latitude, const double& _longitude, unsigned int& i, unsigned int& j);

		/**
		* @brief Converts grid coordinates (i,j) into WGS84 coordinates
		* @note This computation uses IOUtils::VincentyDistance and therefore is up to .5 mm precise.
		* @param _latitude point latitude
		* @param _longitude point longitude
		* @param i matching X coordinate in the grid
		* @param j matching Y coordinate in the grid
		*/
		void grid_to_WGS84(const unsigned int& i, const unsigned int& j, double& _latitude, double& _longitude);

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
			const double& cellsize/*,
			const MapProj& proj=Grid2DObject::NULL_proj*/);
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
			const double& cellsize, const Array2D<double>& grid2D_in/*,
			const MapProj& proj=Grid2DObject::NULL_proj*/); //TODO: const CArray would be better...
		
		/**
		* @brief check if the current Grid2DObject has the same geolocalization attributes
		* as another Grid2DObject (as well as same cells). The grid coordinates (xllcorner & yllcorner) are NOT
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
		static const MapProj NULL_proj;

		void setValues(const unsigned int& ncols, const unsigned int& nrows,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude, const double& cellsize/*,
			const MapProj& proj=Grid2DObject::NULL_proj*/);

		void checkCoordinates(const MapProj& proj);
};

#endif
