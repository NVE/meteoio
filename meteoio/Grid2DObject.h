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

#include <meteoio/Coords.h>
#include <meteoio/Array2D.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <iostream>
#include <cmath>

namespace mio {

/**
 * @class Grid2DObject
 * @brief A class to represent 2D Grids. Typical application as DEM or Landuse Model.
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2008-12-20
 */

#ifdef _POPC_
class Grid2DObjetctDummy {}; //HACK for POPC

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
class Grid2DObject {
#endif
	public:
		///structure to contain the grid coordinates of a point in a 2D grid
		typedef struct GRID_POINT_2D {
			unsigned int ix; ///<grid index along X
			unsigned int iy; ///<grid index along Y
		} grid_point_2d;

		Grid2DObject& operator=(const Grid2DObject&); ///<Assignement operator
		friend std::ostream& operator<<(std::ostream& os, const Grid2DObject& grid);

		/**
		* @brief Default constructor.
		* Initializes all variables to 0, except lat/long which are initialized to IOUtils::nodata
		*/
		Grid2DObject();
		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
		             const double& cellsize, const Coords& _llcorner);

		Grid2DObject(const unsigned int& ncols, const unsigned int& nrows,
		             const double& cellsize, const Coords& _llcorner, const Array2D<double>& grid2D_in);

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
		* @brief Compute the positional parameters that are not already known
		* This means that the Coords::point object that is given either contains geographic coordinates or
		* grid indices. This method will calculate the missing ones (so that (i,j) match with (lat,lon)
		* and (east,north)). If the given point had a "NULL" projection, it will be set to the grid's.
		* @param point coordinate to convert
		* @return false if the given point was invalid or outside the grid (sets (i,j) to closest values within the grid)
		*/
		bool gridify(Coords& point) const;

		/**
		* @brief Compute the positional parameters that are not already known
		* This means that the Coords::point object that is given either contains geographic coordinates or
		* grid indices. This method will calculate the missing ones (so that (i,j) match with (lat,lon)
		* and (east,north)). Any point that is either invalid or outside the grid is removed from the vector.
		* If the given point had a "NULL" projection, it will be set to the grid's.
		* @param vec_points vector containing the coordinates to convert
		* @return false if invalid or external points had to be removed
		*/
		bool gridify(std::vector<Coords>& vec_points) const;

		/**
		* @brief Set all variables in one go.
		* @param ncols (unsigned int) number of colums in the grid2D
		* @param nrows (unsigned int) number of rows in the grid2D
		* @param cellsize (double) value for cellsize in grid2D
		* @param _llcorner lower left corner point
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows,
		         const double& cellsize, const Coords& _llcorner);
		/**
		* @brief Set all variables in one go. Notably the member grid2D of type
		* Array2D\<double\> will be destroyed and recreated to size ncols x nrows.
		* @param ncols (unsigned int) number of colums in the grid2D
		* @param nrows (unsigned int) number of rows in the grid2D
		* @param cellsize (double) value for cellsize in grid2D
		* @param _llcorner lower left corner point
		* @param grid2D_in (CArray\<double\>&) grid to be copied by value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows,
		         const double& cellsize, const Coords& _llcorner, const Array2D<double>& grid2D_in); //TODO: const CArray would be better...

		/**
		* @brief check if the current Grid2DObject has the same geolocalization attributes
		* as another Grid2DObject (as well as same cells). The grid coordinates (xllcorner & yllcorner) are NOT
		* checked as these might be tweaked for convenience (like between input grid and local grid)
		* @param target (Grid2DObject) grid to compare to
		* @return (bool) true if same geolocalization
		*/
		bool isSameGeolocalization(const Grid2DObject& target) const;

		Array2D<double> grid2D; ///<the grid itself (simple 2D table containing the values for each point)
		unsigned int ncols; ///<number of columns in the grid
		unsigned int nrows; ///<number of rows in the grid
		double cellsize; ///<dimension in meters of a cell (considered to be square)
		Coords llcorner; ///<lower left corner of the grid

		/**
		* @brief Partitional algorithm to classify each point of the grid.
		* The classification is given by a list of growing thresholds, the 'clusters' are then a simple
		* range of values. Each cluster comes with an 'id' that replaces the values of the points.
		*
		*
		* @param thresholds (const std::vector<double>&) ordered list of thresholds representing a scale of values. Each level of this scale defines a cluster
		* @param ids (const std::vector<double>&) clusters Ids to be used. clustersId.size()=thresholds.size()+1
		* @return (bool) true if clusturization was succesfull
		*/
		bool clusterization(const std::vector<double>& thresholds, const std::vector<double>& ids);

 protected:
		void setValues(const unsigned int& ncols, const unsigned int& nrows,
		               const double& cellsize, const Coords& _llcorner);
		void setValues(const unsigned int& _ncols, const unsigned int& _nrows,
		               const double& _cellsize);

		/**
		* @brief Converts WGS84 coordinates into grid coordinates (i,j)
		* @param point coordinate to convert
		* @return false if the given point was outside the grid (sets (i,j) to closest values within the grid)
		*/
		bool WGS84_to_grid(Coords& point) const;

		/**
		* @brief Converts grid coordinates (i,j) into WGS84 coordinates
		* @param point coordinate to convert
		* @return false if the given point was invalid (outside the grid or nodata and if possible sets (i,j) to closest values within the grid)
		*/
		bool grid_to_WGS84(Coords& point) const;
};
} //end namespace

#endif
