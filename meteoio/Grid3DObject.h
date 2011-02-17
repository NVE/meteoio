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
#ifndef __GRID3DOBJECT_H__
#define __GRID3DOBJECT_H__

#include <meteoio/Coords.h>
#include <meteoio/Array3D.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/Grid2DObject.h>

#include <cmath>
#include <iostream>

namespace mio {

/**
 * @class Grid3DObject
 * @brief A class to represent 3D Grids. Typical application: wind field
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2009-07-20
 */

#ifdef _POPC_
class Grid3DObjetctDummy {}; //HACK for POPC

class Grid3DObject : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Grid3DObject{
#endif  
	public:
		typedef struct GRID_POINT_3D { //TODO: this potentially conflicts with the definition in Grid2DObject
			unsigned int ix; ///<grid index along X
			unsigned int iy; ///<grid index along Y
			unsigned int iz; ///<grid index along Z
		} grid_point_3d;

		Grid3DObject& operator=(const Grid3DObject&); ///<Assignement operator
		friend std::ostream& operator<<(std::ostream& os, const Grid3DObject& grid);

		/**
		* Default constructor.
		* Initializes all variables to 0, except nodata, which is initialized to -9999.0
		*/
		Grid3DObject();

		/**
		* A constructor that can be used to create a Grid3DObject that is contained in the
		* one passed as _grid3Dobj argument. The resulting Grid3DObject is a by value copy of
		* a subspace of the space spanned by the _grid3Dobj
		*/
		Grid3DObject(const Grid3DObject& _grid3Dobj,
				   const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
				   const unsigned int& _nwidth, const unsigned int& _nheight, const unsigned int& _ndepth);

		Grid3DObject(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
				const double& cellsize, const Coords& _llcorner);

		Grid3DObject(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& cellsize, const Coords& _llcorner, const Array3D<double>& grid3D);

		/**
		* @brief Set all variables in one go.
		* @param ncols (unsigned int&) number of colums in the grid3D (1st dimension)
		* @param nrows (unsigned int&) number of rows in the grid3D (2nd dimension)
		* @param depth (unsigned int&) number of depth in the grid3D (3rd dimension)
		* @param cellsize (double&) value for cellsize in grid3D
		* @param _llcorner lower left corner coordinates
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& depth,
			const double& cellsize, const Coords& _llcorner);
		/**
		* @brief Set all variables in one go. Notably the member grid3D of type Array3D<double> 
		* will be destroyed and recreated to size ncols x nrows.
		*
		* @param ncols (unsigned int&) number of colums in the grid3D
		* @param nrows (unsigned int&) number of rows in the grid3D
		* @param ndepth (unsigned int&) number of depth in the grid3D (3rd dimension)
		* @param cellsize (double&) value for cellsize in grid3D
		* @param _llcorner lower left corner coordinates
		* @param grid3D_in (Array\<double\>&) grid to be copied by value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& cellsize, const Coords& _llcorner, const Array3D<double>& grid3D_in);

		/**
		* @brief Compute the positional parameters that are not already known
		* This means that the Coords::point object that is given either contains geographic coordinates or 
		* grid indices. This method will calculate the missing ones (so that (i,j,k) match with (lat,lon,alt)
		* and (east,north,alt).
		* @param point coordinate to convert
		* @return false if the given point was invalid or outside the grid (sets (i,j) to closest values within the grid)
		*/
		bool gridify(Coords& point) const;

		/**
		* @brief Compute the positional parameters that are not already known
		* This means that the Coords::point object that is given either contains geographic coordinates or 
		* grid indices. This method will calculate the missing ones (so that (i,j) match with (lat,lon)
		* and (east,north)). Any point that is either invalid or outside the grid is removed from the vector.
		* @param vec_points vector containing the coordinates to convert
		* @return false if invalid or external points had to be removed
		*/
		bool gridify(std::vector<Coords>& vec_points) const;

		/**
		* @brief check if the current Grid3DObject has the same geolocalization attributes
		* as another Grid3DObject. The grid coordinates (xllcorner & yllcorner) are NOT
		* checked as these might be tweaked for convenience (like between input grid and local grid)
		* @param target (Grid3DObject) grid to compare to
		* @return (bool) true if same geolocalization
		*/
		bool isSameGeolocalization(const Grid3DObject& target);
		
		void extractLayer(const unsigned int& z, Grid2DObject& layer);

		Array3D<double> grid3D;
		unsigned int ncols, nrows, ndepth;
		double cellsize;
		Coords llcorner;
		//NOTE: the altitude is understood as above sea level,
		//that is we curently don't support altitude as above the local ground
		//std::vector<double> thickness;

 protected:
		void setValues(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& cellsize);
		void setValues(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& cellsize, const Coords& _llcorner);

		/**
		* @brief Converts WGS84 coordinates into grid coordinates (i,j)
		* @param point coordinate to convert
		* @return false if the given point was outside the grid (sets (i,j) to closest values within the grid)
		*/
		bool WGS84_to_grid(Coords point) const;

		/**
		* @brief Converts grid coordinates (i,j) into WGS84 coordinates
		* @param point coordinate to convert
		* @return false if the given point was invalid (outside the grid or nodata and if possible sets (i,j) to closest values within the grid)
		*/
		bool grid_to_WGS84(Coords& point) const;
};
} //end namespace

#endif
