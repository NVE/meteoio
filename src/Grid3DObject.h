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

#include "Array3D.h"
#include "IOExceptions.h"
#include "MapProj.h"

/**
 * @class Grid3DObject
 * @brief A class to represent 3D Grids. Typical application: wind field
 *
 * @author Thomas Egger
 * @date   2009-07-20
 */
#ifdef _POPC_
class Grid3DObject : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Grid3DObject{
#endif  
	public:
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
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize);

		Grid3DObject(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const Array3D<double>& grid3D);

		/**
		* @brief Set all variables in one go.
		* @param ncols (unsigned int&) number of colums in the grid3D (1st dimension)
		* @param nrows (unsigned int&) number of rows in the grid3D (2nd dimension)
		* @param depth (unsigned int&) number of depth in the grid3D (3rd dimension)
		* @param xllcorner (double&) x-coordinate of lower left corner
		* @param yllcorner (double&) y-coordinate of lower left corner
		* @param latitude (double&) decimal latitude
		* @param longitude (double&) decimal longitude
		* @param cellsize (double&) value for cellsize in grid3D
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& depth,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize);
		/**
		* @brief Set all variables in one go. Notably the member grid3D of type Array3D<double> 
		* will be destroyed and recreated to size ncols x nrows.
		*
		* @param ncols (unsigned int&) number of colums in the grid3D
		* @param nrows (unsigned int&) number of rows in the grid3D
		* @param ndepth (unsigned int&) number of depth in the grid3D (3rd dimension)
		* @param xllcorner (double&) x-coordinate of lower left corner
		* @param yllcorner (double&) y-coordinate of lower left corner
		* @param latitude (double&) decimal latitude
		* @param longitude (double&) decimal longitude
		* @param cellsize (double&) value for cellsize in grid3D
		* @param grid3D_in (Array\<double\>&) grid to be copied by value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const Array3D<double>& grid3D_in);

		/**
		* @brief check if the current Grid3DObject has the same geolocalization attributes
		* as another Grid3DObject. The grid coordinates (xllcorner & yllcorner) are NOT
		* checked as these might be tweaked for convenience (like between input grid and local grid)
		* @param target (Grid3DObject) grid to compare to
		* @return (bool) true if same geolocalization
		*/
		bool isSameGeolocalization(const Grid3DObject& target);

		Array3D<double> grid3D;
		unsigned int ncols, nrows, ndepth;
		double xllcorner, yllcorner, cellsize;
		double latitude, longitude;
		//std::vector<double> thickness;

 protected:
		void setValues(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude, const double& cellsize);

		void checkCoordinates(const MapProj& proj);
};

#endif
