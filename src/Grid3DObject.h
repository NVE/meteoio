#ifndef __GRID3DOBJECT_H__
#define __GRID3DOBJECT_H__

#include "Array3D.h"
#include "IOExceptions.h"
#include "IOUtils.h"

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
		* @param ncols (unsigned int&) number of colums in the grid3D
		* @param nrows (unsigned int&) number of rows in the grid3D
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
		* @param xllcorner (double&) x-coordinate of lower left corner
		* @param yllcorner (double&) y-coordinate of lower left corner
		* @param latitude (double&) decimal latitude
		* @param longitude (double&) decimal longitude
		* @param cellsize (double&) value for cellsize in grid3D
		* @param grid3D_in (Array<double>&) grid to be copied by value
		*/
		void set(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude,
			const double& cellsize, const Array3D<double>& grid3D_in);


		Array3D<double> grid3D;
		unsigned int ncols, nrows, ndepth;
		double xllcorner, yllcorner, cellsize;
		double latitude, longitude;
		//std::vector<double> thickness;

 protected:
		void setValues(const unsigned int& ncols, const unsigned int& nrows, const unsigned int& ndepth,
			const double& xllcorner, const double& yllcorner,
			const double& latitude, const double& longitude, const double& cellsize);

		void checkCoordinates();
};

#endif
