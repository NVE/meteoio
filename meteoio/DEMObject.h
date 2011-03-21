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
#ifndef DEMOBJECT_H
#define DEMOBJECT_H

#include <meteoio/Array.h>
#include <meteoio/Array2D.h>
#include <meteoio/Grid2DObject.h>
#include <meteoio/IOUtils.h>

#include <cmath>
#include <limits>

namespace mio {

#ifdef _POPC_
class DEMObjectDummy {}; //HACK for POPC
#endif

/**
 * @class DEMObject
 * @brief A class to represent DEMs: reads elevation grids, computes local slope, azimuth, curvature.
 * The nodata parameter is supposed to be IOUtils::nodata.
 *
 * @ingroup data_str
 * @author Gaël Rosset - Mathias Bavay
 * @date   2009-07-20
 */
class DEMObject : public Grid2DObject {
	public:
		Array2D<double> slope;
		Array2D<double> azi;
		Array2D<double> curvature;
		Array2D<double> Nx, Ny, Nz;
		double min_altitude, max_altitude;
		double min_slope, max_slope;
		double min_curvature, max_curvature;

		///Keywords for slope computation algorithm
		typedef enum SLOPE_TYPE {
			DFLT, ///< whatever algorithm that has been defined as default
			FLEM, ///< four nearest neighbors (Fleming and Hoffer, 1979). It seems to be the same as (Zevenbergen and Thorne, 1987)
			HICK, ///< maximum downhill slope method (Dunn and Hickey, 1998)
			HORN, ///< eight neighbor algorithm (Horn, 1981) as used by ArcGIS. It seems to be the same as (Corripio, 2002) but for border cells.
			CORR, ///< surface normal vector using the two triangle method (Corripio, 2002) and eight-neighbor algorithm (Horn, 1981) for border cells
			D8 ///< discretized azimuth directions (angles for N, NE, etc) and slope rounded to nearest integer
		} slope_type;

		///Keywords for automatic update of parameters. They can be combined with "|"
		typedef enum UPDATE_TYPE {
			NO_UPDATE=0, ///< no updates at all
			SLOPE=1, ///< update the slopes
			NORMAL=2, ///< update the normals
			CURVATURE=4 ///< update the curvatures
		} update_type;
		
		DEMObject(const slope_type& i_algorithm=DFLT);
		
		DEMObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
		          const double& cellsize_in, const Coords& llcorner_in, const slope_type& i_algorithm=DFLT);
	
		DEMObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
		          const double& cellsize_in, const Coords& llcorner_in, const Array2D<double>& altitude_in,
		          const bool& i_update=true, const slope_type& i_algorithm=DFLT);
		
		DEMObject(const Grid2DObject& dem_in, const bool& i_update=true, const slope_type& i_algorithm=DFLT);

		DEMObject (const DEMObject& i_dem,
		           const unsigned int& i_nx, const unsigned int& i_ny, //Point in the plane
		           const unsigned int& i_ncols, const unsigned int& i_nrows, //dimensions of the sub-plane
		           const bool& i_update=true, const slope_type& i_algorithm=DFLT);

		void setDefaultAlgorithm(const slope_type& i_algorithm);
		void setUpdatePpt(const update_type& in_update_flag);
		int getUpdatePpt() const;

		void update(const std::string& algorithm);
		void update(const slope_type& algorithm=DFLT);
		void updateAllMinMax();
		void printFailures();
		void sanitize();
		double horizontalDistance(const double& xcoord1, const double& ycoord1, const double& xcoord2, const double& ycoord2);
		double horizontalDistance(Coords point1, const Coords& point2);
		double terrainDistance(Coords point1, const Coords& point2);
		void getPointsBetween(Coords point1, Coords point2, std::vector<GRID_POINT_2D>& vec_points);
		void getPointsBetween(const Coords& point, const double& bearing, std::vector<GRID_POINT_2D>& vec_points);
		double getHorizon(const Coords& point, const double& bearing);
		void getHorizon(const Coords& point, const double& increment, std::vector<double>& horizon);

	private:
		void CalculateAziSlopeCurve(slope_type algorithm);
		double CalculateAspect(const double& Nx, const double& Ny, const double& Nz, const double& slope, const double no_slope=M_PI);
		void CalculateHick(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz);
		void CalculateFleming(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz);
		void CalculateHorn(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz);
		void CalculateCorripio(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz);
		void (DEMObject::*CalculateSlope)(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz);
		double getCurvature(double A[4][4]);

		double steepestGradient(double A[4][4]);
		double lineGradient(const double& A1, const double& A2, const double& A3);
		double fillMissingGradient(const double& delta1, const double& delta2);
		void surfaceGradient(double& dx_sum, double& dy_sum, double A[4][4]);
		double avgHeight(const double& z1, const double &z2, const double& z3);
		void getNeighbours(const unsigned int i, const unsigned int j, double A[4][4]);
		double safeGet(const int i, const int j);

		int update_flag;
		slope_type dflt_algorithm;
		unsigned int slope_failures; ///<contains the number of points that have an elevation but no slope
		unsigned int curvature_failures; ///<contains the number of points that have an elevation but no curvature

#ifdef _POPC_
	public:
		virtual void Serialize(POPBuffer &buf, bool pack);
#endif
};
} //end namespace

#endif
