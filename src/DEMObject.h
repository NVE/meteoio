#ifndef DEMOBJECT_H
#define DEMOBJECT_H

#include "Array2D.h"
#include "Grid2DObject.h"
#include "IOUtils.h"

/**
 * @class DEMObject
 * @brief A class to represent DEMs: reads elevation grids, computes local slope, azimuth, curvature.
 * The nodata parameter is supposed to be IOUtils::nodata.
 * @author Gaël Rosset - Mathias Bavay
 * @date   2009-07-20
 */
class DEMObject: public Grid2DObject {
	public:
		CArray2D<double> slope;
		CArray2D<double> azi;
		CArray2D<double> curvature;
		CArray2D<double> Nx, Ny, Nz;
		double min_altitude, max_altitude;
		double min_slope, max_slope;
		double min_curvature, max_curvature;
		
		DEMObject();
		
		DEMObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
			const double& xllcorner_in, const double& yllcorner_in,
			const double& latitude_in, const double& longitude_in,
			const double& cellsize_in);
	
		DEMObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
			const double& xllcorner_in, const double& yllcorner_in,
			const double& latitude_in, const double& longitude_in,
			const double& cellsize_in, const CArray2D<double>& altitude_in);
		
		DEMObject(const Grid2DObject& dem_in);

		DEMObject (const DEMObject& dem_in, const unsigned int start_col, const unsigned int nb_cols);
		
		void update();

	private:
		void CalculateAziSlopeCurve(const int& Hick);
		double CalculateAzi(const double& Nx, const double& Ny, const double& Nz, const double& slope);
		void CalculateHickNormal(const unsigned int& i, const unsigned int& j, const double& dx_sum, const double& dy_sum);
		void CalculateCorripioNormal(const unsigned int& i, const unsigned int& j, const double& dx_sum, const double& dy_sum);
		double getCurvature(const unsigned int& i, const unsigned int& j);
	
		void CalculateSurfaceDeltas(const unsigned int& i, const unsigned int& j,double *dx1, double *dx2, double *dx3, double *dy1, double *dy2, double *dy3);
		
		double OppositeDir(const double& z1, const double& z2, const double& z);
		void updateAllMinMax();

#ifdef _POPC_
	public:
		void Serialize(POPBuffer &buf, bool pack);
#endif
};

#endif
