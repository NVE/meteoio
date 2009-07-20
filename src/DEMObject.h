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
		
		DEMObject():Grid2DObject(){};
		DEMObject(Grid2DObject& dem_in);
		
		void Update();

		virtual DEMObject* sub(const unsigned int start_col, const unsigned int nb_cols);

	private:
		void CalculateAziSlopeCurve(const int Hick);
		double CalculateAzi(const double Nx, const double Ny, const double Nz, const double slope);
		void CalculateHickNormal(const unsigned int i, const unsigned int j, const double dx_sum, const double dy_sum);
		void CalculateCorripioNormal(const unsigned int i, const unsigned int j, const double dx_sum, const double dy_sum);
		double getCurvature(const unsigned int i, const unsigned int j);
	
		void CalculateSurfaceDeltas(const unsigned int i, const unsigned int j,double *dx1, double *dx2, double *dx3, double *dy1, double *dy2, double *dy3);
		
		double OppositeDir(const double z1, const double z2, const double z);

#ifdef _POPC_
	public:
		void Serialize(POPBuffer &buf, bool pack);
#endif
};

#endif
