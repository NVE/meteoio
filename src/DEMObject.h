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
		Array2D<double> slope;
		Array2D<double> azi;
		Array2D<double> curvature;
		Array2D<double> Nx, Ny, Nz;
		double min_altitude, max_altitude;
		double min_slope, max_slope;
		double min_curvature, max_curvature;

		typedef enum SLOPE_TYPE {
			DFLT, ///< whatever algorithm that has been defined as default
			HICK, ///< maximum downhill slope method (Dunn and Hickey, 1998)
			CORR ///< surface normal vector using the two triangle method (Corripio, 2002) and eight-neighbor algorithm (Horn, 1981) for border cells
		} slope_type;
		
		DEMObject();
		
		DEMObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
			const double& xllcorner_in, const double& yllcorner_in,
			const double& latitude_in, const double& longitude_in,
			const double& cellsize_in);
	
		DEMObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
			const double& xllcorner_in, const double& yllcorner_in,
			const double& latitude_in, const double& longitude_in,
			const double& cellsize_in, const Array2D<double>& altitude_in, 
			const bool& _update=true);
		
		DEMObject(const Grid2DObject& dem_in, const bool& _update=true);

		DEMObject (const DEMObject& _dem,
				const unsigned int& _nx, const unsigned int& _ny, //Point in the plane
				const unsigned int& _ncols, const unsigned int& _nrows, //dimensions of the sub-plane
				const bool& _update=true);
		
		void update(const string& algorithm);
		void update(const slope_type& algorithm=DFLT);
		void updateAllMinMax();

	private:
		void CalculateAziSlopeCurve(const slope_type& algorithm);
		double CalculateAzi(const double& Nx, const double& Ny, const double& Nz, const double& slope);
		void CalculateHickNormal(const unsigned int& i, const unsigned int& j, const double& dx_sum, const double& dy_sum);
		void CalculateCorripioNormal(const unsigned int& i, const unsigned int& j, const double& dx_sum, const double& dy_sum);
		double getCurvature(const unsigned int& i, const unsigned int& j);
	
		void CalculateSurfaceDeltas(const unsigned int& i, const unsigned int& j,double *dx1, double *dx2, double *dx3, double *dy1, double *dy2, double *dy3);
		
		double OppositeDir(const double& z1, const double& z2, const double& z);
		double safeGet(const long& i, const long& j);

		static const slope_type dflt_algorithm;

#ifdef _POPC_
	public:
		virtual void Serialize(POPBuffer &buf, bool pack);
#endif
};

#endif
