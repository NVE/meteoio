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
/**
 * @file libinterpol2D.h
 * This is the two 2D meteo interpolation statistical library.
 */
#ifndef INTERPOL2D_H
#define INTERPOL2D_H

#include <meteoio/StationData.h>
#include <meteoio/MeteoData.h>
#include <meteoio/Grid2DObject.h>
#include <meteoio/DEMObject.h>
#include <meteoio/Date.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <cmath>
#include <vector>
#include <iomanip>
#include <assert.h>

#define MAX_INPUT_STATIONS 255
#define GRAVITY	9.80665		     // (m s-2)
#define GAS_CONSTANT_AIR 287.	     // (J kg-1 K-1)

namespace mio {

/**
 * @class Interpol2D
 * @brief A class to perform 2D spatial interpolations.
 * Each parameter to be interpolated declares which interpolation method to use. 
 * Then the class computes the interpolation for each 2D grid point,
 * combining the inputs provided by the available data sources.
 * @author Mathias Bavay
 * @date   2009-01-20
 */

typedef double (*LapseRateProjectPtr)(const double& value, const double& altitude, 
                                      const double& new_altitude, const std::vector<double>& coeffs);
 
class Interpol2D {
	public:
		///Keywords for selecting the regression algorithm to use
		typedef enum REG_TYPES {
			R_CST, ///< no elevation dependence (ie: constant)
			R_LIN ///< linear elevation dependence
		} reg_types;
		
		static void stdPressureGrid2DFill(const DEMObject& dem, Grid2DObject& grid);
		static void constantGrid2DFill(const double& value, const DEMObject& dem, Grid2DObject& grid);
		static void constantLapseGrid2DFill(const double& value, const double& altitude, 
                                                    const DEMObject& dem, const std::vector<double>& vecCoefficients,
                                                    const LapseRateProjectPtr& funcptr, Grid2DObject& grid);
		static void LapseIDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                                     const DEMObject& dem, const std::vector<double>& vecCoefficients,
                                     const LapseRateProjectPtr& funcptr,
                                     Grid2DObject& grid);
		static void IDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                                const DEMObject& dem, Grid2DObject& grid);
		static void SimpleDEMWindInterpolate(const DEMObject& dem, Grid2DObject& VW, Grid2DObject& DW);

		//projections functions
		static double ConstProject(const double& value, const double& altitude, const double& new_altitude,
		                           const std::vector<double>& coeffs);
		static double LinProject(const double& value, const double& altitude, const double& new_altitude, 
		                         const std::vector<double>& coeffs);

		static int LinRegression(const std::vector<double>& data_in, 
		                         const std::vector<double>& elevations, std::vector<double>& coeffs);

		//these should be in the libphysicslaws, that still has to be created...!
		static double lw_AirPressure(const double& altitude);
		static double RhtoDewPoint(double RH, double TA, const short int& force_water);
		static double DewPointtoRh(double TD, double TA, const short int& force_water);

		//some consts
		const static double dflt_temperature_lapse_rate; ///<default lapse rate for temperature(elevation)

	private:
		//generic functions
		static double InvHorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2);
		static double HorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2);
		static double HorizontalDistance(const DEMObject& dem, const int& i, const int& j, 
		                                 const double& X2, const double& Y2);
		static double getReferenceAltitude(const DEMObject& dem);
		
		//core methods
		static void LinRegressionCore(const std::vector<double>& X, const std::vector<double>& Y,
		                              const unsigned int imax, double& a, double& b, double& r);
		static double IDWCore(const double& x, const double& y,
		                      const std::vector<double>& vecData_in,
		                      const std::vector<StationData>& vecStations_in);

		//weighting methods
		double weightInvDist(const double& d2);
		double weightInvDistSqrt(const double& d2);
		double weightInvDist2(const double& d2);
		double weightInvDistN(const double& d2);
		double dist_pow; //power for the weighting method weightInvDistN

	private:
		//static members
		const static double wind_ys; ///coefficient for wind dependency on slope
		const static double wind_yc; ///coefficient for wind dependency on curvature
};
} //end namespace

#endif
