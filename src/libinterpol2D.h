//This is the two 2D meteo interpolation library.
#ifndef INTERPOL2D_H
#define INTERPOL2D_H

#include <vector>
#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"
#include "DEMObject.h"

#define MAX_INPUT_STATIONS 255
#define GRAVITY	9.80665		     // (m s-2)
#define GAS_CONSTANT_AIR 287.	     // (J kg-1 K-1)

/**
 * @class Interpol2D
 * @brief A class to perform 2D spatial interpolations.
 * Each parameter to be interpolated declares which interpolation method to use
 * for single and multiple data sources. Then the class computes the interpolation for each 2D grid point,
 * combining the inputs provided by the available data sources.
 * @author Mathias Bavay
 * @date   2009-01-20
 */
 
class Interpol2D {
	public:
		/**
		* @enum INTERP_TYPES
		* keywords for selecting the spatial interpolation algorithm
		*/
		//available methods for single source and multiple sources interpolations
		typedef enum INTERP_TYPES {
			I_PRESS, ///< standard air pressure interpolation
			I_RH, ///< relative humidity interpolation
			I_VW, ///< wind velocity interpolation (using a heuristic terrain effect)
			I_CST, ///< constant fill
			I_IDWK, ///< Inverse Distance Weighting fill
			I_LAPSE_CST, ///< constant fill with an elevation lapse rate
			I_LAPSE_IDWK ///< Inverse Distance Weighting with an elevation lapse rate fill
		} interp_types;
		typedef enum REG_TYPES {
			R_CST, ///< no elevation dependence (ie: constant)
			R_LIN ///< linear elevation dependence
		} reg_types;
		
		Interpol2D(interp_types Isingle, 
			interp_types Imultiple, 
			const std::vector<double>& sourcesData, 
			const std::vector<StationData>& sourcesMeta, 
	    	const DEMObject& dem);
		
		void calculate(Grid2DObject& param_out);
		void calculate(Grid2DObject& param_out, const std::vector<double>& vecExtraStations, Grid2DObject& extra_param_in);

	private:
		//generic functions
		double HorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2);
		double HorizontalDistance(const int& i, const int& j, const double& X2, const double& Y2);
		double AvgSources(const std::vector<double>& data_in);
		void BuildStationsElevations(const std::vector<StationData>& vecStations_in, std::vector<double>& vecElevations);
		
		//regressions
		void LinRegressionCore(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, const int ignore_index);
		int LinRegression(const std::vector<double>& data_in, const std::vector<double>& elevations, std::vector<double>& coeffs);
		
		//projections functions
		double ConstProject(const double& val, const double& alt, const double& new_alt, const std::vector<double>& coeffs);
		double LinProject(const double& value, const double& altitude, const double& new_altitude, const std::vector<double>& coeffs);
		
		///Member function pointer
		double (Interpol2D::*LapseRateProject)(const double& value, 
							const double& altitude, 
							const double& new_altitude, 
							const std::vector<double>& coeffs); 
		
		//filling functions
		void StdPressureFill(Grid2DObject& param, const DEMObject& topoHeight);
		void ConstFill(Grid2DObject& param, const double& value);
		void LapseConstFill(Grid2DObject& param_out, const DEMObject& topoHeight, const double& value, const double& altitude);
		
		
		void LapseIDWKrieging(Grid2DObject& T, const DEMObject& topoHeight,
				const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in);
		double IDWKriegingCore(const double& x, const double& y, 
						   const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations);
		void IDWKrieging(Grid2DObject& T, const std::vector<double>& data_in, const std::vector<StationData>& vecStations);
		void SimpleDEMWindInterpolate(Grid2DObject& VW, Grid2DObject& DW);
		double RhtoDewPoint(double RH, double TA, const short int force_water);
		double DewPointtoRh(double TD, double TA, const short int force_water);
		double lw_AirPressure(const double altitude);

	private:
		//static members
		const static double dflt_temperature_lapse_rate;///default lapse rate for temperature(elevation)
		const static double wind_ys;			///coefficient for wind dependency on slope
		const static double wind_yc;			///coefficient for wind dependency on curvature

		double ref_altitude;				///elevation to use for common elevation reprojection
		
		interp_types single_type, multiple_type;	//interpolations choices for single and multiple sources
		//reg_types reg_type;				//choice of regression methods
		
		unsigned int InputSize;
		
		std::vector<double> vecCoefficients;		///<Regression coefficients 0-3
		const DEMObject& dem;				///<Reference to be initialized in the constructor
		const std::vector<StationData>& InputMeta;	///<Reference to be initialized in the constructor
		const std::vector<double>& inputData;		///<Reference to be initialized in the constructor
};

#endif
