//This is the two 2D meteo interpolation library.
#ifndef INTERPOL2D_H
#define INTERPOL2D_H

#include <vector>
#include "Array2D.h"
#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"

#define MAX_INPUT_STATIONS 255

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
		//available methods for single source and multiple sources interpolations
		typedef enum INTERP_TYPES {I_PRESS, I_RH, I_CST, I_IDWK, I_LAPSE_CST, I_LAPSE_IDWK} interp_types;
		typedef enum REG_TYPES {R_CST, R_LIN} reg_types;
		
		/**
		* @brief Constructor. Sets private members for latter use
		* @param Isingle enum for the type of single data source interpolation strategy
		* @param Imultiple enum for the type of multiple data source interpolation strategy
		* @param sourcesData vector of source data
		* @param sourcesMeta vector of metadata about the sources
		* @param dem Digital Elevation Model object
		*/
		Interpol2D(interp_types Isingle, 
			interp_types Imultiple, 
			const vector<double>& sourcesData, 
			const vector<StationData>& sourcesMeta, 
	     const Grid2DObject& dem);
		
		/**
		* @brief Computes the interpolation using the parameters set by the constructor
		* @param param_out 2D grid containing the interpolated values
		*/
		void calculate(CArray2D<double>& param_out);
		void calculate(CArray2D<double>& param_out, const vector<double>& vecExtraStations, CArray2D<double>& extra_param_in);

	private:
		//generic functions
		double HorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2);
		double HorizontalDistance(const int& i, const int& j, const double& X2, const double& Y2);
		double AvgSources(const vector<double>& data_in);
		void BuildStationsElevations(const vector<StationData>& vecStations_in, vector<double>& vecElevations);
		
		//regressions
		void LinRegressionCore(const vector<double>& X, const vector<double>& Y, double& a, double& b, double& r, const int ignore_index);
		int LinRegression(const vector<double>& data_in, const vector<double>& elevations, vector<double>& coeffs);
		
		//projections functions
		double ConstProject(const double& val, const double& alt, const double& new_alt, const vector<double>& coeffs);
		double LinProject(const double& value, const double& altitude, const double& new_altitude, const vector<double>& coeffs);
		
		///Member function pointer
		double (Interpol2D::*LapseRateProject)(const double& value, 
							const double& altitude, 
							const double& new_altitude, 
							const vector<double>& coeffs); 
		
		//filling functions
		void StdPressureFill(CArray2D<double>& param, const CArray2D<double>& topoheight);
		void ConstFill(CArray2D<double>& param, const double& value);
		void LapseConstFill(CArray2D<double>& param, const double& value, const double& altitude, const CArray2D<double>& topoHeight);
		
		
		void LapseIDWKrieging(CArray2D<double>& T, const CArray2D<double>& topoHeight, 
					const vector<double>& vecData_in, const vector<StationData>& vecStations_in);
		double IDWKriegingCore(const double& x, const double& y, 
					const vector<double>& vecData_in, const vector<StationData>& vecStations);
		void IDWKrieging(CArray2D<double>& T, const vector<double>& data_in, const vector<StationData>& vecStations);
		

	private:
		//static members
		const static double ref_altitude;
		const static double dflt_temperature_lapse_rate;
		
		interp_types single_type, multiple_type;	//interpolations choices for single and multiple sources
		//reg_types reg_type;				//choice of regression methods
		
		double xllcorner, yllcorner, cellsize;	//for more transparent access, to set when instanciating
		unsigned int nx, ny;				//2D grid dimensions
		unsigned int InputSize;
		
		vector<double> vecCoefficients;	      ///<Regression coefficients 0-3
		const Grid2DObject& dem;                    ///<Reference to be initialized in the constructor
		const CArray2D<double>& InputTopo;          ///<Reference to be initialized in the constructor
		const vector<StationData>& InputMeta;       ///<Reference to be initialized in the constructor
		const vector<double>& inputData;            ///<Reference to be initialized in the constructor
};

#endif
