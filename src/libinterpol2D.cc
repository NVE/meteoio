//This is the two 2D meteo interpolation library.
#include <cmath>
#include <vector>
#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"
#include "Date_IO.h"
#include "IOExceptions.h"
#include "Laws.h"
#include "libinterpol2D.h"
#include "IOUtils.h"
#include "DEMObject.h"


const double Interpol2D::dflt_temperature_lapse_rate = -0.0065;
const double Interpol2D::wind_ys = 0.58;
const double Interpol2D::wind_yc = 0.42;

/**
* @brief Constructor. Sets private members for latter use
* @param Isingle enum for the type of single data source interpolation strategy
* @param Imultiple enum for the type of multiple data source interpolation strategy
* @param sourcesData vector of source data
* @param sourcesMeta vector of metadata about the sources
* @param dem Digital Elevation Model object
*/
Interpol2D::Interpol2D(interp_types Isingle, 
		       interp_types Imultiple, 
		       const vector<double>& vecData, 
		       const vector<StationData>& vecMeta, 
		       const DEMObject& dem_in) : dem(dem_in), InputTopo(dem.grid2D), InputMeta(vecMeta), inputData(vecData)
{

	single_type   = Isingle;
	multiple_type = Imultiple;
	
	xllcorner = dem.xllcorner;	//TODO: instead, access by address of DEM member
	yllcorner = dem.yllcorner;
	cellsize  = dem.cellsize;
	nx        = dem.ncols;
	ny        = dem.nrows;
	if(dem.min_altitude!=IOUtils::nodata && dem.max_altitude!=IOUtils::nodata) {
		//we use the median elevation as the reference elevation for reprojections
		ref_altitude = 0.5 * (dem.min_altitude+dem.max_altitude); 
	} else {
		//since there is nothing else that we can do, we use an arbitrary elevation
		ref_altitude = 1500.;
	}

	vecCoefficients.resize(4, 0.0); //Coefficients 0-3 TODO!! make it more dynamic!
	InputSize = (unsigned int)inputData.size();

}

//Usefull functions
/**
* @brief Computes the horizontal distance between points, given by coordinates in a geographic grid
* @param X1 (const double) first point's X coordinate
* @param Y1 (const double) first point's Y coordinate
* @param X2 (const double) second point's X coordinate
* @param Y2 (const double) second point's Y coordinate
* @return (double) distance in m
*/
double Interpol2D::HorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2)
{
	//This function computes the horizontaldistance between two points
	//coordinates are given in a square, metric grid system
	return sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) );
}

/**
* @brief Computes the horizontal distance between points, given by their cells indexes
* @param X1 (const double) first point's i index
* @param Y1 (const double) first point's j index
* @param X2 (const double) second point's X coordinate
* @param Y2 (const double) second point's Y coordinate
* @return (double) distance in m
*/
double Interpol2D::HorizontalDistance(const int& i, const int& j, const double& X2, const double& Y2)
{
	//This function computes the horizontal distance between two points
	//coordinates are given in a square, metric grid system
	//for grid points toward real coordinates
	const double X1 = (xllcorner+i*cellsize);
	const double Y1 = (yllcorner+j*cellsize);
	
	return sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) );
}

/**
* @brief Returns the average of data contained in a vector
* @param data_in (vector\<double\>) vector of data
* @return (double) average of the vector
*/
double Interpol2D::AvgSources(const vector<double>& data_in)
{
	//This function computes the average of all the data sources
	double avg=0;

	for (unsigned int i=0; i<(unsigned int)data_in.size(); i++) {
		avg += data_in[i];
	}
	return (avg/(double)data_in.size());
}

/**
* @brief From a vector of StationData, builds a vector of elevations
* @param vecStations_in (vector\<StationData\>) vector of StationData objects
* @param vecElevations (vector\<double\>) vector of elevations (same order as the input)
*/
void Interpol2D::BuildStationsElevations(const vector<StationData>& vecStations_in, vector<double>& vecElevations)
{
	//from a vector of stations meta data, builds a 1D array
	for (unsigned int i=0; i<(unsigned int)vecStations_in.size(); i++) {
		vecElevations.push_back(vecStations_in[i].getAltitude());
	}
}

//Data regression models
/**
* @brief Computes the linear regression coefficients fitting the points given as X and Y in two vectors
* the linear regression has the form Y = aX + b with a regression coefficient r
* @param X (vector\<double\>) vector of X coordinates
* @param Y (vector\<double\>) vector of Y coordinates (same order as X)
* @param a (double) slope of the linear regression
* @param b (double) origin of the linear regression
* @param r (double) linear regression coefficient
* @param ignore_index (const int) if positive, index of a point to exclude from the regression
*/
void Interpol2D::LinRegressionCore(const vector<double>& X, const vector<double>& Y, double& a, double& b, double& r, const int ignore_index)
{
	//finds the linear regression for points (x,y,z,Value)
	//if ignore_index>=0, ignore given index (as a way to remopve a point from the interpolation)
	double x_avg=0, y_avg=0;
	double sx=0, sy=0, sxy=0;
	const int imax = (int)X.size();

	//computing x_avg and y_avg
	for (int i=0; i<imax; i++) {
		if (i!=ignore_index) {
			x_avg += X[i];
			y_avg += Y[i];
		}
	}
	if (ignore_index>=0) {
		x_avg /= (double)(imax - 1);
		y_avg /= (double)(imax - 1);
	} else {
		x_avg /= (double)imax;
		y_avg /= (double)imax;
	}

	//computing sx, sy, sxy
	for (int i=0; i<imax; i++) {
		if (i!=ignore_index) {
			sx += (X[i]-x_avg) * (X[i]-x_avg);
			sy += (Y[i]-y_avg) * (Y[i]-y_avg);
			sxy += (X[i]-x_avg) * (Y[i]-y_avg);
		}
	}

	//computing the regression line
	a = sxy / sx;			//a
	b = y_avg - a*x_avg;		//b
	r = sxy / sqrt(sx*sy);		//r
}

/**
* @brief Computes the linear regression coefficients fitting the points given as X and Y in two vectors
* the linear regression has the form Y = aX + b with a regression coefficient r. If the regression coefficient is not good enough, a bad point is looked removed.
* @param X (vector\<double\>) vector of X coordinates
* @param Y (vector\<double\>) vector of Y coordinates (same order as X)
* @param coeffs (vector\<double\>) a,b,r coefficients in a vector
* @return (int) EXIT_SUCCESS or EXIT_FAILURE
*/
int Interpol2D::LinRegression(const vector<double>& X, const vector<double>& Y, vector<double>& coeffs)
{
	//finds the linear regression for points (x,y,z,Value)
	const double r_thres=0.7;
	double a,b,r;
	
	if ((unsigned int)X.size()==2) {
		printf("[W] only two points for linear regression!\n"); //HACK
	}
	if((unsigned int)X.size()<2) { //this should not be needed, we should have refrained from calling LinRegression in such a case
		printf("Not enough data point for linear regression!\n"); //HACK
		coeffs[1]=0.;
		coeffs[2]=X[1];
		coeffs[3]=1.;
		return EXIT_FAILURE;
	}

	LinRegressionCore(X, Y, coeffs[1], coeffs[2], coeffs[3],-1);
	if(coeffs[3]<r_thres && (unsigned int)X.size()>=4) { //is r good enough? (and we need at least 3 points for a linear regression)
		//if r is not good enough, we try by removing 1 point and keep the best result from these N-1 calculations...
		for (unsigned int i=0; i<(unsigned int)X.size(); i++) {
			LinRegressionCore(X, Y, a, b, r,(int)i);
			if (fabs(r)>fabs(coeffs[3])) {
				coeffs[1]=a;
				coeffs[2]=b;
				coeffs[3]=r;
			}
		}
	}
	//check if r is reasonnable
	if (fabs(coeffs[3])<r_thres) {
		printf("Poor regression coefficient: %g\n",coeffs[3]); //HACK
		//for (unsigned int i=0; i<X.size(); i++){
		// printf("%g %g\n",X[i],Y[i]);
		// }
		//throw IOException("Poor linear regression coefficient", AT);
	}

	return EXIT_SUCCESS;
}


//Now, the core interpolation functions: they project a given parameter to a reference altitude, given a constant lapse rate
//example: Ta projected to 1500m with a rate of -0.0065K/m
/**
* @brief Projects a given parameter to another elevation: 
* This implementation keeps the value constant as a function of the elevation
* @param value (const double) original value
* @param altitude (const double) altitude of the original value (unused)
* @param new_altitude (const double) altitude of the reprojected value (unused)
* @param coeffs (const vector\<double\>) coefficients to use for the projection (unused)
* @return (double) reprojected value
*/
double Interpol2D::ConstProject(const double& value, const double& altitude, const double& new_altitude, const vector<double>& coeffs)
{
	(void) altitude;	//to avoid the compiler seeing them as unused
	(void) new_altitude;
	(void) coeffs;
	return value;
}

/**
* @brief Projects a given parameter to another elevation: 
* This implementation assumes a linear dependency of the value as a function of the elevation
* @param value (const double) original value
* @param altitude (const double) altitude of the original value
* @param new_altitude (const double) altitude of the reprojected value
* @param coeffs (const vector\<double\>) coefficients to use for the projection
* @return (double) reprojected value
*/
double Interpol2D::LinProject(const double& value, const double& altitude, const double& new_altitude, const vector<double>& coeffs)
{
	//linear lapse: coeffs must have been already computed
	if (coeffs.size()<1) {
		throw IOException("Linear regression coefficients not initialized", AT);
	}
	return (value + coeffs[1] * (new_altitude - altitude));
}

//Filling Functions
/**
* @brief Grid filling function: 
* This implementation builds a standard air pressure as a function of the elevation
* @param param (Array2D\<double\>) 2D array to fill
* @param topoheight (Array2D\<double\>) array of elevations (dem)
*/
void Interpol2D::StdPressureFill(Array2D<double>& param, const Array2D<double>& topoheight) {
	//provide each point with an altitude dependant pressure... it is worth what it is...
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			if (topoheight(i,j)!=IOUtils::nodata) {
				param(i,j) = lw_AirPressure(topoheight(i,j));
			} else {
				param(i,j) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief Grid filling function: 
* This implementation fills the grid with a constant value
* @param param (Array2D\<double\>) 2D array to fill
* @param value (double) value to put in the grid
*/
void Interpol2D::ConstFill(Array2D<double>& param, const double& value)
{
	//fills a data table with constant values
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			param(i,j) = value;	//should we here write nodata when it is nodata in dem?
		}
	}
}

/**
* @brief Grid filling function: 
* This implementation fills a flat grid with a constant value, and then reproject it to the terrain's elevation.
* for example, the air temperature measured at one point at 1500m would be given as value, the 1500m as altitude and the dem would allow to reproject this temperature on the full DEM using a default lapse rate.
* @param param_out (Array2D\<double\>) 2D array to fill
* @param value (double) value to put in the grid
* @param altitude (double) altitude of the "value"
* @param topoheight (Array2D\<double\>) array of elevations (dem)
*/
void Interpol2D::LapseConstFill(Array2D<double>& param_out, const double& value, const double& altitude, const Array2D<double>& topoHeight)
{
	//fills a data table with constant values and then reprojects it to the DEM's elevation from a given altitude
	//the laspe rate parameters must have been set before
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			if (topoHeight(i,j)!=IOUtils::nodata) {
				param_out(i,j) = (this->*LapseRateProject)(value, altitude,topoHeight(i,j), vecCoefficients);
			} else {
				param_out(i,j) = IOUtils::nodata;
			}
		}
	}
}

double Interpol2D::IDWKriegingCore(const double& x, const double& y, 
				   const vector<double>& vecData_in, const vector<StationData>& vecStations)
{
	//The value at any given cell is the sum of the weighted contribution from each source
	double parameter=0., norm=0., weight;
	
	for (unsigned int i=0; i<(unsigned int)vecStations.size(); i++) {
		weight=1./(HorizontalDistance(x, y, vecStations[i].getEasting(), vecStations[i].getNorthing()) + 1e-6);
		parameter += weight*vecData_in[i];
		norm += weight;
	}
	return (parameter/norm);	//normalization
}


void Interpol2D::LapseIDWKrieging(Array2D<double>& T, const Array2D<double>& topoHeight, 
				  const vector<double>& vecData_in, const vector<StationData>& vecStations_in)
{
	//multiple source stations: lapse rate projection, IDW Krieging, re-projection
	vector<double> vecTref(vecStations_in.size(), 0.0); // init to 0.0
	
	for (unsigned int i=0; i<(unsigned int)vecStations_in.size(); i++) {
		vecTref[i] = (this->*LapseRateProject)(vecData_in[i], vecStations_in[i].getAltitude(), ref_altitude, vecCoefficients);
	}
	
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			if (topoHeight(i,j)!=IOUtils::nodata) {
				T(i,j) = IDWKriegingCore((xllcorner+i*cellsize), (yllcorner+j*cellsize),vecTref, vecStations_in);
				T(i,j) = (this->*LapseRateProject)(T(i,j), ref_altitude, topoHeight(i,j), vecCoefficients);
			} else {
				T(i,j) = IOUtils::nodata;
			}
		}
	}
}

void Interpol2D::IDWKrieging(Array2D<double>& T, const vector<double>& vecData_in, const vector<StationData>& vecStations)
{
	//multiple source stations: simple IDW Krieging
	for (unsigned int i=0; i<nx; i++) {
		for(unsigned int j=0; j<ny; j++) {	//should we write nodata when dem=nodata?
			T(i,j) = IDWKriegingCore((xllcorner+i*cellsize), (yllcorner+j*cellsize), vecData_in, vecStations);
		}
	}
}

/**
* @brief Computes the interpolation using the parameters set by the constructor
* @param param_out 2D grid containing the interpolated values
*/
void Interpol2D::calculate(Array2D<double>& param_out)
{
	unsigned short int flag_ok=0;
	vector<double> vecStationElevations;
	
	if (InputSize==0) {	//no data
		if (single_type == I_PRESS) {
			StdPressureFill(param_out, InputTopo);
			flag_ok=1;
		}
	
		if (flag_ok==0) {
                  	throw IOException("Wrong interpolation type for zero data source", AT);
		}
	}
	
	if (InputSize==1) { //single data source
		if (single_type==I_CST) {
			ConstFill(param_out, inputData[0]);
			flag_ok=1;
	
		} else if (single_type==I_LAPSE_CST) {
			vecCoefficients[1] = dflt_temperature_lapse_rate;//HACK, it should depend on a user given value!
			LapseRateProject = &Interpol2D::LinProject;
			LapseConstFill(param_out, inputData[0], InputMeta[0].getAltitude(), InputTopo);
			flag_ok=1;
		}
	
		if (flag_ok==0) {
			throw IOException("Wrong interpolation type for single data source", AT);
		}
	}
	
	if (InputSize>1) { //multiple data sources
		if (multiple_type==I_CST) {
			ConstFill(param_out, AvgSources(inputData));
			flag_ok=1;
		} else if (multiple_type==I_IDWK) {
			IDWKrieging(param_out, inputData, InputMeta);
			flag_ok=1;
		} else if (multiple_type==I_LAPSE_CST) {
			BuildStationsElevations(InputMeta, vecStationElevations);
			vecCoefficients[1] = dflt_temperature_lapse_rate;//HACK, it should depend on a user given value!
			LapseRateProject = &Interpol2D::ConstProject;
			LapseConstFill(param_out, AvgSources(inputData), AvgSources(vecStationElevations), InputTopo);
			flag_ok=1;
		} else if (multiple_type==I_LAPSE_IDWK) {
			BuildStationsElevations(InputMeta, vecStationElevations);
			LinRegression(vecStationElevations, inputData, vecCoefficients);
			LapseRateProject = &Interpol2D::LinProject;
			LapseIDWKrieging(param_out, InputTopo, inputData, InputMeta);
			flag_ok=1;
		}

		if (flag_ok==0) {
			throw IOException("Wrong interpolation type for multiple data sources", AT);
		}
	}
}

void Interpol2D::calculate(Array2D<double>& param_out, const vector<double>& vecExtraData, Array2D<double>& extra_param_in) {
	unsigned short int flag_ok=0;
	vector<double> vecStationElevations;
	
	if (InputSize==0) { //no data
		if (flag_ok==0) {
			throw IOException("Wrong interpolation type for zero data source", AT);
		}
	}
	
	if (InputSize==1) {
		//single data source
		if (flag_ok==0) {
			throw IOException("Wrong interpolation type for single data source", AT);
		}
	}

	if (InputSize>1) {
		//multiple data sources
		if (multiple_type==I_RH) {
			//here, RH->Td, interpolations, Td->RH
			vector<double> vecTdStations(inputData.size(), 0.0); // init to 0.0

			//Compute dew point temperatures at stations
			for (unsigned int i=0; i<(unsigned int)inputData.size(); i++) {
				vecTdStations[i] = RhtoDewPoint(inputData[i],vecExtraData[i], 1);
			}
			
			//Krieging on Td
			BuildStationsElevations(InputMeta, vecStationElevations);
			LinRegression(vecStationElevations, vecTdStations, vecCoefficients);
			LapseRateProject = &Interpol2D::LinProject;
			LapseIDWKrieging(param_out, InputTopo, vecTdStations, InputMeta);

			//Recompute Rh from the interpolated td
			for (unsigned int i=0;i<nx;i++) {
				for (unsigned int j=0;j<ny;j++) {
					param_out(i,j) = DewPointtoRh(param_out(i,j),extra_param_in(i,j), 1);
				}
			}
			
			flag_ok=1;
		} else if (multiple_type==I_VW) {
			//Krieging
			BuildStationsElevations(InputMeta, vecStationElevations);
			LinRegression(vecStationElevations, inputData, vecCoefficients);
			LapseRateProject = &Interpol2D::LinProject;
			LapseIDWKrieging(param_out, InputTopo, inputData, InputMeta);
			SimpleDEMWindInterpolate(param_out, extra_param_in);//HACK: extra_param_in is DW, how do we bring it out?
			
			flag_ok=1;
		}

		if (flag_ok==0) {
			throw IOException("Wrong interpolation type for multiple data sources", AT);
		}
	}
}

void Interpol2D::SimpleDEMWindInterpolate(Array2D<double>& VW, Array2D<double>& DW)
{
//This method computes the speed of the wind and returns a table in 2D with this values
//This Wind interpolation is similar to Liston and Elder (2006)
	double speed;			// Wind speed (m s-1)
	double dir;			// Wind direction
	double u;			// Zonal component u (m s-1)
	double v;			// Meridional component v (m s-1)
	double beta;			// Terrain slope
	double azi;			// Topographic slope azimuth
	double curvature;		// Topographic curvature
	double slopeDir;		// Slope in the direction of the wind
	double Ww;			// Wind weighting
	double Od;			// Diverting factor
	
	// For each cell
	for (unsigned int i=0;i<nx-1;i++) {
		for (unsigned int j=0;j<ny-1;j++){
			// Get data
			speed = VW(i,j);
			dir = DW(i,j) * ((M_PI) / 180.);		

			//Speed and direction converted to zonal et meridional
			//components 
			u = (-1.) * (speed * sin(dir));
			v = (-1.) * (speed * cos(dir));

			// Converted back to speed and direction
			speed = sqrt(u*u + v*v);
			dir = (1.5 * M_PI) - atan(v/u);
			
			// Get the slope of this cell
			beta = dem.slope(i, j);
			
			// Get the slope azimuth of this cell
			azi = dem.azi(i, j);
			
			// Get the curvature of this cell
			curvature = dem.curvature(i, j);

			//normalize curvature and beta. 
			//Note: it should be slopeDir instead of beta, but beta is more efficient
			//to compute (only once for each dem) and it should not be that different...
			beta = (beta - dem.min_slope)/(dem.max_slope - dem.min_slope) - 0.5;
			curvature = (curvature - dem.min_curvature)/(dem.max_curvature - dem.min_curvature) - 0.5;

			// Calculate the slope in the direction of the wind
			slopeDir = beta * cos(dir - azi);
	
			// Calculate the wind weighting factor
			Ww = 1. + wind_ys * slopeDir + wind_yc * curvature;

			// Calculate the terrain-modified wind speed
			VW(i, j) = Ww * speed;

			// Modify the wind direction by a diverting factor
			Od = -0.5 * slopeDir * sin(2.*(azi - dir));

			// Add this factor to the wind direction and
			// transform this in degrees
			DW(i, j) = (dir + Od) * (180. / (M_PI));
			if( DW(i, j)>360. ) {
				DW(i, j) -= 360.;
			}
		}
	}
}
