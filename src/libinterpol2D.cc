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


const double Interpol2D::dflt_temperature_lapse_rate = -0.0065;
const double Interpol2D::ref_altitude = 1500.;

Interpol2D::Interpol2D(interp_types Isingle, 
		       interp_types Imultiple, 
		       const vector<double>& vecData, 
		       const vector<StationData>& vecMeta, 
		       const Grid2DObject& dem_in) : dem(dem_in), InputTopo(dem.grid2D), InputMeta(vecMeta), inputData(vecData){

	single_type   = Isingle;
	multiple_type = Imultiple;
	
	xllcorner = dem.xllcorner;	//TODO: instead, access DEM members in functions
	yllcorner = dem.yllcorner;
	cellsize  = dem.cellsize;
	nx        = dem.ncols;
	ny        = dem.nrows;

	vecCoefficients.resize(4, 0.0); //Coefficients 0-3 TODO!! make it more dynamic!
	InputSize = (unsigned int)inputData.size();

	//Debugging output:
	//std::cerr << AT << "\n    inputData.size()==" << inputData.size() 
	//	  << "\n    dim: " << nx << " x " << ny << std::endl;
}

//Usefull functions
double Interpol2D::HorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2)
{
	//This function computes the horizontaldistance between two points
	//coordinates are given in the Swiss grid system
	return sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) );
}

double Interpol2D::HorizontalDistance(const int& i, const int& j, const double& X2, const double& Y2)
{
	//This function computes the horizontaldistance between two points
	//for grid points toward real coordinates
	const double X1 = (xllcorner+i*cellsize);
	const double Y1 = (yllcorner+j*cellsize);
	
	return sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) );
}

double Interpol2D::AvgSources(const vector<double>& data_in)
{
	//This function computes the average of all the data sources
	double avg=0;

	for (unsigned int i=0; i<(unsigned int)data_in.size(); i++) {
		avg += data_in[i];
	}
	return (avg/(double)data_in.size());
}

void Interpol2D::BuildStationsElevations(const vector<StationData>& vecStations_in, vector<double>& vecElevations)
{
	//from a vector of stations meta data, builds a 1D array
	for (unsigned int i=0; i<(unsigned int)vecStations_in.size(); i++) {
		vecElevations.push_back(vecStations_in[i].getAltitude());
	}
}

//Data regression models
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
		//THROW IOException("Poor linear regression coefficient", AT);
	}

	return EXIT_SUCCESS;
}


//Now, the core interpolation functions: they project a given parameter to a reference altitude, given a constant lapse rate
//example: Ta projected to 1500m with a rate of -0.0065K/m
double Interpol2D::ConstProject(const double& value, const double& altitude, const double& new_altitude, const vector<double>& coeffs)
{
	(void) altitude;	//to avoid the compiler seeing them as unused
	(void) new_altitude;
	(void) coeffs;
	return value;
}

double Interpol2D::LinProject(const double& value, const double& altitude, const double& new_altitude, const vector<double>& coeffs)
{
	//linear lapse: coeffs must have been already computed
	if (coeffs.size()<1) {
		THROW IOException("Linear regression coefficients not initialized", AT);
	}
	return (value + coeffs[1] * (new_altitude - altitude));
}

//Filling Functions
void Interpol2D::StdPressureFill(CArray2D<double>& param, const CArray2D<double>& topoheight) {
	//provide each point with an altitude dependant pressure... it is worth what it is...
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			if (topoheight(i,j)!=dem.nodata) {
				param(i,j) = lw_AirPressure(topoheight(i,j));
			} else {
				param(i,j) = dem.nodata;
			}
		}
	}
}

void Interpol2D::ConstFill(CArray2D<double>& param, const double& value)
{
	//fills a data table with constant values
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			param(i,j) = value;	//should we here write nodata when it is nodata in dem?
		}
	}
}

void Interpol2D::LapseConstFill(CArray2D<double>& param_out, const double& value, const double& altitude, const CArray2D<double>& topoHeight)
{
	//fills a data table with constant values and then reprojects it to the DEM's elevation from a given altitude
	//the laspe rate parameters must have been set before
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			if (topoHeight(i,j)!=dem.nodata) {
				param_out(i,j) = (this->*LapseRateProject)(value, altitude,topoHeight(i,j), vecCoefficients);
			} else {
				param_out(i,j) = dem.nodata;
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


void Interpol2D::LapseIDWKrieging(CArray2D<double>& T, const CArray2D<double>& topoHeight, 
				  const vector<double>& vecData_in, const vector<StationData>& vecStations_in)
{
	//multiple source stations: lapse rate projection, IDW Krieging, re-projection
	vector<double> vecTref(vecStations_in.size(), 0.0); // init to 0.0
	
	for (unsigned int i=0; i<(unsigned int)vecStations_in.size(); i++) {
		vecTref[i] = (this->*LapseRateProject)(vecData_in[i], vecStations_in[i].getAltitude(), ref_altitude, vecCoefficients);
	}
	
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			if (topoHeight(i,j)!=dem.nodata) {
				T(i,j) = IDWKriegingCore((xllcorner+i*cellsize), (yllcorner+j*cellsize),vecTref, vecStations_in);
				T(i,j) = (this->*LapseRateProject)(T(i,j), ref_altitude, topoHeight(i,j), vecCoefficients);
			} else {
				T(i,j) = dem.nodata;
			}
		}
	}
}


void Interpol2D::IDWKrieging(CArray2D<double>& T, const vector<double>& vecData_in, const vector<StationData>& vecStations)
{
	//multiple source stations: simple IDW Krieging
	for (unsigned int i=0; i<nx; i++) {
		for(unsigned int j=0; j<ny; j++) {	//should we write nodata when dem=nodata?
			T(i,j) = IDWKriegingCore((xllcorner+i*cellsize), (yllcorner+j*cellsize), vecData_in, vecStations);
		}
	}
}


void Interpol2D::calculate(CArray2D<double>& param_out)
{
	unsigned short int flag_ok=0;
	vector<double> vecStationElevations;
	
	if (InputSize==0) {	//no data
		if (single_type == I_PRESS) {
			StdPressureFill(param_out, InputTopo);
			flag_ok=1;
		}
	
		if (flag_ok==0) {
                  	THROW IOException("Wrong interpolation type for zero data source", AT);
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
			THROW IOException("Wrong interpolation type for single data source", AT);
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
		} else if (multiple_type==I_VW) { 
			//Welcome Gaël Rosset!! This is your part!
			//here, do the wind interpolation
			//since this is a very specific method, it makes sense to simply call a function from here
			//and have all the processing managed by it

			//SimpleDEMWindInterpolate(param_out, InputTopo, inputData, InputMeta);
			
			flag_ok=1;
		}

		if (flag_ok==0) {
			THROW IOException("Wrong interpolation type for multiple data sources", AT);
		}
	}
}

void Interpol2D::calculate(CArray2D<double>& param_out, const vector<double>& vecExtraData, CArray2D<double>& extra_param_in) {
	unsigned short int flag_ok=0;
	vector<double> vecStationElevations;
	
	if (InputSize==0) { //no data
		if (flag_ok==0) {
			THROW IOException("Wrong interpolation type for zero data source", AT);
		}
	}
	
	if (InputSize==1) {
		//single data source
		if (flag_ok==0) {
			THROW IOException("Wrong interpolation type for single data source", AT);
		}
	}

	if (InputSize>1) {
		//multiple data sources
		if (multiple_type==I_RH) {
			//here, RH->Td, interpolations, Td->RH
			vector<double> vecTdStations(inputData.size(), 0.0); // init to 0.0

			//Compute dew point temperatures at stations
			for (unsigned int i=0; i<(unsigned int)inputData.size(); i++) {
				vecTdStations[i] = RhtoDewPoint(inputData[i],vecExtraData[i]);
			}
			
			//Krieging on Td
			BuildStationsElevations(InputMeta, vecStationElevations);
			LinRegression(vecStationElevations, vecTdStations, vecCoefficients);
			LapseRateProject = &Interpol2D::LinProject;
			LapseIDWKrieging(param_out, InputTopo, vecTdStations, InputMeta);

			//Recompute Rh from the interpolated td
			for (unsigned int i=0;i<nx;i++) {
				for (unsigned int j=0;j<ny;j++) {
					param_out(i,j) = DewPointtoRh(param_out(i,j),extra_param_in(i,j));
				}
			}
			
			flag_ok=1;
		}

		if (flag_ok==0) {
			THROW IOException("Wrong interpolation type for multiple data sources", AT);
		}
	}
}

