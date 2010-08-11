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
//This is the two 2D meteo interpolation library.
#include <meteoio/libinterpol2D.h>

using namespace std;

namespace mio {

const double Interpol2D::dflt_temperature_lapse_rate = -0.0065;
const double Interpol2D::wind_ys = 0.58;
const double Interpol2D::wind_yc = 0.42;

double Interpol2D::getReferenceAltitude(const DEMObject& dem)
{
	double ref_altitude = 1500.;

	if(dem.min_altitude!=IOUtils::nodata && dem.max_altitude!=IOUtils::nodata) {
		//we use the median elevation as the reference elevation for reprojections
		ref_altitude = 0.5 * (dem.min_altitude+dem.max_altitude); 
	} else {
		//since there is nothing else that we can do, we use an arbitrary elevation
		ref_altitude = 1500.;
	}
	return ref_altitude;
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
double Interpol2D::HorizontalDistance(const DEMObject& dem, const int& i, const int& j, const double& X2, const double& Y2)
{
	//This function computes the horizontal distance between two points
	//coordinates are given in a square, metric grid system
	//for grid points toward real coordinates
	const double X1 = (dem.llcorner.getEasting()+i*dem.cellsize);
	const double Y1 = (dem.llcorner.getNorthing()+j*dem.cellsize);
	
	return sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) );
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
void Interpol2D::LinRegressionCore(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, const int ignore_index)
{
	//finds the linear regression for points (x,y,z,Value)
	//if ignore_index>=0, ignore given index (as a way to remopve a point from the interpolation)
	double x_avg=0., y_avg=0.;
	double sx=0., sy=0., sxy=0.;
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
	a = sxy / sx;          //pbl if X=cst & Y!=cst
	b = y_avg - a*x_avg;   //pbl if X=cst & Y!=cst
	if(sy==0) {
		//horizontal line: all y's are equals
		r = 1.;
	} else {
		//any other line
		r = sxy / sqrt(sx*sy);
	}
}

/**
* @brief Computes the linear regression coefficients fitting the points given as X and Y in two vectors
* the linear regression has the form Y = aX + b with a regression coefficient r. If the regression coefficient is not good enough, a bad point is looked removed.
* @param X (vector\<double\>) vector of X coordinates
* @param Y (vector\<double\>) vector of Y coordinates (same order as X)
* @param coeffs (vector\<double\>) a,b,r coefficients in a vector
* @return (int) EXIT_SUCCESS or EXIT_FAILURE
*/
int Interpol2D::LinRegression(const std::vector<double>& X, const std::vector<double>& Y, std::vector<double>& coeffs)
{
	//finds the linear regression for points (x,y,z,Value)
	const double r_thres=0.7;
	double a,b,r;

	if ((unsigned int)X.size()==2) {
		std::cout << "[W] only two points for linear regression!" << std::endl;
	}
	if((unsigned int)X.size()<2) { //this should not be needed, we should have refrained from calling LinRegression in such a case
		std::cerr << "[E] Not enough data point for linear regression!" << std::endl;
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
		std::cout << "[W] Poor regression coefficient: " << std::setprecision(4) << coeffs[3] << std::endl;
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
* This implementation keeps the value constant as a function of the elevation.
* This interface has to follow the interface of *LapseRateProjectPtr
* @param value original value
* @param altitude altitude of the original value
* @param new_altitude altitude of the reprojected value
* @param coeffs coefficients to use for the projection
* @return reprojected value
*/
double Interpol2D::ConstProject(const double& value, const double&, const double&, const std::vector<double>&)
{
	return value;
}

/**
* @brief Projects a given parameter to another elevation: 
* This implementation assumes a linear dependency of the value as a function of the elevation.
* This interface has to follow the interface of *LapseRateProjectPtr
* @param value original value
* @param altitude altitude of the original value
* @param new_altitude altitude of the reprojected value
* @param coeffs coefficients to use for the projection
* @return reprojected value
*/
double Interpol2D::LinProject(const double& value, const double& altitude, const double& new_altitude, const std::vector<double>& coeffs)
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
* @param dem array of elevations (dem)
* @param grid 2D array to fill
*/
void Interpol2D::stdPressureGrid2DFill(const DEMObject& dem, Grid2DObject& grid) {
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);

	//provide each point with an altitude dependant pressure... it is worth what it is...
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			if (dem.grid2D(i,j)!=IOUtils::nodata) {
				grid.grid2D(i,j) = lw_AirPressure(dem.grid2D(i,j));
			} else {
				grid.grid2D(i,j) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief Grid filling function:
* This implementation fills the grid with a constant value
* @param value value to put in the grid
* @param dem array of elevations (dem). This is needed in order to know if a point is "nodata"
* @param grid 2D array to fill
*/
void Interpol2D::constantGrid2DFill(const double& value, const DEMObject& dem, Grid2DObject& grid)
{
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);

	//fills a data table with constant values
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			if (dem.grid2D(i,j)!=IOUtils::nodata) {
				grid.grid2D(i,j) = value;
			} else {
				grid.grid2D(i,j) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief Grid filling function: 
* This implementation fills a flat grid with a constant value and then reprojects it to the terrain's elevation.
* for example, the air temperature measured at one point at 1500m would be given as value, the 1500m as altitude and the dem would allow to reproject this temperature on the full DEM using the detrending function provided as pointer (with its previously calculated coefficients).
* @param value value to put in the grid
* @param altitude altitude of the "value"
* @param dem array of elevations (dem)
* @param vecCoefficients vector of detrending coefficients
* @param funcptr detrending function pointer (that uses the detrending coefficients)
* @param grid 2D array to fill
*/
void Interpol2D::constantLapseGrid2DFill(const double& value, const double& altitude, 
                                         const DEMObject& dem, const std::vector<double>& vecCoefficients,
                                         const LapseRateProjectPtr& funcptr, Grid2DObject& grid)
{
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);

	//fills a data table with constant values and then reprojects it to the DEM's elevation from a given altitude
	//the laspe rate parameters must have been set before
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			if (dem.grid2D(i,j)!=IOUtils::nodata) {
				grid.grid2D(i,j) = funcptr(value, altitude, dem.grid2D(i,j), vecCoefficients);
			} else {
				grid.grid2D(i,j) = IOUtils::nodata;
			}
		}
	}
}

double Interpol2D::IDWCore(const double& x, const double& y, const std::vector<double>& vecData_in,
                           const std::vector<StationData>& vecStations_in)
{
	//The value at any given cell is the sum of the weighted contribution from each source
	double parameter=0., norm=0., weight;
	
	for (unsigned int i=0; i<(unsigned int)vecStations_in.size(); i++) {
		weight=1./(HorizontalDistance(x, y, vecStations_in[i].position.getEasting(),
		       vecStations_in[i].position.getNorthing()) + 1e-6);
		parameter += weight*vecData_in[i];
		norm += weight;
	}
	return (parameter/norm); //normalization
}

/**
* @brief Grid filling function: 
* This implementation fills a flat grid using Inverse Distance Weighting and then reproject it to the terrain's elevation.
* for example, the air temperatures measured at several stations would be given as values, the stations altitude and positions
* as positions and projected to a flat grid. Afterward, the grid would be reprojected to the correct elevation as given
* by the dem would using the detrending function provided as pointer (with its previously calculated coefficients).
* @param vecData_in input values to use for the IDW
* @param vecStations_in position of the "values" (altitude and coordinates)
* @param dem array of elevations (dem)
* @param vecCoefficients vector of detrending coefficients
* @param funcptr detrending function pointer (that uses the detrending coefficients)
* @param grid 2D array to fill
*/
void Interpol2D::LapseIDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                          const DEMObject& dem, const std::vector<double>& vecCoefficients,
                          const LapseRateProjectPtr& funcptr,
                          Grid2DObject& grid)
{	//multiple source stations: lapse rate projection, IDW Krieging, re-projection
	const double ref_altitude = getReferenceAltitude(dem);

	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);
	std::vector<double> vecTref(vecStations_in.size(), 0.0); // init to 0.0
	
	for (unsigned int i=0; i<(unsigned int)vecStations_in.size(); i++) {
		vecTref[i] = funcptr(vecData_in[i], vecStations_in[i].position.getAltitude(), 
		                     ref_altitude, vecCoefficients);
	}

	const double xllcorner = dem.llcorner.getEasting();
	const double yllcorner = dem.llcorner.getNorthing();
	const double cellsize = dem.cellsize;
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			if (dem.grid2D(i,j)!=IOUtils::nodata) {
				grid.grid2D(i,j) = IDWCore((xllcorner+i*cellsize),
				                           (yllcorner+j*cellsize), vecTref, vecStations_in);
				grid.grid2D(i,j) = funcptr(grid.grid2D(i,j), ref_altitude,
				                           dem.grid2D(i,j), vecCoefficients);
			} else {
				grid.grid2D(i,j) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief Grid filling function: 
* This implementation fills a grid using Inverse Distance Weighting.
* for example, the air temperatures measured at several stations would be given as values, the stations positions
* as positions and projected to a grid. No elevation detrending is performed, the DEM is only used for checking if a grid point is "nodata".
* @param vecData_in input values to use for the IDW
* @param vecStations_in position of the "values" (altitude and coordinates)
* @param dem array of elevations (dem). This is needed in order to know if a point is "nodata"
* @param grid 2D array to fill
*/
void Interpol2D::IDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                     const DEMObject& dem, Grid2DObject& grid)
{
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);

	//multiple source stations: simple IDW Krieging
	const double xllcorner = dem.llcorner.getEasting();
	const double yllcorner = dem.llcorner.getNorthing();
	const double cellsize = dem.cellsize;
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			if (dem.grid2D(i,j)!=IOUtils::nodata) {
				grid.grid2D(i,j) = IDWCore((xllcorner+i*cellsize), (yllcorner+j*cellsize),
				                           vecData_in, vecStations_in);
			} else {
				grid.grid2D(i,j) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief Grid filling function: 
* This implementation fills a grid using a curvature and slope algorithm, as described in "A Meteorological
* Distribution System for High-Resolution Terrestrial Modeling (MicroMet)", Liston and Elder, 2006.
* @param dem array of elevations (dem). The slope must have been updated as it is required for the DEM analysis.
* @param VW 2D array of Wind Velocity to fill
* @param DW 2D array of Wind Direction to fill
*/
void Interpol2D::SimpleDEMWindInterpolate(const DEMObject& dem, Grid2DObject& VW, Grid2DObject& DW)
{
	if ((!VW.isSameGeolocalization(DW)) || (!VW.isSameGeolocalization(dem))){
		throw IOException("Requested grid VW and grid DW don't match the geolocalization of the DEM", AT);
	}

	//This method computes the speed of the wind and returns a table in 2D with this values
	double speed;		// Wind speed (m s-1)
	double dir;		// Wind direction
	double u;		// Zonal component u (m s-1)
	double v;		// Meridional component v (m s-1)
	double beta;		// Terrain slope
	double azi;		// Topographic slope azimuth
	double curvature;	// Topographic curvature
	double slopeDir;	// Slope in the direction of the wind
	double Ww;		// Wind weighting
	double Od;		// Diverting factor
	
	for (unsigned int j=0;j<VW.nrows-1;j++) {
		for (unsigned int i=0;i<VW.ncols-1;i++){
			// Get input data
			speed = VW.grid2D(i,j);
			dir = DW.grid2D(i,j);
			beta = dem.slope(i, j);
			azi = dem.azi(i, j);
			curvature = dem.curvature(i, j);

			if(speed==IOUtils::nodata || dir==IOUtils::nodata || beta==IOUtils::nodata || azi==IOUtils::nodata || curvature==IOUtils::nodata) {
				VW.grid2D(i, j) = IOUtils::nodata;
				DW.grid2D(i, j) = IOUtils::nodata;
			} else {
				//convert direction to rad
				dir *= ((M_PI) / 180.);
				//Speed and direction converted to zonal et meridional
				//components 
				u = (-1.) * (speed * sin(dir));
				v = (-1.) * (speed * cos(dir));

				// Converted back to speed and direction
				speed = sqrt(u*u + v*v);
				dir = (1.5 * M_PI) - atan(v/u);

				//normalize curvature and beta. 
				//Note: it should be slopeDir instead of beta, but beta is more efficient
				//to compute (only once for each dem) and it should not be that different...
				beta = (beta - dem.min_slope)/(dem.max_slope - dem.min_slope) - 0.5;
				curvature = (curvature - dem.min_curvature)/(dem.max_curvature - dem.min_curvature) - 0.5;

				// Calculate the slope in the direction of the wind
				slopeDir = beta * cos(dir - azi);
		
				// Calculate the wind weighting factor
				Ww = 1. + wind_ys * slopeDir + wind_yc * curvature;

				// Modify the wind direction by a diverting factor
				Od = -0.5 * slopeDir * sin(2.*(azi - dir));

				// Calculate the terrain-modified wind speed
				VW.grid2D(i, j) = Ww * speed;

				// Add the diverting factor to the wind direction and convert to degrees
				DW.grid2D(i, j) = (dir + Od) * (180. / (M_PI));
				if( DW.grid2D(i, j)>360. ) {
					DW.grid2D(i, j) -= 360.;
				}
			}
		}
	}
}

double Interpol2D::RhtoDewPoint(double RH, double TA, const short int& force_water)
{
	//Convert a Relative Humidity into a dew point temperature
	//TA is in Kelvins, RH between 0 and 1, returns Td in Kelvins
	TA = K_TO_C(TA);
	double Es, E, Tdw, Tdi; //saturation and current water vapor pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97;	//parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55;	//parameters for ice
	const double Tfreeze = 0.;			//freezing temperature
	const double Tnucl = -16.0;			//nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);		//distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6);	//distance to pure water

	//in order to avoid getting NaN if RH=0
	RH += 0.0001;
	assert(RH>0.);
	if (TA >= Tfreeze || force_water==1) {//above freezing point, water
		Es = Aw * exp( (Bw * TA) / (Cw + TA) );
		E = RH * Es;
		Tdw = ( Cw * log(E / Aw) ) / ( Bw - log(E / Aw) );
		return C_TO_K(Tdw);
	}
	if (TA < Tnucl) { //below nucleation, ice
		Es = Ai * exp( (Bi * TA) / (Ci + TA) );
		E = RH * Es;
		Tdi = ( Ci * log(E / Ai) ) / ( Bi - log(E / Ai) );
		return C_TO_K(Tdi);
	}

	//no clear state, we do a smooth interpolation between water and ice
	Es = Ai * exp( (Bi*TA) / (Ci + TA) );
	E = RH * Es;
	Tdi = ( Ci * log(E / Ai) ) / ( Bi - log(E / Ai) );

	Es = Aw * exp( (Bw * TA) / (Cw + TA) );
	E = RH * Es;
	Tdw = ( Cw * log(E / Aw) ) / ( Bw - log(E / Aw) );

	return C_TO_K( (di / (di + dw) * Tdi + dw / (di + dw) * Tdw) );
}

double Interpol2D::DewPointtoRh(double TD, double TA, const short int& force_water)
{
	//Convert a dew point temperature into a Relative Humidity
	//TA, TD are in Kelvins, RH is returned between 0 and 1
	TA = K_TO_C(TA);
	TD = K_TO_C(TD);
	double Es, E, Rhi, Rhw, Rh;			//saturation and current water vapro pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97;	//parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55;	//parameters for ice
	const double Tfreeze = 0.;			//freezing temperature
	const double Tnucl = -16.0;			//nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);		//distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6);	//distance to pure water

	if (TA >= Tfreeze || force_water==1) {
		//above freezing point, water
		Es = Aw * exp( (Bw * TA) / (Cw + TA) );
		E  = Aw * exp( (Bw * TD) / (Cw + TD) );
		Rhw = (E / Es);
		if (Rhw > 1.) {
			return 1.;
		} else {
			return Rhw;
		}
	}
	if (TA < Tnucl) {
		//below nucleation, ice
		Es = Ai * exp( (Bi * TA) / (Ci + TA) );
		E  = Ai * exp( (Bi * TD) / (Ci + TD) );
		Rhi = (E / Es);
		if (Rhi > 1.) {
			return 1.;
		} else {
			return Rhi;
		}
	}

	//no clear state, we do a smooth interpolation between water and ice
	Es = Ai * exp( (Bi * TA) / (Ci + TA) );
	E  = Ai * exp( (Bi * TD) / (Ci + TD) );
	Rhi = E / Es;

	Es = Aw * exp( (Bw * TA) / (Cw + TA) );
	E  = Aw * exp( (Bw * TD) / (Cw + TD) );
	Rhw = E / Es;

	Rh = (di / (di + dw) * Rhi + dw / (di + dw) * Rhw);
	if(Rh > 1.) {
		return 1.;
	} else {
		return Rh;
	}
}

double Interpol2D::lw_AirPressure(const double& altitude)
{
	double p;
	const double p0 = 101325.; 		// Air and standard pressure in Pa
	const double lapse_rate = 0.0065;	// K m-1
	const double sea_level_temp = 288.15;	// K
	const double expo = GRAVITY / (lapse_rate * GAS_CONSTANT_AIR);
	const double R0 = 6356766.0;		// Earth's radius in m
	
	p = p0 * pow( 1. - ( (lapse_rate * R0 * altitude) / (sea_level_temp * (R0 + altitude)) ), expo );
	
	return(p);
}

} //namespace
