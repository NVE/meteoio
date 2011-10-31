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
#include <cmath>

#include <meteoio/meteostats/libinterpol2D.h>
#include <meteoio/meteolaws/Atmosphere.h>

using namespace std;

namespace mio {

//Quake3 fast 1/x² approximation
// For Magic Derivation see: Chris Lomont http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf
// Credited to Greg Walsh.
// 32  Bit float magic number
#define SQRT_MAGIC_D 0x5f3759df
#define SQRT_MAGIC_F 0x5f375a86

#ifdef _MSC_VER
#pragma warning( push ) //for Visual C++
#pragma warning(disable:4244) //Visual C++ righhtfully complains... but this behavior is what we want!
#endif
//maximum relative error is <1.7% while computation time for sqrt is <1/4. At 0, returns a large number
//on a large scale interpolation test on TA, max relative error is 1e-6
inline float invSqrt(const float x) {
	const float xhalf = 0.5f*x;

	union {
		// get bits for floating value
		float x;
		int i;
	} u;
	u.x = x;
	u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
	return u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
}

inline double invSqrt(const double x) {
	const double xhalf = 0.5f*x;

	union {
		// get bits for floating value
		float x;
		int i;
	} u;
	u.x = x;
	u.i = SQRT_MAGIC_D - (u.i >> 1);  // gives initial guess y0
	return u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
}
#ifdef _MSC_VER
#pragma warning( pop ) //for Visual C++, restore previous warnings behavior
#endif

inline float fastSqrt_Q3(const float x) {
	return x * invSqrt(x);
}

inline double fastSqrt_Q3(const double x) {
	return x * invSqrt(x);
}

const double Interpol2D::wind_ys = 0.58;
const double Interpol2D::wind_yc = 0.42;
const double Interpol2D::bilin_inflection = 1200.;

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
	const double DX=(X1-X2), DY=(Y1-Y2);
	return sqrt( DX*DX + DY*DY );
}

/**
* @brief Computes the 1/horizontal distance between points, given by coordinates in a geographic grid
* @param X1 (const double) first point's X coordinate
* @param Y1 (const double) first point's Y coordinate
* @param X2 (const double) second point's X coordinate
* @param Y2 (const double) second point's Y coordinate
* @return (double) 1/distance in m
*/
double Interpol2D::InvHorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2)
{
	//This function computes 1/horizontaldistance between two points
	//coordinates are given in a square, metric grid system
	const double DX=(X1-X2), DY=(Y1-Y2);
	return invSqrt( DX*DX + DY*DY ); //we use the optimized approximation for 1/sqrt
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
{//TODO: store DEM easting/northing, etc as private members
	//This function computes the horizontal distance between two points
	//coordinates are given in a square, metric grid system
	//for grid points toward real coordinates
	const double X1 = (dem.llcorner.getEasting()+i*dem.cellsize);
	const double Y1 = (dem.llcorner.getNorthing()+j*dem.cellsize);
	const double DX=(X1-X2), DY=(Y1-Y2);
	return sqrt( DX*DX + DY*DY );
}

/**
* @brief Build the list of (stations index, distance to grid cell) ordered by their distance to a grid cell
* @param x x coordinate of cell
* @param y y coordinate of cell
* @param list list of pairs (stations index, distance to grid cell)
*/
void Interpol2D::getNeighbors(const double& x, const double& y,
                              const std::vector<StationData>& vecStations,
                              std::vector< std::pair<double, unsigned int> >& list)
{
	if(list.size()>0) list.clear();

	for(unsigned int i=0; i<vecStations.size(); i++) {
		const Coords& position = vecStations[i].position;
		const double DX = x-position.getEasting();
		const double DY = y-position.getNorthing();
		const double d2 = (DX*DX + DY*DY);
		const std::pair <double, unsigned int> tmp(d2,i);
		list.push_back(tmp);
	}

	sort (list.begin(), list.end());
}

//these weighting functions take the square of a distance as an argument and return a weight
double Interpol2D::weightInvDist(const double& d2)
{
	return invSqrt( d2 ); //we use the optimized approximation for 1/sqrt
}
double Interpol2D::weightInvDistSqrt(const double& d2)
{
	return fastSqrt_Q3( invSqrt(d2) ); //we use the optimized approximation for 1/sqrt
}
double Interpol2D::weightInvDist2(const double& d2)
{
	return 1./d2; //we use the optimized approximation for 1/sqrt
}
double Interpol2D::weightInvDistN(const double& d2)
{
	return pow( invSqrt(d2) , dist_pow); //we use the optimized approximation for 1/sqrt
}

//Data regression models
/**
* @brief Computes the linear regression coefficients fitting the points given as X and Y in two vectors
* It relies on Interpol1D::NoisyLinRegression.
* @param in_X vector of X coordinates
* @param in_Y vector of Y coordinates (same order as X)
* @param coeffs a,b,r coefficients in a vector
* @return EXIT_SUCCESS or EXIT_FAILURE
*/
int Interpol2D::LinRegression(const std::vector<double>& in_X, const std::vector<double>& in_Y, std::vector<double>& coeffs)
{
	std::stringstream mesg;
	const int code = Interpol1D::NoisyLinRegression(in_X, in_Y, coeffs[1], coeffs[2], coeffs[3], mesg);
	cout << mesg.str();
	return code;
}

//temporary solution while we migrate the regression classes
int Interpol2D::BiLinRegression(const std::vector<double>& in_X, const std::vector<double>& in_Y, std::vector<double>& coeffs)
{
	std::vector<double> params;
	const int code = Interpol1D::twoLinRegression(in_X, in_Y, Interpol2D::bilin_inflection, params);
	coeffs[1] = params[1]; coeffs[2] = params[2]; coeffs[3] = 1.;
	coeffs[4] = params[3]; coeffs[5] = params[4]; coeffs[6] = 1.;
	return code;
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

/**
* @brief Projects a given parameter to another elevation:
* This implementation assumes a 2 segments linear dependency of the value as a function of the elevation.
* It uses Interpol2D::bilin_inflection as the inflection point altitude. This interface has to follow the interface of *LapseRateProjectPtr
* @param value original value
* @param altitude altitude of the original value
* @param new_altitude altitude of the reprojected value
* @param coeffs coefficients to use for the projection
* @return reprojected value
*/
double Interpol2D::BiLinProject(const double& value, const double& altitude, const double& new_altitude, const std::vector<double>& coeffs)
{
	//linear lapse: coeffs must have been already computed
	if (coeffs.size()<1) {
		throw IOException("Linear regression coefficients not initialized", AT);
	}
	if(altitude<=bilin_inflection)
		return (value + coeffs[1] * (new_altitude - altitude));
	else
		return (value + coeffs[4] * (new_altitude - altitude));
}

/**
* @brief Projects a given parameter to another elevation:
* This implementation assumes that the value increases by a given fraction as a function of the elevation.
* This interface has to follow the interface of *LapseRateProjectPtr
* @param value original value
* @param altitude altitude of the original value
* @param new_altitude altitude of the reprojected value
* @param coeffs coefficients to use for the projection
* @return reprojected value
*/
double Interpol2D::FracProject(const double& value, const double& altitude, const double& new_altitude, const std::vector<double>& coeffs)
{
	//linear lapse: coeffs must have been already computed
	if (coeffs.size()<1) {
		throw IOException("Linear regression coefficients not initialized", AT);
	}
	return (value * (1. + coeffs[1] * (new_altitude - altitude)));
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
			const double& cell_altitude=dem.grid2D(i,j);
			if (cell_altitude!=IOUtils::nodata) {
				grid.grid2D(i,j) = Atmosphere::stdAirPressure(cell_altitude);
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
			const double cell_altitude=dem.grid2D(i,j);
			if (cell_altitude!=IOUtils::nodata) {
				grid.grid2D(i,j) = funcptr(value, altitude, cell_altitude, vecCoefficients);
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
	const size_t n_stations=vecStations_in.size();
	double parameter=0., norm=0.;
	const double scale = 1.e6;

	for (size_t i=0; i<n_stations; i++) {
		/*const double weight=1./(HorizontalDistance(x, y, vecStations_in[i].position.getEasting(),
		       vecStations_in[i].position.getNorthing()) + 1e-6);*/
		const Coords& position = vecStations_in[i].position;
		const double DX=x-position.getEasting(); //optimization: precompute and store in a vector?
		const double DY=y-position.getNorthing();
		const double weight = invSqrt( DX*DX + DY*DY + scale ); //use the optimized 1/sqrt approximation
		//const double weight = weightInvDistSqrt( (DX*DX + DY*DY) );
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
	const size_t n_stations=vecStations_in.size();

	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);
	std::vector<double> vecTref(vecStations_in.size(), 0.0); // init to 0.0

	for (size_t i=0; i<n_stations; i++) {
		vecTref[i] = funcptr(vecData_in[i], vecStations_in[i].position.getAltitude(),
		                     ref_altitude, vecCoefficients);
	}

	const double xllcorner = dem.llcorner.getEasting();
	const double yllcorner = dem.llcorner.getNorthing();
	const double cellsize = dem.cellsize;
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			const double cell_altitude=dem.grid2D(i,j);
			if (cell_altitude!=IOUtils::nodata) {
				grid.grid2D(i,j) = IDWCore((xllcorner+i*cellsize),
				                           (yllcorner+j*cellsize), vecTref, vecStations_in);
				grid.grid2D(i,j) = funcptr(grid.grid2D(i,j), ref_altitude,
				                           cell_altitude, vecCoefficients);
			} else {
				grid.grid2D(i,j) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief Grid filling function:
* Similar to Interpol2D::LapseIDW but using a limited number of stations for each cell. We also assume a two segments regression for altitude detrending with
* a fixed 1200m above sea level inflection point.
* @param vecData_in input values to use for the IDW
* @param vecStations_in position of the "values" (altitude and coordinates)
* @param dem array of elevations (dem)
* @param nrOfNeighbors number of neighboring stations to use for each pixel
* @param grid 2D array to fill
* @param r2 average r² coefficient of the lapse rate regressions
*/
void Interpol2D::LocalLapseIDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                               const DEMObject& dem, const unsigned int& nrOfNeighbors,
                               Grid2DObject& grid, double& r2)
{
	unsigned int count=0;
	double sum=0;
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);

	//run algorithm
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			//LL_IDW_pixel returns nodata when appropriate
			double r;
			const double value = LLIDW_pixel(i,j,vecData_in, vecStations_in, dem, nrOfNeighbors, r);
			grid.grid2D(i,j) = value;
			if(value!=IOUtils::nodata) {
				sum += fabs(r);
				count++;
			}
		}
	}
	if(count>0)
		r2 = sum/(double)count;
	else
		r2 = 0.;
}

//calculate a local pixel for LocalLapseIDW
double Interpol2D::LLIDW_pixel(const unsigned int& i, const unsigned int& j,
                                const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                                const DEMObject& dem, const unsigned int& nrOfNeighbors, double& r2)
{
	const double& cell_altitude=dem.grid2D(i,j);
	if(cell_altitude==IOUtils::nodata)
		return IOUtils::nodata;

	std::vector< std::pair<double, unsigned int> > list;
	std::vector<double> X, Y, coeffs;

	//fill vectors with appropriate neighbors
	const double x = dem.llcorner.getEasting()+i*dem.cellsize;
	const double y = dem.llcorner.getNorthing()+j*dem.cellsize;
	getNeighbors(x, y, vecStations_in, list);
	for(unsigned int st=0; st<nrOfNeighbors; st++) {
		const unsigned int st_index=list[st].second;
		const double value = vecData_in[st_index];
		const double alt = vecStations_in[st_index].position.getAltitude();
		if ((value != IOUtils::nodata) && (alt != IOUtils::nodata)) {
			X.push_back( alt );
			Y.push_back( value );
		}
	}

	//compute lapse rate
	if(X.size()==0)
		return IOUtils::nodata;
	std::stringstream mesg;
	//coeffs.resize(4, 0.0);
	//Interpol1D::NoisyLinRegression(X, Y, coeffs[1], coeffs[2], coeffs[3], mesg);
	coeffs.resize(7,0.);
	BiLinRegression(X, Y, coeffs);
	r2=coeffs[3]*coeffs[6]; //Is it correct?

	//compute local pixel value
	unsigned int count=0;
	double pixel_value=0., norm=0.;
	const double scale=0.;
	for(unsigned int st=0; st<nrOfNeighbors; st++) {
		const unsigned int st_index=list[st].second;
		const double value = vecData_in[st_index];
		const double alt = vecStations_in[st_index].position.getAltitude();
		if ((value != IOUtils::nodata) && (alt != IOUtils::nodata)) {
			//const double contrib = LinProject(value, alt, cell_altitude, coeffs);
			const double contrib = BiLinProject(value, alt, cell_altitude, coeffs);
			const double weight = invSqrt( list[st].first + scale + 1.e-6 );
			pixel_value += weight*contrib;
			norm += weight;
			count++;
		}
	}

	if(count>0)
		return (pixel_value/norm);
	else
		return IOUtils::nodata;
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

	const double to_rad = M_PI/180.;
	const double dem_min_slope=dem.min_slope*to_rad, dem_range_slope=(dem.max_slope-dem_min_slope)*to_rad;
	const double dem_min_curvature=dem.min_curvature, dem_range_curvature=(dem.max_curvature-dem_min_curvature);

	for (unsigned int j=0;j<VW.nrows-1;j++) {
		for (unsigned int i=0;i<VW.ncols-1;i++){
			// Get input data
			speed = VW.grid2D(i,j);
			dir = DW.grid2D(i,j);
			beta = dem.slope(i, j)*to_rad;
			azi = dem.azi(i, j)*to_rad;
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
				beta = (beta - dem_min_slope)/dem_range_slope - 0.5;
				curvature = (curvature - dem_min_curvature)/dem_range_curvature - 0.5;

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

/**
* @brief Distribute precipitation in a way that reflects snow redistribution on the ground, according to (Huss, 2008)
* This method modifies the solid precipitation distribution according to the local slope and curvature. See
* <i>"Quantitative evaluation of different hydrological modelling approaches in a partly glacierized Swiss watershed"</i>, Magnusson et All., Hydrological Processes, 2010, under review.
* and
* <i>"Modelling runoff from highly glacierized alpine catchments in a changing climate"</i>, Huss et All., Hydrological Processes, <b>22</b>, 3888-3902, 2008.
* @param dem array of elevations (dem). The slope must have been updated as it is required for the DEM analysis.
* @param ta array of air temperatures used to determine if precipitation is rain or snow
* @param grid 2D array of precipitation to fill
* @author Florian Kobierska, Jan Magnusson and Mathias Bavay
*/
void Interpol2D::PrecipSnow(const DEMObject& dem, const Grid2DObject& ta, Grid2DObject& grid)
{
	if(!grid.isSameGeolocalization(dem)) {
		throw IOException("Requested grid does not match the geolocalization of the DEM", AT);
	}
	const double dem_max_curvature=dem.max_curvature, dem_range_curvature=(dem.max_curvature-dem.min_curvature);

	for (unsigned int j=0;j<grid.nrows;j++) {
		for (unsigned int i=0;i<grid.ncols;i++) {
			// Get input data
			const double slope = dem.slope(i, j);
			const double curvature = dem.curvature(i, j);
			double val = grid.grid2D(i, j);

			if(ta.grid2D(i, j)<=273.15) {
				//we only modify the grid of precipitations if air temperature
				//at this point is below or at freezing
				if(slope==IOUtils::nodata || curvature==IOUtils::nodata) {
					val = IOUtils::nodata;
				} else if (slope>60.) {
					//No snow precipitation happens for these slopes
					val = 0.;
				} else if (slope>40.) {
					//Linear transition from no snow to 100% snow
					val *= (60.-slope)/20.;
				} //else: unchanged

				if(val!=IOUtils::nodata && dem_range_curvature!=0.) {
					//cf Huss
					grid.grid2D(i, j) = val*(0.5-(curvature-dem_max_curvature)/dem_range_curvature);
				}
			}
		}
	}
}

/**
* @brief Ordinary Kriging matrix formulation
* This implements the matrix formulation of Ordinary Kriging, as shown (for example) in
* <i>"Statistics for spatial data"</i>, Noel A. C. Cressie, John Wiley & Sons, revised edition, 1993, pp122.
* @param vecData vector containing the values as measured at the stations
* @param vecStations vector of stations
* @param dem digital elevation model
* @param variogram variogram regression model
* @param grid 2D array of precipitation to fill
* @author Mathias Bavay
*/
void Interpol2D::ODKriging(const std::vector<double>& vecData, const std::vector<StationData>& vecStations, const DEMObject& dem, const Fit1D& variogram, Grid2DObject& grid)
{
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);
	size_t nrOfMeasurments = vecStations.size();

	Matrix G((unsigned int)nrOfMeasurments, (unsigned int)nrOfMeasurments);
	Matrix gamma(nrOfMeasurments, (unsigned int)1);
	const Matrix One(nrOfMeasurments, (unsigned int)1, 1.);
	const Matrix One_T = One.getT();

	//precompute various coordinates in the grid
	const double llcorner_x = grid.llcorner.getEasting();
	const double llcorner_y = grid.llcorner.getNorthing();
	const double cellsize = grid.cellsize;

	//fill the G matrix
	//HACK: are we filling with the proper values? A covariance matrix would be different...
	for(size_t j=1; j<=nrOfMeasurments; j++) {
		const Coords& st1 = vecStations[j-1].position;
		const double x1 = st1.getEasting();
		const double y1 = st1.getNorthing();

		for(size_t i=1; i<=j; i++) {
			//compute distance between stations
			const Coords& st2 = vecStations[i-1].position;
			const double DX = x1-st2.getEasting();
			const double DY = y1-st2.getNorthing();
			const double distance = fastSqrt_Q3(DX*DX + DY*DY);
			G(i,j) = variogram.f(distance);
		}
		//G(j,j)=1.; //HACK what should we put on the diagonal?
	}
	//fill the upper half (an exact copy of the lower half)
	for(size_t j=1; j<=nrOfMeasurments; j++) {
		for(size_t i=j+1; i<=nrOfMeasurments; i++) {
			G(i,j) = G(j,i);
		}
	}

	//G inverse matrix
	const Matrix Ginv = G.getInv();

	//calculate constant denominator
	const Matrix OneT_Ginv = One_T * Ginv;
	const double denom = Matrix::scalar( OneT_Ginv * One );

	//now, calculate each point
	for(size_t j=0; j<grid.nrows; j++) {
		for(size_t i=0; i<grid.ncols; i++) {
			const double x = llcorner_x+i*cellsize;
			const double y = llcorner_y+j*cellsize;

			//fill gamma
			for(size_t st=0; st<nrOfMeasurments; st++) {
				//compute distance between cell and each station
				const Coords& position = vecStations[st].position;
				const double DX = x-position.getEasting();
				const double DY = y-position.getNorthing();
				const double distance = fastSqrt_Q3(DX*DX + DY*DY);

				gamma(st+1,1) = variogram.f(distance); //matrix starts at 1
			}

			const Matrix lambdaT = Matrix::T(gamma + One * ((1. - Matrix::scalar(OneT_Ginv*gamma)) / denom) ) * Ginv;

			//calculate local parameter interpolation
			double p = 0.;
			for(size_t st=0; st<nrOfMeasurments; st++) {
				p += lambdaT(1,st+1) * vecData[st]; //matrix starts at 1
			}
			grid.grid2D(i,j) = p;
		}
	}
}

} //namespace
