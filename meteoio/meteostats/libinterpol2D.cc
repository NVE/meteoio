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
#include <meteoio/meteolaws/Meteoconst.h> //for math constants
#include <meteoio/MathOptim.h> //math optimizations

using namespace std;

namespace mio {
const double Interpol2D::wind_ys = 0.58;
const double Interpol2D::wind_yc = 0.42;

//Usefull functions
/**
* @brief Computes the horizontal distance between points, given by coordinates in a geographic grid
* @param X1 (const double) first point's X coordinate
* @param Y1 (const double) first point's Y coordinate
* @param X2 (const double) second point's X coordinate
* @param Y2 (const double) second point's Y coordinate
* @return (double) distance in m
*/
inline double Interpol2D::HorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2)
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
inline double Interpol2D::InvHorizontalDistance(const double& X1, const double& Y1, const double& X2, const double& Y2)
{
	//This function computes 1/horizontaldistance between two points
	//coordinates are given in a square, metric grid system
	const double DX=(X1-X2), DY=(Y1-Y2);
	return Optim::invSqrt( DX*DX + DY*DY ); //we use the optimized approximation for 1/sqrt
}

/**
* @brief Computes the horizontal distance between points, given by their cells indexes
* @param X1 (const double) first point's i index
* @param Y1 (const double) first point's j index
* @param X2 (const double) second point's X coordinate
* @param Y2 (const double) second point's Y coordinate
* @return (double) distance in m
*/
inline double Interpol2D::HorizontalDistance(const DEMObject& dem, const int& i, const int& j, const double& X2, const double& Y2)
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
* @brief Build the list of (distance to grid cell, stations index) ordered by their distance to a grid cell
* @param x x coordinate of cell
* @param y y coordinate of cell
* @param list list of pairs (distance to grid cell, stations index)
*/
void Interpol2D::getNeighbors(const double& x, const double& y,
                              const std::vector<StationData>& vecStations,
                              std::vector< std::pair<double, size_t> >& list)
{
	list.resize(vecStations.size());

	for(size_t i=0; i<vecStations.size(); i++) {
		const Coords& position = vecStations[i].position;
		const double DX = x-position.getEasting();
		const double DY = y-position.getNorthing();
		const double d2 = (DX*DX + DY*DY);
		const std::pair <double, size_t> tmp(d2,i);
		list[i] = tmp;
	}

	sort (list.begin(), list.end());
}

//convert a vector of stations into two vectors of eastings and northings
void Interpol2D::buildPositionsVectors(const std::vector<StationData>& vecStations, std::vector<double>& vecEastings, std::vector<double>& vecNorthings)
{
	vecEastings.resize( vecStations.size() );
	vecNorthings.resize( vecStations.size() );
	for (size_t i=0; i<vecStations.size(); i++) {
		const Coords& position = vecStations[i].position;
		vecEastings[i] = position.getEasting();
		vecNorthings[i] = position.getNorthing();
	}
}

//these weighting functions take the square of a distance as an argument and return a weight
inline double Interpol2D::weightInvDist(const double& d2)
{
	return Optim::invSqrt( d2 ); //we use the optimized approximation for 1/sqrt
}
inline double Interpol2D::weightInvDistSqrt(const double& d2)
{
	return Optim::fastSqrt_Q3( Optim::invSqrt(d2) ); //we use the optimized approximation for 1/sqrt
}
inline double Interpol2D::weightInvDist2(const double& d2)
{
	return 1./d2; //we use the optimized approximation for 1/sqrt
}
inline double Interpol2D::weightInvDistN(const double& d2)
{
	return pow( Optim::invSqrt(d2) , dist_pow); //we use the optimized approximation for 1/sqrt
}

//Filling Functions
/**
* @brief Grid filling function:
* This implementation builds a standard air pressure as a function of the elevation
* @param dem array of elevations (dem)
* @param grid 2D array to fill
*/
void Interpol2D::stdPressure(const DEMObject& dem, Grid2DObject& grid) {
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
void Interpol2D::constant(const double& value, const DEMObject& dem, Grid2DObject& grid)
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

double Interpol2D::IDWCore(const double& x, const double& y, const std::vector<double>& vecData_in,
                           const std::vector<double>& vecEastings, const std::vector<double>& vecNorthings)
{
	//The value at any given cell is the sum of the weighted contribution from each source
	const size_t n_stations=vecEastings.size();
	double parameter=0., norm=0.;
	const double scale = 1.e6;

	for (size_t i=0; i<n_stations; i++) {
		const double DX=x-vecEastings[i];
		const double DY=y-vecNorthings[i];
		const double weight = Optim::invSqrt( DX*DX + DY*DY + scale ); //use the optimized 1/sqrt approximation
		parameter += weight*vecData_in[i];
		norm += weight;
	}
	return (parameter/norm); //normalization
}

/*
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
/*void Interpol2D::LocalLapseIDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                               const DEMObject& dem, const size_t& nrOfNeighbors,
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
			const double value = LLIDW_pixel(i,j,vecData_in, vecStations_in, dem, nrOfNeighbors, r); //TODO: precompute in vectors
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
                                const DEMObject& dem, const size_t& nrOfNeighbors, double& r2)
{
	const double& cell_altitude=dem.grid2D(i,j);
	if(cell_altitude==IOUtils::nodata)
		return IOUtils::nodata;

	std::vector< std::pair<double, size_t> > list;
	std::vector<double> X, Y, coeffs;

	//fill vectors with appropriate neighbors
	const double x = dem.llcorner.getEasting()+i*dem.cellsize;
	const double y = dem.llcorner.getNorthing()+j*dem.cellsize;
	getNeighbors(x, y, vecStations_in, list);
	const size_t max_stations = std::min(list.size(), nrOfNeighbors);
	for(size_t st=0; st<max_stations; st++) {
		const size_t st_index=list[st].second;
		const double value = vecData_in[st_index];
		const double alt = vecStations_in[st_index].position.getAltitude();
		if ((value != IOUtils::nodata) && (alt != IOUtils::nodata)) {
			X.push_back( alt );
			Y.push_back( value );
		}
	}

	//compute lapse rate
	if(X.empty())
		return IOUtils::nodata;
	coeffs.resize(7,0.);
	//BiLinRegression(X, Y, coeffs);
	r2=coeffs[3]*coeffs[6]; //Is it correct?

	//compute local pixel value
	unsigned int count=0;
	double pixel_value=0., norm=0.;
	const double scale=0.;
	for(size_t st=0; st<max_stations; st++) {
		const size_t st_index=list[st].second;
		const double value = vecData_in[st_index];
		const double alt = vecStations_in[st_index].position.getAltitude();
		if ((value != IOUtils::nodata) && (alt != IOUtils::nodata)) {
			//const double contrib = LinProject(value, alt, cell_altitude, coeffs);
			//const double contrib = BiLinProject(value, alt, cell_altitude, coeffs);
			const double weight = Optim::invSqrt( list[st].first + scale + 1.e-6 );
			//pixel_value += weight*contrib;
			norm += weight;
			count++;
		}
	}

	if(count>0)
		return (pixel_value/norm);
	else
		return IOUtils::nodata;
}*/

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
	std::vector<double> vecEastings, vecNorthings;
	buildPositionsVectors(vecStations_in, vecEastings, vecNorthings);

	//multiple source stations: simple IDW Krieging
	const double xllcorner = dem.llcorner.getEasting();
	const double yllcorner = dem.llcorner.getNorthing();
	const double cellsize = dem.cellsize;
	for (unsigned int j=0; j<grid.nrows; j++) {
		for (unsigned int i=0; i<grid.ncols; i++) {
			if (dem.grid2D(i,j)!=IOUtils::nodata) {
				grid.grid2D(i,j) = IDWCore((xllcorner+i*cellsize), (yllcorner+j*cellsize),
				                           vecData_in, vecEastings, vecNorthings);
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
* @param i_dem array of elevations (dem). The slope must have been updated as it is required for the DEM analysis.
* @param VW 2D array of Wind Velocity to fill
* @param DW 2D array of Wind Direction to fill
*/
void Interpol2D::SimpleDEMWindInterpolate(const DEMObject& i_dem, Grid2DObject& VW, Grid2DObject& DW)
{
	if ((!VW.isSameGeolocalization(DW)) || (!VW.isSameGeolocalization(i_dem))){
		throw IOException("Requested grid VW and grid DW don't match the geolocalization of the DEM", AT);
	}

	const bool recomputeDEM = i_dem.curvature.isEmpty();
	DEMObject *intern_dem = NULL;
	if(recomputeDEM) {
		std::cerr << "[W] WIND_CURV spatial interpolations algorithm selected but no dem curvature available! Computing it...\n";
		intern_dem = new DEMObject(i_dem);
		intern_dem->setUpdatePpt((DEMObject::update_type)(DEMObject::SLOPE|DEMObject::CURVATURE));
		intern_dem->update();
	}
	const DEMObject *dem = (recomputeDEM)? intern_dem : &i_dem;

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

	const double dem_min_slope=dem->min_slope*Cst::to_rad;
	const double dem_min_curvature=dem->min_curvature;
	double dem_range_slope=(dem->max_slope-dem_min_slope)*Cst::to_rad;
	double dem_range_curvature=(dem->max_curvature-dem_min_curvature);
	if(dem_range_slope==0.) dem_range_slope = 1.; //to avoid division by zero below
	if(dem_range_curvature==0.) dem_range_curvature = 1.; //to avoid division by zero below

	for (unsigned int j=0;j<VW.nrows;j++) {
		for (unsigned int i=0;i<VW.ncols;i++){
			speed = VW.grid2D(i,j);
			if(speed==0.) continue; //we can not apply any correction factor!
			dir = DW.grid2D(i,j);
			beta = dem->slope(i, j)*Cst::to_rad;
			azi = dem->azi(i, j)*Cst::to_rad;
			curvature = dem->curvature(i, j);

			if(speed==IOUtils::nodata || dir==IOUtils::nodata || beta==IOUtils::nodata || azi==IOUtils::nodata || curvature==IOUtils::nodata) {
				VW.grid2D(i, j) = IOUtils::nodata;
				DW.grid2D(i, j) = IOUtils::nodata;
			} else {
				//convert direction to rad
				dir *= Cst::to_rad;
				//Speed and direction converted to zonal et meridional
				//components
				u = (-1.) * (speed * sin(dir));
				v = (-1.) * (speed * cos(dir));

				// Converted back to speed and direction
				speed = sqrt(u*u + v*v);
				dir = (1.5 * Cst::PI) - atan(v/u);

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
				DW.grid2D(i, j) = (dir + Od) * Cst::to_deg;
				if( DW.grid2D(i, j)>360. ) {
					DW.grid2D(i, j) -= 360.;
				}
			}
		}
	}

	if(intern_dem!=NULL) delete (intern_dem);
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
	//precompute various coordinates in the grid
	const double llcorner_x = grid.llcorner.getEasting();
	const double llcorner_y = grid.llcorner.getNorthing();
	const double cellsize = grid.cellsize;

	Matrix Ginv(nrOfMeasurments+1, nrOfMeasurments+1);

	//fill the Ginv matrix
	for(size_t j=1; j<=nrOfMeasurments; j++) {
		const Coords& st1 = vecStations[j-1].position;
		const double x1 = st1.getEasting();
		const double y1 = st1.getNorthing();

		for(size_t i=1; i<=j; i++) {
			//compute distance between stations
			const Coords& st2 = vecStations[i-1].position;
			const double DX = x1-st2.getEasting();
			const double DY = y1-st2.getNorthing();
			const double distance = Optim::fastSqrt_Q3(DX*DX + DY*DY);
			Ginv(i,j) = variogram.f(distance);
		}
		Ginv(j,j)=1.; //HACK diagonal should contain the nugget...
		Ginv(nrOfMeasurments+1,j) = 1.; //last line filled with 1s
	}
	//fill the upper half (an exact copy of the lower half)
	for(size_t j=1; j<=nrOfMeasurments; j++) {
		for(size_t i=j+1; i<=nrOfMeasurments; i++) {
			Ginv(i,j) = Ginv(j,i);
		}
	}
	//add last column of 1's and a zero
	for(size_t i=1; i<=nrOfMeasurments; i++) Ginv(i,nrOfMeasurments+1) = 1.;
	Ginv(nrOfMeasurments+1,nrOfMeasurments+1) = 0.;
	//invert the matrix
	Ginv.inv();

	Matrix G0(nrOfMeasurments+1, (size_t)1);
	//now, calculate each point
	for(size_t j=0; j<grid.nrows; j++) {
		for(size_t i=0; i<grid.ncols; i++) {
			const double x = llcorner_x+static_cast<double>(i)*cellsize;
			const double y = llcorner_y+static_cast<double>(j)*cellsize;

			//fill gamma
			for(size_t st=0; st<nrOfMeasurments; st++) {
				//compute distance between cell and each station
				const Coords& position = vecStations[st].position;
				const double DX = x-position.getEasting();
				const double DY = y-position.getNorthing();
				const double distance = Optim::fastSqrt_Q3(DX*DX + DY*DY);

				G0(st+1,1) = variogram.f(distance); //matrix starts at 1, not 0
			}
			G0(nrOfMeasurments+1,1) = 1.; //last value is always 1

			const Matrix lambda = Ginv*G0;

			//calculate local parameter interpolation
			double p = 0.;
			for(size_t st=0; st<nrOfMeasurments; st++) {
				p += lambda(st+1,1) * vecData[st]; //matrix starts at 1, not 0
			}
			grid.grid2D(i,j) = p;
		}
	}
}

} //namespace
