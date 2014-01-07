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
#include <algorithm>

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
{
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
	for (size_t j=0; j<grid.nrows; j++) {
		for (size_t i=0; i<grid.ncols; i++) {
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
	for (size_t j=0; j<grid.nrows; j++) {
		for (size_t i=0; i<grid.ncols; i++) {
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

/** @brief Grid filling function:
* Similar to Interpol2D::LapseIDW but using a limited number of stations for each cell.
* @param vecData_in input values to use for the IDW
* @param vecStations_in position of the "values" (altitude and coordinates)
* @param dem array of elevations (dem)
* @param nrOfNeighbors number of neighboring stations to use for each pixel
* @param grid 2D array to fill
*/
void Interpol2D::LocalLapseIDW(const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                               const DEMObject& dem, const size_t& nrOfNeighbors,
                               Grid2DObject& grid)
{
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);

	//run algorithm
	for (size_t j=0; j<grid.nrows; j++) {
		for (size_t i=0; i<grid.ncols; i++) {
			//LL_IDW_pixel returns nodata when appropriate
			grid.grid2D(i,j) = LLIDW_pixel(i, j, vecData_in, vecStations_in, dem, nrOfNeighbors); //TODO: precompute in vectors
		}
	}
}

//calculate a local pixel for LocalLapseIDW
double Interpol2D::LLIDW_pixel(const size_t& i, const size_t& j,
                               const std::vector<double>& vecData_in, const std::vector<StationData>& vecStations_in,
                               const DEMObject& dem, const size_t& nrOfNeighbors)
{
	const double cell_altitude = dem.grid2D(i,j);
	if(cell_altitude==IOUtils::nodata)
		return IOUtils::nodata;

	std::vector< std::pair<double, size_t> > list;
	std::vector<double> X, Y;

	//fill vectors with appropriate neighbors
	const double x = dem.llcorner.getEasting()+i*dem.cellsize;
	const double y = dem.llcorner.getNorthing()+j*dem.cellsize;
	getNeighbors(x, y, vecStations_in, list);
	const size_t max_stations = std::min(list.size(), nrOfNeighbors);
	for(size_t st=0; st<max_stations; st++) {
		const size_t st_index = list[st].second;
		const double value = vecData_in[st_index];
		const double alt = vecStations_in[st_index].position.getAltitude();
		if ((value != IOUtils::nodata) && (alt != IOUtils::nodata)) {
			X.push_back( alt );
			Y.push_back( value );
		}
	}

	//compute lapse rate
	if(X.empty()) return IOUtils::nodata;
	const Fit1D trend(Fit1D::NOISY_LINEAR, X, Y);

	//compute local pixel value
	unsigned int count=0;
	double pixel_value=0., norm=0.;
	const double scale=0.;
	for(size_t st=0; st<max_stations; st++) {
		const size_t st_index=list[st].second;
		const double alt = vecStations_in[st_index].position.getAltitude();
		const double value = vecData_in[st_index];
		if ((value != IOUtils::nodata) && (alt != IOUtils::nodata)) {
			const double contrib = value - trend(alt);
			const double weight = Optim::invSqrt( list[st].first + scale + 1.e-6 );
			pixel_value += weight*contrib;
			norm += weight;
			count++;
		}
	}

	if(count>0)
		return (pixel_value/norm) + trend(cell_altitude);
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
	std::vector<double> vecEastings, vecNorthings;
	buildPositionsVectors(vecStations_in, vecEastings, vecNorthings);

	//multiple source stations: simple IDW Krieging
	const double xllcorner = dem.llcorner.getEasting();
	const double yllcorner = dem.llcorner.getNorthing();
	const double cellsize = dem.cellsize;
	for (size_t jj=0; jj<grid.nrows; jj++) {
		for (size_t ii=0; ii<grid.ncols; ii++) {
			if (dem.grid2D(ii,jj)!=IOUtils::nodata) {
				grid.grid2D(ii,jj) = IDWCore((xllcorner+double(ii)*cellsize), (yllcorner+double(jj)*cellsize),
				                           vecData_in, vecEastings, vecNorthings);
			} else {
				grid.grid2D(ii,jj) = IOUtils::nodata;
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

	for (size_t j=0;j<VW.nrows;j++) {
		for (size_t i=0;i<VW.ncols;i++){
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
* @author Florian Kobierska, Jan Magnusson, Rob Spence and Mathias Bavay
*/
void Interpol2D::CurvatureCorrection(DEMObject& dem, const Grid2DObject& ta, Grid2DObject& grid)
{
	if(!grid.isSameGeolocalization(dem)) {
		throw IOException("Requested grid does not match the geolocalization of the DEM", AT);
	}
	const double dem_max_curvature = dem.max_curvature, dem_range_curvature=(dem.max_curvature-dem.min_curvature);
	if(dem_range_curvature==0.) return;

	const double orig_mean = grid.grid2D.getMean();

	for (size_t j=0;j<grid.nrows;j++) {
		for (size_t i=0;i<grid.ncols;i++) {
			if(ta.grid2D(i,j)>273.15) continue; //modify the grid of precipitations only if air temperature is below or at freezing

			const double slope = dem.slope(i, j);
			const double curvature = dem.curvature(i, j);
			if(slope==IOUtils::nodata || curvature==IOUtils::nodata) continue;

			double& val = grid.grid2D(i, j);
			if(val!=IOUtils::nodata && dem_range_curvature!=0.) { //cf Huss
				val *= 0.5-(curvature-dem_max_curvature) / dem_range_curvature;
			}
		}
	}

	//HACK: correction for precipitation sum over the whole domain
	//this is a cheap/crappy way of compensating for the spatial redistribution of snow on the slopes
	const double new_mean = grid.grid2D.getMean();
	if(new_mean!=0.) grid.grid2D *= orig_mean/new_mean;

}

void Interpol2D::steepestDescentDisplacement(const DEMObject& dem, const Grid2DObject& grid, const size_t& ii, const size_t& jj, short &d_i_dest, short &d_j_dest)
{
	double max_slope = 0.;
	d_i_dest = 0, d_j_dest = 0;

	//loop around all adjacent cells to find the cell with the steepest downhill slope
	for(short d_i=-1; d_i<=1; d_i++) {
		for(short d_j=-1; d_j<=1; d_j++) {
			const double elev_pt1 = dem.grid2D(ii, jj);
			const double elev_pt2 = dem.grid2D(ii + d_i, jj + d_j);
			const double precip_1 = grid.grid2D(ii, jj);
			const double precip_2 = grid.grid2D(ii + d_i, jj + d_j);
			const double height_ratio = (elev_pt1+precip_1) / (elev_pt2+precip_2);
			const double new_slope = dem.slope(ii + d_i, jj + d_j);

			if ((new_slope>max_slope) && (height_ratio>1.)){
				max_slope = new_slope;
				d_i_dest = d_i;
				d_j_dest = d_j;
			}
		}
	}
}

double Interpol2D::depositAroundCell(const DEMObject& dem, const size_t& ii, const size_t& jj, const double& precip, Grid2DObject &grid)
{
	//else add precip to the cell and remove the same amount from the precip variable
	grid.grid2D(ii, jj) += precip;
	double distributed_precip = precip;

	for(short d_i=-1;d_i<=1;d_i++){
		for(short d_j=-1;d_j<=1;d_j++){
			const double elev_pt1 = dem.grid2D(ii, jj);
			const double elev_pt2 = dem.grid2D(ii + d_i, jj + d_j);
			const double precip_1 = grid.grid2D(ii, jj);
			const double precip_2 = grid.grid2D(ii + d_i, jj + d_j);
			const double height_ratio = (elev_pt1+precip_1) / (elev_pt2+precip_2);

			if ((d_i!=0)||(d_j!=0)){
				if (height_ratio>1.){
					grid.grid2D(ii + d_i, jj + d_j) += precip;
					distributed_precip += precip;
				}
			}
		}
	}

	return distributed_precip;
}

/**
 * @brief redistribute precip from steeper slopes to gentler slopes by following the steepest path from top to bottom
 * and gradually depositing precip during descent
 * @param dem array of elevations (dem). The slope must have been updated as it is required for the DEM analysis.
 * @param ta array of air temperatures used to determine if precipitation is rain or snow
 * @param grid 2D array of precipitation to fill
 * @author Rob Spence and Mathias Bavay
 */
void Interpol2D::SteepSlopeRedistribution(const DEMObject& dem, const Grid2DObject& ta, Grid2DObject& grid)
{
	for (size_t jj=1; jj<(grid.nrows-1); jj++) {
		for (size_t ii=1; ii<(grid.ncols-1); ii++) {
			if(grid.grid2D(ii,jj)==IOUtils::nodata) continue;
			if(ta.grid2D(ii, jj)>273.15) continue; //modify precipitation only for air temperatures at or below freezing

			const double slope = dem.slope(ii, jj);
			const double curvature = dem.curvature(ii, jj);
			if(slope==IOUtils::nodata || curvature==IOUtils::nodata) continue;
			if(slope<=40.) continue; //redistribution only above 40 degrees

			//remove all precip above 60 deg or linearly decrease it
			double precip = (slope>60.)? grid.grid2D(ii, jj) : grid.grid2D(ii, jj) * ((40.-slope)/-30.);
			grid.grid2D(ii, jj) -= precip; //we will redistribute the precipitation in a different way

			const double increment = precip / 50.; //break removed precip into smaller amounts to be redistributed
			double counter = 0.5;                  //counter will determine amount of precip deposited

			size_t ii_dest = ii, jj_dest = jj;
			while(precip>0.) {
				short d_i, d_j;
				steepestDescentDisplacement(dem, grid, ii_dest, jj_dest, d_i, d_j);
				//move to the destination cell
				ii_dest += d_i;
				jj_dest += d_j;

				//if (((ii_dest+d_i)<0) || ((jj_dest+d_j)<0) || ((ii_dest+d_i)>=grid.ncols)|| ((jj_dest+d_j)>=grid.nrows)) {
				if ((ii_dest==0) || (jj_dest==0) || (ii_dest==(grid.ncols-1))|| (jj_dest==(grid.nrows-1))){
					//we are getting out of the domain: deposit local contribution
					grid.grid2D(ii_dest, jj_dest) += counter*increment;
					break;
				}
				if(d_i==0 && d_j==0) {
					//local minimum, everything stays here...
					grid.grid2D(ii_dest, jj_dest) += precip;
					break;
				}

				precip -= depositAroundCell(dem, ii_dest, jj_dest, counter*increment, grid);
				counter += 0.25; //greater amount of precip is deposited as we move down the slope
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
