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
#include <math.h>
#include <limits>
#include "DEMObject.h"
#include "Grid2DObject.h"
#include "Array.h"

/**
* @file DEMObject.cc
* @brief implementation of the DEMBoject class
*/

/**
* @brief Default constructor.
* Initializes all variables to 0, except lat/long which are initialized to IOUtils::nodata
* @param _algorithm specify the default algorithm to use for slope computation (default=DFLT)
*/
DEMObject::DEMObject(const slope_type& _algorithm) : Grid2DObject(), slope(), azi(), curvature(), Nx(), Ny(), Nz()
{
	min_altitude = min_slope = min_curvature = std::numeric_limits<double>::max();
	max_altitude = max_slope = max_curvature = -std::numeric_limits<double>::max();
	slope_failures = curvature_failures = 0;
	setDefaultAlgorithm(_algorithm);
}

/**
* @brief Constructor that sets variables.
* @param _ncols (unsigned int) number of colums in the grid2D
* @param _nrows (unsigned int) number of rows in the grid2D
* @param _xllcorner (double) x-coordinate of lower left corner
* @param _yllcorner (double) y-coordinate of lower left corner
* @param _latitude (double) decimal latitude, can be IOUtils::nodata
* @param _longitude (double) decimal longitude, can be IOUtils::nodata
* @param _cellsize (double) value for cellsize in grid2D
* @param _algorithm specify the default algorithm to use for slope computation (default=DFLT)
*/
DEMObject::DEMObject(const unsigned int& _ncols, const unsigned int& _nrows,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude,
				const double& _cellsize, const slope_type& _algorithm)
	: Grid2DObject(_ncols, _nrows, _xllcorner, _yllcorner, _latitude, _longitude, _cellsize),
	  slope(), azi(), curvature(), Nx(), Ny(), Nz()
{
	min_altitude = min_slope = min_curvature = std::numeric_limits<double>::max();
	max_altitude = max_slope = max_curvature = -std::numeric_limits<double>::max();
	slope_failures = curvature_failures = 0;
	setDefaultAlgorithm(_algorithm);
}

/**
* @brief Constructor that sets variables.
* @param _ncols (unsigned int) number of colums in the grid2D
* @param _nrows (unsigned int) number of rows in the grid2D
* @param _xllcorner (double) x-coordinate of lower left corner
* @param _yllcorner (double) y-coordinate of lower left corner
* @param _latitude (double) decimal latitude, can be IOUtils::nodata
* @param _longitude (double) decimal longitude, can be IOUtils::nodata
* @param _cellsize (double) value for cellsize in grid2D
* @param _altitude (Array2D\<double\>&) grid2D of elevations
* @param _update (bool) also update slope/normals/curvatures and their min/max? (default=true)
* @param _algorithm specify the default algorithm to use for slope computation (default=DFLT)
*/
DEMObject::DEMObject(const unsigned int& _ncols, const unsigned int& _nrows,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude,
				const double& _cellsize, const Array2D<double>& _altitude,
				const bool& _update, const slope_type& _algorithm)
	: Grid2DObject(_ncols, _nrows, _xllcorner, _yllcorner, _latitude, _longitude, _cellsize, _altitude),
	  slope(), azi(), curvature(), Nx(), Ny(), Nz()
{
	slope_failures = curvature_failures = 0;
	setDefaultAlgorithm(_algorithm);
	if(_update==false) {
		updateAllMinMax();
	} else {
		update(_algorithm);
	}
}

/**
* @brief Constructor that sets variables from a Grid2DObject
* @param _dem (Grid2DObject&) grid contained in a Grid2DObject
* @param _update (bool) also update slope/normals/curvatures and their min/max? (default=true)
* @param _algorithm specify the default algorithm to use for slope computation (default=DFLT)
*/
DEMObject::DEMObject(const Grid2DObject& _dem, const bool& _update, const slope_type& _algorithm)
  : Grid2DObject(_dem.ncols, _dem.nrows, _dem.xllcorner, _dem.yllcorner, _dem.latitude, _dem.longitude, _dem.cellsize, _dem.grid2D), 
    slope(), azi(), curvature(), Nx(), Ny(), Nz()
{
	slope_failures = curvature_failures = 0;
	setDefaultAlgorithm(_algorithm);
	if(_update==false) {
		updateAllMinMax();
	} else {
		update(_algorithm);
	}
}

/**
* @brief Constructor that sets variables from a subset of another DEMObject, 
* given an origin (X,Y) (first index being 0) and a number of columns and rows
* @param _dem (DEMObject&) dem contained in a DEMDObject
* @param _nx (unsigned int&) X coordinate of the new origin
* @param _ny (unsigned int&) Y coordinate of the new origin
* @param _ncols (unsigned int&) number of columns for the subset dem
* @param _nrows (unsigned int&) number of rows for the subset dem
* @param _update (bool) also update slope/normals/curvatures and their min/max? (default=true)
* @param _algorithm specify the default algorithm to use for slope computation (default=DFLT)
*/
DEMObject::DEMObject(const DEMObject& _dem, const unsigned int& _nx, const unsigned int& _ny, 
				 const unsigned int& _ncols, const unsigned int& _nrows, const bool& _update, const slope_type& _algorithm)
	: Grid2DObject(_dem, _nx,_ny, _ncols,_nrows), slope(_dem.slope,_nx,_ny, _ncols,_nrows),
	  azi(_dem.azi,_nx,_ny, _ncols,_nrows), curvature(_dem.curvature,_nx,_ny, _ncols,_nrows),
	  Nx(_dem.Nx,_nx,_ny, _ncols,_nrows), Ny(_dem.Ny,_nx,_ny, _ncols,_nrows), Nz(_dem.Nz,_nx,_ny, _ncols,_nrows)
{
	if ((_ncols==0) || (_nrows==0)) {
		throw InvalidArgumentException("requesting a subset of 0 columns or rows for DEMObject", AT);
	}

	slope_failures = curvature_failures = 0;
	setDefaultAlgorithm(_algorithm);
	if(_update==false) {
		updateAllMinMax();
	} else {
		update(_algorithm);
	}
}

/**
* @brief Force the computation of the local slope, azimuth, normal vector and curvature.
* It has to be called manually since it can require some time to compute. Without this call, 
* the above mentionned parameters are NOT up to date.
* @param algorithm (slope_type&) algorithm to use for computing slope, azimuth and normals
*/
void DEMObject::update(const slope_type& algorithm) {
//This method recomputes the attributes that are not read as parameters 
//(such as slope, azimuth, normal vector)

	// Creating tables
	slope.resize(ncols, nrows);
	azi.resize(ncols, nrows);
	curvature.resize(ncols, nrows);
	Nx.resize(ncols, nrows);
	Ny.resize(ncols, nrows);
	Nz.resize(ncols, nrows);

	CalculateAziSlopeCurve(algorithm);
	updateAllMinMax();
}

/**
* @brief Force the computation of the local slope, azimuth, normal vector and curvature.
* It has to be called manually since it can require some time to compute. Without this call, 
* the above mentionned parameters are NOT up to date.
* @param algorithm (const string&) algorithm to use for computing slope, azimuth and normals
* it is either:
* - HICK that uses the maximum downhill slope method (Dunn and Hickey, 1998)
* - CORRIPIO that uses the surface normal vector using the two triangle method given in Corripio (2002)
* and the eight-neighbor algorithm of Horn (1981) for border cells.
* - CARDINAL uses CORRIPIO but discretizes the resulting azimuth to 8 cardinal directions and the slope is rounded to the nearest degree. Curvature and normals are left untouched.
* 
* The azimuth is always computed using the Hodgson (1998) algorithm.
*/
void DEMObject::update(const string& algorithm) {
//This method recomputes the attributes that are not read as parameters 
//(such as slope, azimuth, normal vector)
	slope_type type;

	if(algorithm.compare("HICK")==0) {
		type=HICK;
	} else if(algorithm.compare("CORRIPIO")==0) {
		type=CORR;
	} else if(algorithm.compare("CARDINAL")==0) {
		type=CARD;
	} else if(algorithm.compare("DEFAULT")==0) {
		type=DFLT;
	} else {
		throw InvalidArgumentException("Chosen slope algorithm " + algorithm + " not available", AT);
	}
	
	update(type);
}

/**
* @brief Sets the default slope calculation algorithm
* @param _algorithm specify the default algorithm to use for slope computation
*/
void DEMObject::setDefaultAlgorithm(const slope_type& _algorithm) {
//This method MUST be called by each constructor!
	if(_algorithm==DFLT) {
		dflt_algorithm = CORR;
	} else {
		dflt_algorithm = _algorithm;
	}
}

/**
* @brief Recomputes the min/max of altitude, slope and curvature
* It return +/- std::numeric_limits\<double\>::max() for a given parameter if its grid was empty/undefined
*/
void DEMObject::updateAllMinMax() {
//updates the min/max parameters of all 2D tables
	min_altitude = grid2D.getMin(IOUtils::PARSE_NODATA);
	max_altitude = grid2D.getMax(IOUtils::PARSE_NODATA);
	min_slope = slope.getMin(IOUtils::PARSE_NODATA);
	max_slope = slope.getMax(IOUtils::PARSE_NODATA);
	min_curvature = curvature.getMin(IOUtils::PARSE_NODATA);
	max_curvature = curvature.getMax(IOUtils::PARSE_NODATA);
}

/**
* @brief Prints the list of points that have an elevation different than nodata but no slope or curvature
* Such points can happen if they are surrounded by too many points whose elevation is nodata
* If no such points exist, it prints nothing.
*/
void DEMObject::printFailures() {
	bool header=true;

	for ( unsigned int i = 0; i < ncols; i++ ) {
		for ( unsigned int j = 0; j < nrows; j++ ) {
			if(((slope(i,j)==IOUtils::nodata) || (curvature(i,j)==IOUtils::nodata)) && (grid2D(i,j)!=IOUtils::nodata)) {
				if(header==true) {
					std::cout << "[i] DEM slope/curvature could not be computed at the following points \n";
					std::cout << "[i]\tGrid Point\tElevation\tSlope\tCurvature\n";
					header=false;
				}
				std::cout << "[i]\t(" << i << "," << j << ")" << "\t\t" << grid2D(i,j) << "\t\t" << slope(i,j) << "\t" << curvature(i,j) << "\n";
			}
		}
	}
	if(header==false) {
		std::cout << endl;
	}
}

/**
* @brief Clean up the DEM Object
* When computing the slope and curvature, it is possible to get points where the elevation is known
* but where no slope/azimuth/normals/curvature could be computed. This method sets the elevation to nodata for such points,
* so that latter use of the DEM would be simpler (simply test the elevation in order to know if the point can be used
* and it guarantees that all other informations are available).If the slope/azimuth/normals/curvature tables were manually updated, this method will NOT perform any work (it requires the count of slopes/curvature failures to be greater than zero)
*
* IMPORTANT: calling this method DOES change the table of elevations!
*/
void DEMObject::sanitize() {
	if(slope_failures>0 || curvature_failures>0) {
		for ( unsigned int i = 0; i < ncols; i++ ) {
			for ( unsigned int j = 0; j < nrows; j++ ) {
				if(((slope(i,j)==IOUtils::nodata) || (curvature(i,j)==IOUtils::nodata)) && (grid2D(i,j)!=IOUtils::nodata)) {
					grid2D(i,j) = IOUtils::nodata;
				}
			}
		}
	}
}

/**
* @brief Computes the horizontal distance between two points in a metric grid
* @param xcoord1 east coordinate of the first point
* @param ycoord1 north coordinate of the first point
* @param xcoord2 east coordinate of the second point
* @param ycoord2 north coordinate of the second point
* @return horizontal distance in meters
*
*/
double DEMObject::horizontalDistance(const double& xcoord1, const double& ycoord1, const double& xcoord2, const double& ycoord2)
{
	return sqrt( pow2(xcoord2-xcoord1) + pow2(ycoord2-ycoord1) );
}


/**
* @brief Returns the distance *following the terrain* between two coordinates
* @param xcoord1 east coordinate of the first point
* @param ycoord1 north coordinate of the first point
* @param xcoord2 east coordinate of the second point
* @param ycoord2 north coordinate of the second point
* @return distance following the terrain in meters
*
*/
double DEMObject::terrainDistance(const double& xcoord1, const double& ycoord1, const double& xcoord2, const double& ycoord2) {
	std::vector<POINT> vec_points;
	double distance=0.;
	unsigned int last_point=0; //point 0 is always the starting point

	getPointsBetween(xcoord1, ycoord1, xcoord2, ycoord2, vec_points);
	if(vec_points.size()<=1) {
		return 0.;
	}

	for(unsigned int ii=1; ii<vec_points.size(); ii++) {
		const int ix1=vec_points[last_point].ix;
		const int iy1=vec_points[last_point].iy;
		const int ix2=vec_points[ii].ix;
		const int iy2=vec_points[ii].iy;

		if(grid2D(ix2,iy2)!=IOUtils::nodata) {
			if(grid2D(ix1,iy1)!=IOUtils::nodata) {
				//distance += sqrt( pow2((ix2-ix1)*cellsize) + pow2((iy2-iy1)*cellsize) + pow2(grid2D(ix2,iy2)-grid2D(ix1,iy1)) );
				const double z1=grid2D(ix1,iy1);
				const double z2=grid2D(ix2,iy2);
				const double tmpx=pow2((double)(ix2-ix1)*cellsize);
				const double tmpy=pow2((double)(iy2-iy1)*cellsize);
				const double tmpz=pow2(z2-z1);
				distance += sqrt( tmpx + tmpy + tmpz );
			}
			last_point = ii;
		}
	}

	return distance;
}

/**
* @brief Returns a list of grid points that are on the straight line between two coordinates
* @param xcoord1 east coordinate of the first point
* @param ycoord1 north coordinate of the first point
* @param xcoord2 east coordinate of the second point
* @param ycoord2 north coordinate of the second point
* @param vec_points vector of points that are in between
*
*/
void DEMObject::getPointsBetween(double xcoord1, double ycoord1, double xcoord2, double ycoord2, vector<POINT>& vec_points) {

	if(xcoord1>xcoord2) {
		//we want xcoord1<xcoord2, so we swap the two points
		const double tmpx = xcoord1, tmpy= ycoord1;
		xcoord1 = xcoord2; ycoord1 = ycoord2;
		xcoord2 = tmpx, ycoord2 = tmpy;
	}

	//extension of the line segment (pts1, pts2) along the X axis
	const int ix1 = (int)floor( (xcoord1-xllcorner)/cellsize );
	const int iy1 = (int)floor( (ycoord1-yllcorner)/cellsize );
	const int ix2 = (int)floor( (xcoord2-xllcorner)/cellsize );
	const int iy2 = (int)floor( (ycoord2-yllcorner)/cellsize );

	if(ix1==ix2) {
		//special case of vertical alignement
		for(int iy=MIN(iy1,iy2); iy<=MAX(iy1,iy2); iy++) {
			POINT pts;
			pts.ix = ix1;
			pts.iy = iy;
			vec_points.push_back(pts);
		}
	} else {
		//normal case
		//equation of the line between the two points
		const double a = ((double)(iy2-iy1)) / ((double)(ix2-ix1));
		const double b = (double)iy1 - a * (double)ix1;

		for(int ix=ix1; ix<=ix2; ix++) {
			//extension of the line segment (ix, ix+1) along the Y axis
			int y1 = (int)floor( a*(double)ix+b );
			//const int y2 = MIN( (int)floor( a*((double)ix+1)+b ) , iy2);
			int y2 = (int)floor( a*((double)ix+1)+b );
			if(ix==ix2 && y1==iy2) {
				//we don't want to overshoot when reaching the target cell
				y2 = y1;
			}

			if(y1>y2) {
				//we want y1<y2, so we swap the two coordinates
				const double ytemp=y1;
				y1=y2; y2=ytemp;
			}

			for(int iy=y1; iy<=y2; iy++) {
				POINT pts;
				pts.ix = ix;
				pts.iy = iy;
				//make sure we only return points within the dem
				if(ix>0 && ix<(signed)ncols && iy>0 && iy<(signed)nrows) {
					vec_points.push_back(pts);
				}
			}
		}
	}
}

void DEMObject::CalculateAziSlopeCurve(slope_type algorithm) {
//This computes the slope and the aspect at a given cell as well as the x and y components of the normal vector
	double A[4][4]; //table to store neigbouring heights: 3x3 matrix but we want to start at [1][1]

	if(algorithm==DFLT) {
		algorithm = dflt_algorithm;
	}

	slope_failures = curvature_failures = 0;
	if(algorithm==HICK) {
		for ( unsigned int i = 0; i < ncols; i++ ) {
			for ( unsigned int j = 0; j < nrows; j++ ) {
				if( grid2D(i,j) == IOUtils::nodata ) {
					Nx(i,j) = Ny(i,j) = Nz(i,j) = slope(i,j) = IOUtils::nodata;
					azi(i,j) = IOUtils::nodata;
					curvature(i,j) = IOUtils::nodata;
				} else {
					getNeighbours(i, j, A);
					CalculateHick(A, slope(i,j), Nx(i,j), Ny(i,j), Nz(i,j));
					azi(i,j) = CalculateAspect(Nx(i,j), Ny(i,j), Nz(i,j), slope(i,j));
					curvature(i,j) = getCurvature(A);
				}
			}
		}
	} else if(algorithm==CORR) {
		for ( unsigned int i = 0; i < ncols; i++ ) {
			for ( unsigned int j = 0; j < nrows; j++ ) {
				if( grid2D(i,j) == IOUtils::nodata ) {
					Nx(i,j) = Ny(i,j) = Nz(i,j) = slope(i,j) = IOUtils::nodata;
					azi(i,j) = IOUtils::nodata;
					curvature(i,j) = IOUtils::nodata;
				} else {
					getNeighbours(i, j, A);
					CalculateCorripio(A, slope(i,j), Nx(i,j), Ny(i,j), Nz(i,j));
					azi(i,j) = CalculateAspect(Nx(i,j), Ny(i,j), Nz(i,j), slope(i,j));
					curvature(i,j) = getCurvature(A);
				}
			}
		}
	} else if(algorithm==FLEM) {
		for ( unsigned int i = 0; i < ncols; i++ ) {
			for ( unsigned int j = 0; j < nrows; j++ ) {
				if( grid2D(i,j) == IOUtils::nodata ) {
					Nx(i,j) = Ny(i,j) = Nz(i,j) = slope(i,j) = IOUtils::nodata;
					azi(i,j) = IOUtils::nodata;
					curvature(i,j) = IOUtils::nodata;
				} else {
					getNeighbours(i, j, A);
					CalculateFleming(A, slope(i,j), Nx(i,j), Ny(i,j), Nz(i,j));
					azi(i,j) = CalculateAspect(Nx(i,j), Ny(i,j), Nz(i,j), slope(i,j));
					curvature(i,j) = getCurvature(A);
				}
			}
		}
	} else if(algorithm==CARD) {
		for ( unsigned int i = 0; i < ncols; i++ ) {
			for ( unsigned int j = 0; j < nrows; j++ ) {
				if( grid2D(i,j) == IOUtils::nodata ) {
					Nx(i,j) = Ny(i,j) = Nz(i,j) = slope(i,j) = IOUtils::nodata;
					azi(i,j) = IOUtils::nodata;
					curvature(i,j) = IOUtils::nodata;
				} else {
					getNeighbours(i, j, A);
					CalculateHick(A, slope(i,j), Nx(i,j), Ny(i,j), Nz(i,j));
					azi(i,j) = CalculateAspect(Nx(i,j), Ny(i,j), Nz(i,j), slope(i,j));
					//TODO: instead of PI for flat, get nodata and process it by an extra algorithm
					curvature(i,j) = getCurvature(A);
					if(azi(i,j)!=IOUtils::nodata)
						azi(i,j) = fmod(floor( (azi(i,j)+22.5)/45. )*45., 360.);
					if(slope(i,j)!=IOUtils::nodata)
						slope(i,j) = floor( slope(i,j)+0.5 );
				}
			}
		}
	} else {
		throw InvalidArgumentException("Chosen slope algorithm not available", AT);
	}

	//Inform the user is some points have unexpectidly not been computed
	//(ie: there was an altitude but some parameters could not be computed)
	if(slope_failures>0 || curvature_failures>0) {
		std::cout << "[W] DEMObject: " << slope_failures << " point(s) have an elevation but no slope, " << curvature_failures << " point(s) have an elevation but no curvature." << std::endl;
	}

} // end of CalculateAziSlope

double DEMObject::CalculateAspect(const double& Nx, const double& Ny, const double& Nz, const double& slope) {
//Calculates the aspect at a given point knowing its normal vector and slope
//(direction of the normal pointing out of the surface, clockwise from north)
//This azimuth calculation is similar to Hodgson (1998)

	if(Nx==IOUtils::nodata || Ny==IOUtils::nodata || Nz==IOUtils::nodata || slope==IOUtils::nodata) {
		return IOUtils::nodata;
	}

	if ( slope > 0. ) { //there is some slope
		if ( Nx == 0. ) { //no E-W slope, so it is purely N-S
			if ( Ny < 0. ) {
				return (M_PI);	  // south facing
			} else {
				return (0.);	  // north facing
			}
		} else { //there is a E-W slope
			const double N_norm = sqrt( Nx*Nx + Ny*Ny + Nz*Nz );
			if ( Nx > 0. ) {
				return (M_PI * 0.5 - atan( (Ny/N_norm) / (Nx/N_norm) ));
			} else {
				return (M_PI * 1.5 - atan( (Ny/N_norm) / (Nx/N_norm) ));
			}
		}
	} else { // if slope = 0
		return (M_PI);          // undefined or plain surface
	}
}


void DEMObject::CalculateHick(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz) {
//This calculates the surface normal vector using the steepest slope method:
//the steepest slope found in the eight cells surrounding (i,j) is given to be the slope in (i,j)
//Beware, sudden steps could happen
//The supposedly better maximum downhill slope method (Dunn and Hickey, 1998)
	const double smax = steepestGradient(A); //steepest local gradient

	if(smax==IOUtils::nodata) {
		slope = IOUtils::nodata;
		Nx = IOUtils::nodata;
		Ny = IOUtils::nodata;
		Nz = IOUtils::nodata;
		slope_failures++;
	} else {
		slope = atan( smax );

		//Nx and Ny: x and y components of the normal pointing OUT of the surface
		if ( smax > 0. ) { //ie: there is some slope
			double dx_sum, dy_sum;
			surfaceGradient(dx_sum, dy_sum, A);
			if(dx_sum==IOUtils::nodata || dy_sum==IOUtils::nodata) {
				Nx = IOUtils::nodata;
				Ny = IOUtils::nodata;
				Nz = IOUtils::nodata;
				slope_failures++;
			} else {
				Nx = -1.0 * dx_sum / (2. * cellsize);	//Nx=-dz/dx
				Ny = -1.0 * dy_sum / (2. * cellsize);	//Ny=-dz/dy
				Nz = 1.;				//Nz=1 (normalized by definition of Nx and Ny)
			}
		} else { //ie: there is no slope
			Nx = 0.;
			Ny = 0.;
			Nz = 1.;
		}
	}
}

void DEMObject::CalculateFleming(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz) {
//This calculates the surface normal vector using method by Fleming and Hoffer (1979)

	if(A[2][1]!=IOUtils::nodata && A[2][3]!=IOUtils::nodata && A[3][2]!=IOUtils::nodata && A[1][2]!=IOUtils::nodata) {
		Nx = 0.5 * (A[2][1] - A[2][3]) / cellsize;
		Ny = 0.5 * (A[3][2] - A[1][2]) / cellsize;
		Nz = 1.;
		slope = atan( sqrt(Nx*Nx+Ny*Ny) );
	} else {
		CalculateHick(A, slope, Nx, Ny, Nz);
	}
}

void DEMObject::CalculateCorripio(double A[4][4], double& slope, double& Nx, double& Ny, double& Nz) {
//This calculates the surface normal vector using the two triangle method given in Corripio (2002)
	if ( A[1][2]!=IOUtils::nodata && A[1][3]!=IOUtils::nodata && A[2][2]!=IOUtils::nodata && A[2][3]!=IOUtils::nodata) {
		// See Corripio (2002), knowing that here we normalize the result (divided by Nz=cellsize*cellsize)
		Nx = (0.5 * (A[2][2] - A[2][3] + A[1][2] - A[1][3]) ) / cellsize;
		Ny = (0.5 * (A[2][2] + A[2][3] - A[1][2] - A[1][3]) ) / cellsize;
		Nz = 1.;
		//There is no difference between slope = acos(n_z/|n|) and slope = atan(sqrt(sx*sx+sy*sy))
		//slope = acos( (Nz / sqrt( Nx*Nx + Ny*Ny + Nz*Nz )) );
		slope = atan( sqrt(Nx*Nx+Ny*Ny) );
	} else {
		// eight-neighbor algorithm of Horn (1981)
		CalculateHick(A, slope, Nx, Ny, Nz);
	}
}

double DEMObject::getCurvature(double A[4][4]) {
//This methode computes the curvature of a specific cell
	if(A[2][2]!=IOUtils::nodata) {
		const double Zwe   = avgHeight(A[2][1], A[2][2], A[2][3]);
		const double Zsn   = avgHeight(A[1][2], A[2][2], A[3][2]);
		const double Zswne = avgHeight(A[3][1], A[2][2], A[1][3]);
		const double Znwse = avgHeight(A[1][1], A[2][2], A[3][3]);

		const double sqrt2 = sqrt(2.);
		double sum=0.;
		unsigned int count=0;

		if(Zwe!=IOUtils::nodata) {
			sum += 0.5*(A[2][2]-Zwe);
			count++;
		}
		if(Zsn!=IOUtils::nodata) {
			sum += 0.5*(A[2][2]-Zsn);
			count++;
		}
		if(Zswne!=IOUtils::nodata) {
			sum += 0.5*(A[2][2]-Zswne)/sqrt2;
			count++;
		}
		if(Znwse!=IOUtils::nodata) {
			sum += 0.5*(A[2][2]-Znwse)/sqrt2;
			count++;
		}

		if(count != 0.) return 1./count * sum;
	}
	curvature_failures++;
	return IOUtils::nodata;
}

double DEMObject::steepestGradient(double A[4][4]) {
//best effort to calculate the local steepest gradient
	double smax=-1.;		//maximum slope of all neighboring slopes
	const double sqrt2=sqrt(2.);	//the weight of the 4 corner cells is increased by sqrt(2)

	if(A[2][2]!=IOUtils::nodata) {
		if(A[1][1]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[1][1])/(cellsize*sqrt2) );
		if(A[1][2]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[1][2])/(cellsize) );
		if(A[1][3]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[1][3])/(cellsize*sqrt2) );
		if(A[2][1]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[2][1])/(cellsize) );
		if(A[2][3]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[2][3])/(cellsize) );
		if(A[3][1]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[3][1])/(cellsize*sqrt2) );
		if(A[3][2]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[3][2])/(cellsize) );
		if(A[3][3]!=IOUtils::nodata)
			smax = MAX( smax, fabs(A[2][2] - A[3][3])/(cellsize*sqrt2) );
	}

	if(smax<0.)
		return IOUtils::nodata;
	return smax;
}

double DEMObject::lineGradient(const double& A1, const double& A2, const double& A3) {
//best effort to calculate the local gradient
	double delta = IOUtils::nodata;

	if(A3!=IOUtils::nodata && A1!=IOUtils::nodata) {
		delta = A3 - A1;
	} else {
		if(A2!=IOUtils::nodata) {
			if(A3!=IOUtils::nodata) 
				delta = (A3 - A2) * 2.;
			if(A1!=IOUtils::nodata) 
				delta = (A2 - A1) * 2.;
		}
	}

	return delta;
}

double DEMObject::fillMissingGradient(const double& delta1, const double& delta2) {
//If a gradient could not be computed, try to fill it with some neighboring value
	double dx = IOUtils::nodata;

	if(delta1!=IOUtils::nodata && delta2!=IOUtils::nodata) {
		dx = 0.5 * (delta1+delta2);
	} else {
		if(delta1!=IOUtils::nodata) dx = delta1;
		if(delta2!=IOUtils::nodata) dx = delta2;
	}

	return dx;
}

void DEMObject::surfaceGradient(double& dx_sum, double& dy_sum, double A[4][4]) {
//Compute the gradient for a given cell (i,j) accross its eight surrounding cells (Horn, 1981)
	double dx1, dx2, dx3;
	double dy1, dy2, dy3;

	//general calculation
	dx1 = lineGradient(A[3][1], A[3][2], A[3][3]);
	dx2 = lineGradient(A[2][1], A[2][2], A[2][3]);
	dx3 = lineGradient(A[1][1], A[1][2], A[1][3]);

	dy1 = lineGradient(A[3][1], A[2][1], A[1][1]);
	dy2 = lineGradient(A[3][2], A[2][2], A[1][2]);
	dy3 = lineGradient(A[3][3], A[2][3], A[1][3]);

	//now trying to fill whatever could not be filled...
	if(dx1==IOUtils::nodata) dx1 = fillMissingGradient(dx2, dx3);
	if(dx2==IOUtils::nodata) dx2 = fillMissingGradient(dx1, dx3);
	if(dx3==IOUtils::nodata) dx3 = fillMissingGradient(dx1, dx2);
	if(dy1==IOUtils::nodata) dy1 = fillMissingGradient(dy2, dy3);
	if(dy2==IOUtils::nodata) dy2 = fillMissingGradient(dy1, dy3);
	if(dy3==IOUtils::nodata) dy3 = fillMissingGradient(dy1, dy2);

	if(dx1!=IOUtils::nodata && dy1!=IOUtils::nodata) {
		// principal axis twice to emphasize height difference in that direction
		dx_sum = (dx1 + 2.*dx2 + dx3) * 0.25;
		dy_sum = (dy1 + 2.*dy2 + dy3) * 0.25;
	} else {
		//if dx1==nodata, this also means that dx2==nodata and dx3==nodata
		//(otherwise, dx1 would have received a copy of either dx2 or dx3)
		dx_sum = IOUtils::nodata;
		dy_sum = IOUtils::nodata;
	}
}

double DEMObject::avgHeight(const double& z1, const double &z2, const double& z3) {
//this safely computes the average height accross a vector

	if(z1!=IOUtils::nodata && z3!=IOUtils::nodata) {
		return 0.5*(z1+z3);
	}
	if(z1!=IOUtils::nodata && z2!=IOUtils::nodata) {
		return 0.5*(z1+z2);
	}
	if(z3!=IOUtils::nodata && z2!=IOUtils::nodata) {
		return 0.5*(z3+z2);
	}

	return IOUtils::nodata;
}

void DEMObject::getNeighbours(const unsigned int i, const unsigned int j, double A[4][4]) {
//this fills a 3x3 table containing the neighbouring values
		A[1][1] = safeGet(i-1, j+1);
		A[1][2] = safeGet(i, j+1);
		A[1][3] = safeGet(i+1, j+1);
		A[2][1] = safeGet(i-1, j);
		A[2][2] = safeGet(i, j);
		A[2][3] = safeGet(i+1, j);
		A[3][1] = safeGet(i-1, j-1);
		A[3][2] = safeGet(i, j-1);
		A[3][3] = safeGet(i+1, j-1);
}

double DEMObject::safeGet(const int i, const int j)
{//this function would allow reading the value of *any* point, 
//that is, even for coordinates outside of the grid (where it would return nodata)
//this is to make implementing the slope/curvature computation easier for edges, holes, etc

	if(i<0 || i>=(signed)ncols) {
		return IOUtils::nodata;
	}
	if(j<0 || j>=(signed)nrows) {
		return IOUtils::nodata;
	}
	
	return grid2D((unsigned)i, (unsigned)j);
}

#ifdef _POPC_
#include "marshal_meteoio.h"
void DEMObject::Serialize(POPBuffer &buf, bool pack)
{
	if (pack)
	{
		buf.Pack(&ncols,1);
		buf.Pack(&nrows,1);
		buf.Pack(&xllcorner,1);
		buf.Pack(&yllcorner,1);
		buf.Pack(&latitude,1);
		buf.Pack(&longitude,1);
		buf.Pack(&cellsize,1);
		buf.Pack(&min_altitude,1);
		buf.Pack(&max_altitude,1);
		buf.Pack(&min_slope,1);
		buf.Pack(&max_slope,1);
		buf.Pack(&min_curvature,1);
		buf.Pack(&max_curvature,1);
		buf.Pack(&slope_failures,1);
		buf.Pack(&curvature_failures,1);
		//unsigned int x,y;
		//grid2D.size(x,y);
		marshal_slope_type(buf, dflt_algorithm, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, grid2D, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, slope, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, azi, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, curvature, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, Nx, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, Ny, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, Nz, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&ncols,1);
		buf.UnPack(&nrows,1);
		buf.UnPack(&xllcorner,1);
		buf.UnPack(&yllcorner,1);
		buf.UnPack(&latitude,1);
		buf.UnPack(&longitude,1);
		buf.UnPack(&cellsize,1);
		buf.UnPack(&min_altitude,1);
		buf.UnPack(&max_altitude,1);
		buf.UnPack(&min_slope,1);
		buf.UnPack(&max_slope,1);
		buf.UnPack(&min_curvature,1);
		buf.UnPack(&max_curvature,1);
		buf.UnPack(&slope_failures,1);
		buf.UnPack(&curvature_failures,1);
		grid2D.clear();//if(grid2D!=NULL)delete(grid2D);
		slope.clear();
		azi.clear();
		curvature.clear();
		Nx.clear();
		Ny.clear();
		Nz.clear();
		marshal_slope_type(buf, dflt_algorithm, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, grid2D, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, slope, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, azi, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, curvature, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, Nx, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, Ny, 0, !FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, Nz, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif
