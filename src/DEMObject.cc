#include <math.h>
#include "DEMObject.h"


/**
* @brief Constructor that sets variables from a Grid2DObject
* @param dem_in (Grid2DObject&) grid contained in a Grid2DObject
*/
DEMObject::DEMObject(const Grid2DObject& dem_in)
{
	set(dem_in.ncols, dem_in.nrows, dem_in.xllcorner, dem_in.yllcorner, dem_in.latitude, dem_in.longitude, dem_in.cellsize, dem_in.grid2D);
}

/**
* @brief Extract a subset of the dem, given a starting column (first index being 0) and a number of columns
* @param start_col (unsigned int) starting column index
* @param nb_cols (unsigned int) number of columns
* @return new DEMObject
*/
DEMObject* DEMObject::sub(const unsigned int& start_col, const unsigned int& nb_cols)
{
	if(nb_cols==0) {
		throw InvalidArgumentException("requesting a subset of 0 columns for DEMObject", AT);
	}

	//allocate memory and get a pointer to it for the sub_dem
	DEMObject* sub_dem = new DEMObject(nb_cols, nrows, (xllcorner+start_col*cellsize), yllcorner, IOUtils::nodata, IOUtils::nodata, cellsize);

	//DEMObject* sub_dem = new DEMObject();
	//sub_dem->set(nb_cols, nrows, (xllcorner+start_col*cellsize), yllcorner, IOUtils::nodata, IOUtils::nodata, cellsize);

	//filling the grids
	for (unsigned int i=0; i<nb_cols; i++) {
		for (unsigned int j=0; j < nrows; j++){
			sub_dem->grid2D(i,j) = grid2D(i+start_col,j);
			sub_dem->slope(i,j) = slope(i+start_col,j);
			sub_dem->azi(i,j) = azi(i+start_col,j);
			sub_dem->curvature(i,j) = curvature(i+start_col,j);
			sub_dem->Nx(i,j) = Nx(i+start_col,j);
			sub_dem->Ny(i,j) = Ny(i+start_col,j);
			sub_dem->Nz(i,j) = Nz(i+start_col,j);
		}
	}

	//returning the pointer
	return sub_dem;
}

/**
* @brief Force the computation of the local slope, azimuth, normal vector and curvature.
* It has to be called manually since it can require some time to compute. Without this call, 
* the above mentionned parameters are NOT up to date.
*/
void DEMObject::Update() {
//This method recomputes the attributes that are not read as parameters 
//(such as slope, azimuth, normal vector)

	// Creating tables
	slope.Create(ncols, nrows);
	azi.Create(ncols, nrows);
	curvature.Create(ncols, nrows);
	Nx.Create(ncols, nrows);
	Ny.Create(ncols, nrows);
	Nz.Create(ncols, nrows);

	// Calculate the slope and the slope azimuth
	CalculateAziSlopeCurve(0);
}

void DEMObject::CalculateAziSlopeCurve(const int& Hick) {
//This computes the slope and the aspect at a given cell as well as the x and y components of the normal vector
//****        DEFINITIONS           ****//
// coordinate:                          //
//     x (m) asc from west to east      //
//     y (m) asc from south to north    //
//     z (m) elevation                  //
//                                      //
// index:                               //
//     [0][0]        lower left corner  //
//     [0][ny-1]     upper left corner  //
//     [nx-1][0]     lower right corner //
//     [nx-1][ny-1]  upper right corner //
//**************************************//
	double dx1, dx2, dx3, dy1, dy2, dy3;	// dxi = dz_i / dx_i and dyi = dz_i / dy_i
	double dx_sum, dy_sum;		// sum of dxi resp. dyi
		
	for ( unsigned int i = 0; i < ncols; i++ ) {
		for ( unsigned int j = 0; j < nrows; j++ ) {
			if( grid2D(i,j) == IOUtils::nodata ) {
				Nx(i,j) = Ny(i,j) = Nz(i,j) = slope(i,j) = IOUtils::nodata;
				azi(i,j) = IOUtils::nodata;
				curvature(i,j) = IOUtils::nodata;
			} else {
				CalculateSurfaceDeltas(i, j, &dx1, &dx2, &dx3, &dy1, &dy2, &dy3);
				dx_sum = (dx1 + dx2 + dx3) * 0.25;
				dy_sum = (dy1 + dy2 + dy3) * 0.25;
		
				//calculate slope and normal vector's components using one of two algorithms 
				//(knowing that Hick is the worst performing)
				if ( Hick ) {
					CalculateHickNormal(i, j, dx_sum, dy_sum);
				} else {
					CalculateCorripioNormal(i, j, dx_sum, dy_sum);
				}
		
				azi(i,j) = CalculateAzi(Nx(i,j), Ny(i,j), Nz(i,j), slope(i,j));
				curvature(i,j) = getCurvature(i,j);
			}
		}
	}

} // end of CalculateAziSlope

double DEMObject::CalculateAzi(const double& Nx, const double& Ny, const double& Nz, const double& slope) {
//Calculates the aspect at a given point knowing its normal vector and slope
//(direction of the normal pointing out of the surface, clockwise from north)
//This azimuth calculation is similar to Hodgson (1998)
	const double N_norm = sqrt( Nx*Nx + Ny*Ny + Nz*Nz );

	if ( slope > 0. ) { //there is some slope
		if ( Nx == 0. ) { //no E-W slope, so it is purely N-S
			if ( Ny < 0. ) {
				return (M_PI);	  // south facing
			} else {
				return (0.);	  // north facing
			}
		} else { //there is a E-W slope
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


void DEMObject::CalculateHickNormal(const unsigned int& i, const unsigned int& j, const double& dx_sum, const double& dy_sum) {
//This calculates the surface normal vector using the steepest slope method:
//the steepest slope found in the eight cells surrounding (i,j) is given to be the slope in (i,j)
//Beware, sudden steps could happen
//The supposedly better maximum downhill slope method (Dunn and Hickey, 1998)
	double smax=0.;			//maximum slope of all neighboring slopes
	double dmax;			//largest flat distance of a grid cell
	const double sqrt2=sqrt(2.);	//the weight of the 4 corner cells is increased by sqrt(2)
	int k,l;

	for ( k = -1; k < 2; k++ ) {
		for ( l= -1; l < 2; l++ ) {
			if ( k == 0 || l == 0 ) dmax = cellsize;
			else dmax = cellsize * sqrt2;

			//this could be done using the dx1, dx2, dx3... 
			smax = MAX( smax,fabs( grid2D(i,j) - 
					grid2D(MIN( ncols-1,MAX( 0,i+k ) ), MIN( nrows-1,MAX( 0,j+l ) )) ) / dmax );
		}
	}
	
	slope(i,j) = atan( smax );

	//Nx and Ny: x and y components of the normal pointing OUT of the surface
	if ( smax > 0. ) { //ie: there is some slope
		Nx(i,j) = -1.0 * dx_sum / (2. * cellsize);	//Nx=-dz/dx
		Ny(i,j) = -1.0 * dy_sum / (2. * cellsize);	//Ny=-dz/dy
		Nz(i,j) = 1.;				//Nz=1 (normalized by definition of Nx and Ny)
	} else { //ie: there is no slope
		Nx(i,j) = 0.;
		Ny(i,j) = 0.;
		Nz(i,j) = 1.;
	}
}

void DEMObject::CalculateCorripioNormal(const unsigned int& i, const unsigned int& j, const double& dx_sum, const double& dy_sum) {
//This calculates the surface normal vector using the two triangle method given in Corripio (2002)

	if ( i > 0 && i < ncols - 1 && j > 0 && j < nrows - 1 ) { // for not border grid cells
		// See Corripio (2002), knowing that here we normalize the result (divided by Nz=cellsize*cellsize)
		Nx(i,j) = (0.5 * ((grid2D(i,j) - grid2D(i+1,j) + grid2D(i,j+1) - grid2D(i+1,j+1) ))) / cellsize;
		Ny(i,j) = (0.5 * ((grid2D(i,j) + grid2D(i+1,j) - grid2D(i,j+1) - grid2D(i+1,j+1) ))) / cellsize;
		Nz(i,j) = 1.;
	} else { // for border grid cells
		// eight-neighbor algorithm of Horn (1981) 
		// normal vector n (not normalized!) from the slopes in x and y direction: "n=(1,0,mx)x(0,1,my)"
		// the conventional friendly neighbor method
		Nx(i,j) = -1.0 * dx_sum / (2. * cellsize);
		Ny(i,j) = -1.0 * dy_sum / (2. * cellsize);
		Nz(i,j) = 1.;
	}

	//There is no difference between slope = acos(n_z/|n|) and slope = atan(sqrt(sx*sx+sy*sy))
	//slope = acos( (Nz / sqrt( Nx*Nx + Ny*Ny + Nz*Nz )) );
	slope(i,j) = atan( sqrt( (Nx(i,j))*(Nx(i,j)) + (Ny(i,j))*(Ny(i,j)) ) );
		
	if ( slope(i,j) == 0. ) {
		Nx(i,j) = 0.;
		Ny(i,j) = 0.;
		Nz(i,j) = 1.;
	}
}

void DEMObject::CalculateSurfaceDeltas(const unsigned int& i, const unsigned int& j, double *dx1, double *dx2, double *dx3, double *dy1, double *dy2, double *dy3) {
//Compute the deltaZ for a given cell (i,j) accross its eight surrounding cells
// principal axis twice to emphasize height difference in that direction
// gradient estimation method of Horn (1981)
	if ( i > 0 && i < ncols - 1 && j > 0 && j < nrows - 1 ) { //away from the borders
		// (*dx1+*dx2+*dx3) is an approximation of dz/*dx 
		*dx1 =       grid2D(i+1,j-1) - grid2D(i-1,j-1);
		*dx2 = 2. * (grid2D(i+1,j  ) - grid2D(i-1,j  ));
		*dx3 =       grid2D(i+1,j+1) - grid2D(i-1,j+1);
		
		// (*dy1+*dy2+*dy3) is an approximation of dz/*dy
		*dy1 =       grid2D(i-1,j+1) - grid2D(i-1,j-1);
		*dy2 = 2. * (grid2D(i  ,j+1) - grid2D(i  ,j-1));
		*dy3 =       grid2D(i+1,j+1) - grid2D(i+1,j-1);
	}

	//borders
	else if ( i == 0 && j > 0 && j < nrows - 1 ) {
		*dx1 =       grid2D(i+1,j-1) - grid2D(i  ,j-1);
		*dx2 = 2. * (grid2D(i+1,j  ) - grid2D(i  ,j  )); 
		*dx3 =       grid2D(i+1,j+1) - grid2D(i  ,j+1);	 
		*dy1 =       grid2D(i  ,j+1) - grid2D(i  ,j-1);
		*dy2 = 2. *  (*dy1);
		*dy3 =       grid2D(i+1,j+1) - grid2D(i+1,j-1);
	}
	else if ( j == 0 && i > 0 && i < ncols - 1 ) {
		*dx1 =       grid2D(i+1,j  ) - grid2D(i-1,j  );
		*dx2 = 2. *  (*dx1);
		*dx3 =       grid2D(i+1,j+1) - grid2D(i-1,j+1);
		*dy1 =       grid2D(i-1,j+1) - grid2D(i-1,j  );
		*dy2 = 2. * (grid2D(i  ,j+1) - grid2D(i  ,j  ));
		*dy3 =       grid2D(i+1,j+1) - grid2D(i+1,j  );
	}
	else if ( i == ncols - 1 && j > 0 && j < nrows - 1 ) {
		*dx1 =       grid2D(i  ,j-1) - grid2D(i-1,j-1);
		*dx2 = 2. * (grid2D(i  ,j  ) - grid2D(i-1,j  )); 
		*dx3 =       grid2D(i  ,j+1) - grid2D(i-1,j+1);
		*dy1 =       grid2D(i-1,j+1) - grid2D(i-1,j-1);
		*dy2 = 2. * (grid2D(i  ,j+1) - grid2D(i  ,j-1));
		*dy3 =       (*dy2) * 0.5;
	}
	else if ( j == nrows - 1 && i > 0 && i < ncols - 1 ) {
		*dx1 =       grid2D(i+1,j-1) - grid2D(i-1,j-1);
		*dx2 = 2. * (grid2D(i+1,j  ) - grid2D(i-1,j  ));
		*dx3 =       (*dx2) * 0.5;
		*dy1 =       grid2D(i-1,j  ) - grid2D(i-1,j-1);
		*dy2 = 2. * (grid2D(i  ,j  ) - grid2D(i  ,j-1));
		*dy3 =       grid2D(i+1,j  ) - grid2D(i+1,j-1);
	}

	// four edges
	else if ( i == 0 && j == 0 ) {
		*dx2 = 2. * (grid2D(i+1,j  ) - grid2D(i  ,j  ));
		*dx1 =       (*dx2) * 0.5;
		*dx3 =       grid2D(i+1,j+1) - grid2D(i  ,j+1);
		*dy1 =       grid2D(i  ,j+1) - grid2D(i  ,j  );
		*dy2 = 2. *  (*dy1);
		*dy3 =       grid2D(i+1,j+1) - grid2D(i+1,j  );
	}
	else if ( i == ncols - 1 && j == nrows - 1 ) {
		*dx1 =       grid2D(i  ,j-1) - grid2D(i-1,j-1);
		*dx2 = 2. * (grid2D(i  ,j  ) - grid2D(i-1,j  )); 
		*dx3 =       (*dx2) * 0.5;
		*dy1 =       grid2D(i-1,j  ) - grid2D(i-1,j-1);
		*dy3 =       grid2D(i  ,j  ) - grid2D(i  ,j-1);
		*dy2 = 2. *  (*dy3);
	}
	else if ( i == 0 && j == nrows - 1 ) {
		*dx1 =       grid2D(i+1,j-1) - grid2D(i  ,j-1);
		*dx2 = 2. * (grid2D(i+1,j  ) - grid2D(i  ,j  ));
		*dx3 =       (*dx2) * 0.5;
		*dy1 =       grid2D(i  ,j  ) - grid2D(i  ,j-1);
		*dy2 = 2. * (*dy1);
		*dy3 =       grid2D(i+1,j  ) - grid2D(i+1,j-1);
	}
	
	//if ( j == 0 && i == ncols - 1 ) 
	else {
		*dx1 =       grid2D(i  ,j  ) - grid2D(i-1,j  ); 
		*dx2 = 2. * (*dx1); 
		*dx3 =       grid2D(i  ,j+1) - grid2D(i-1,j+1);
		*dy1 =       grid2D(i-1,j+1) - grid2D(i-1,j  );
		*dy3 =       grid2D(i  ,j+1) - grid2D(i  ,j  );
		*dy2 = 2. *  (*dy3);
	}
}
	
double DEMObject::OppositeDir(const double& z1, const double& z2, const double& z) {
//This method chooses the right topographic height to make the calcul
//and returns the sum of the opposites
	if ( z1 != 0. && z2 != 0.) {
		return z1 + z2;
	}
	// If the two values equals 0, is taken as if the landscape was flat
	else if ( z1 == 0. && z2 == 0. ) {
		return z * 2.;
	}
	else if ( z1 == 0. ) {
		return z2 + z;
	}
	// if ( z2 == 0 )
	else {
		return z1 + z;
	}
}

double DEMObject::getCurvature(const unsigned int& i, const unsigned int& j) {
//This methode computes and returns the curvature of a specific cell

	// Declaration
	double Zn, Zs, Ze, Zw;
	double Zne, Znw, Zse, Zsw;	
	double Zwe, Zsn, Zswne, Znwse;

	//center of the grid
	if ( i > 0 && i < ncols - 1 && j > 0 && j < nrows - 1 ){
		Zs  = grid2D(i  , j+1);
		Ze  = grid2D(i+1, j  );
		Zse = grid2D(i+1, j+1);
		Znw = grid2D(i-1, j-1);
		Zsw = grid2D(i-1, j+1);
		Zw  = grid2D(i-1, j  );
		Zne = grid2D(i+1, j-1);
		Zn  = grid2D(i  , j-1);
	}
	
	// Four edges 
	else if ( i == 0 && j == 0 ){
		Zs  = grid2D(i  , j+1);		
		Ze  = grid2D(i+1, j);
		Zse = grid2D(i+1, j+1);
		Znw = 0.;
		Zsw = 0.;
		Zw  = 0.;
		Zne = 0.;		
		Zn  = 0.;
	} 
	else if ( i == ncols-1 && j == 0 ){
		Zs  = grid2D(i  , j+1);
		Ze  = 0.;
		Zse = 0.;
		Znw = 0.;
		Zsw = grid2D(i-1, j+1);
		Zw  = grid2D(i-1, j);
		Zne = 0.;
		Zn  = 0.;
	} 
	else if ( i == 0 && j == nrows-1 ){
		Zs  = 0.;
		Ze  = grid2D(i+1, j  );
		Zse = 0.;
		Znw = 0.;
		Zsw = 0.;
		Zw  = 0.;
		Zne = grid2D(i+1, j-1);
		Zn  = grid2D(i  , j-1);
	} 
	else if ( i == ncols-1 && j == nrows-1 ){
		Zs  = 0.;
		Ze  = 0.;
		Zse = 0.;
		Znw = grid2D(i-1, j-1);
		Zsw = 0.;
		Zw  = grid2D(i-1, j  );
		Zne = 0.;
		Zn  = grid2D(i  , j-1);
	} 

	// Borders
	else if ( i > 0 && i < ncols-1 && j == 0 ) {
		Zs  = grid2D(i  , j+1);
		Ze  = grid2D(i+1, j  );
		Zse = grid2D(i+1, j+1);
		Znw = 0.;
		Zsw = grid2D(i-1, j+1);
		Zw  = grid2D(i-1, j  );
		Zne = 0.;
		Zn  = 0.;		
	} 
	else if ( i == 0 && j > 0 && j < nrows-1 ) {
		Zs  = grid2D(i  , j+1);
		Ze  = grid2D(i+1, j  );
		Zse = grid2D(i+1, j+1);
		Znw = 0.;
		Zsw = 0.;
		Zw  = 0.;
		Zne = grid2D(i+1, j-1);
		Zn  = grid2D(i  , j-1);
	} 
	else if ( i > 0 && i < ncols-1 && j == nrows-1 ){
		Zs  = 0.;
		Ze  = grid2D(i+1, j  );
		Zse = 0.;
		Znw = grid2D(i-1, j-1);
		Zsw = 0.;
		Zw  = grid2D(i-1, j  );
		Zne = grid2D(i+1, j-1);
		Zn  = grid2D(i  , j-1);
	
	}
	// if(i == ncols-1 && j > 0 && j < nrows-1) 
	else {
		Zs  = grid2D(i  , j+1);
		Ze  = 0.;
		Zse = 0.;
		Znw = grid2D(i-1, j-1);
		Zsw = grid2D(i-1, j+1);
		Zw  = grid2D(i-1, j  );
		Zne = 0.;
		Zn  = grid2D(i  , j-1);
	}

	// Calculate the opposite direction W-E, N-S, SW-NE, SE-NW
	double z = grid2D(i, j);
	Zwe   = 0.5 * OppositeDir(Zw, Ze, z);
	Zsn   = 0.5 * OppositeDir(Zs, Zn, z);
	Zswne = 0.5 * OppositeDir(Zsw, Zne, z);
	Znwse = 0.5 * OppositeDir(Znw, Zse, z);
	
	// Return the curvature
	const double sqrt2 = sqrt(2.);
	return 0.25 * ( 0.5*(z-Zwe) + 0.5*(z-Zsn) + (0.5*(z-Zswne)/sqrt2) + (0.5*(z-Znwse)/sqrt2) );
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
		int x,y;
		grid2D.GetSize(x,y);
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
		grid2D.Destroy();//if(grid2D!=NULL)delete(grid2D);
		slope.Destroy();
		azi.Destroy();
		curvature.Destroy();
		Nx.Destroy();
		Ny.Destroy();
		Nz.Destroy();
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
