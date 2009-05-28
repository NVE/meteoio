#include "Grid2DObject.h"

/*
 * Default constructor. 
 * grid2D attribute is initialized by CArray2D default constructor.
 */
Grid2DObject::Grid2DObject() : grid2D()
{
	ncols = 0;
	nrows = 0;
	xllcorner = 0.0;
	yllcorner = 0.0;
	cellsize = 0.0;
	nodata = -9999.0;
}

Grid2DObject::Grid2DObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
			   const double& xllcorner_in, const double& yllcorner_in,
			   const double& cellsize_in, const double& nodata_in)
{
	set(ncols_in, nrows_in, xllcorner_in, yllcorner_in, cellsize_in, nodata_in);
}

Grid2DObject::Grid2DObject(const unsigned int& ncols_in, const unsigned int& nrows_in,
				const double& xllcorner_in, const double& yllcorner_in,
				const double& cellsize_in, const double& nodata_in, CArray2D<double>& grid2D_in)
{
	set(ncols_in, nrows_in, xllcorner_in, yllcorner_in, cellsize_in, nodata_in, grid2D_in);
}

void Grid2DObject::set(const unsigned int& ncols_in, const unsigned int& nrows_in,
		       const double& xllcorner_in, const double& yllcorner_in,
		       const double& cellsize_in, const double& nodata_in)
{
	ncols = ncols_in;
	nrows = nrows_in;
	xllcorner = xllcorner_in;
	yllcorner = yllcorner_in;
	cellsize = cellsize_in;
	nodata = nodata_in;

	grid2D.Create(ncols, nrows); //This function destroys the grid2D and rebuilds it with size ncols x nrows
}

void Grid2DObject::set(const unsigned int& ncols_in, const unsigned int& nrows_in,
		       const double& xllcorner_in, const double& yllcorner_in,
		       const double& cellsize_in, const double& nodata_in, CArray2D<double>& grid2D_in)
{
	ncols = ncols_in;
	nrows = nrows_in;
	xllcorner = xllcorner_in;
	yllcorner = yllcorner_in;
	cellsize = cellsize_in;
	nodata = nodata_in;

	//Test for equality in size: Only compatible CArray2D<double> grids are permitted
	int nx, ny;
	grid2D_in.GetSize(nx, ny);
	if (((int)ncols != nx) || ((int)nrows != ny)) {
		THROW IOException("Mismatch in size of CArray2D<double> parameter grid2D_in and size of Grid2DObject", AT);
	}

	grid2D = grid2D_in; //Copy by value, after destroying the old grid
}

