#include "Grid2DObject.h"
#include "IOUtils.h"

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
	latitude = IOUtils::nodata;
	longitude = IOUtils::nodata;
	cellsize = 0.0;
}

Grid2DObject::Grid2DObject(const unsigned int ncols_in, const unsigned int nrows_in,
				const double xllcorner_in, const double yllcorner_in,
				const double latitude_in, const double longitude_in,
				const double cellsize_in)
{
	set(ncols_in, nrows_in, xllcorner_in, yllcorner_in, latitude_in, longitude_in, cellsize_in);
}

Grid2DObject::Grid2DObject(const unsigned int ncols_in, const unsigned int nrows_in,
				const double xllcorner_in, const double yllcorner_in,
				const double latitude_in, const double longitude_in,
				const double cellsize_in, CArray2D<double>& grid2D_in)
{
	set(ncols_in, nrows_in, xllcorner_in, yllcorner_in, latitude_in, longitude_in, cellsize_in, grid2D_in);
}

void Grid2DObject::set(const unsigned int ncols_in, const unsigned int nrows_in,
			const double xllcorner_in, const double yllcorner_in,
			const double latitude_in, const double longitude_in,
			const double cellsize_in)
{
	ncols = ncols_in;
	nrows = nrows_in;
	grid2D.Create(ncols, nrows, IOUtils::nodata);
	cellsize = cellsize_in;
	
	xllcorner = xllcorner_in;
	yllcorner = yllcorner_in;
	latitude = latitude_in;
	longitude = longitude_in;
	
	//calculate/check coordinates if necessary
	if(latitude==IOUtils::nodata || longitude==IOUtils::nodata) {
		if(xllcorner==IOUtils::nodata || yllcorner==IOUtils::nodata) {
			throw InvalidArgumentException("missing positional parameters (xll,yll) or (lat,long) for Grid2DObject", AT);
		}
		IOUtils::CH1903_to_WGS84(xllcorner, yllcorner, latitude, longitude); //HACK: replace by local_to_WGS84
	} else {
		if(xllcorner==IOUtils::nodata || yllcorner==IOUtils::nodata) {
			IOUtils::WGS84_to_CH1903(latitude, longitude, xllcorner, yllcorner);  //HACK: replace by WGS84_to_local
		} else {
			double tmp_lat, tmp_lon;
			IOUtils::CH1903_to_WGS84(xllcorner, yllcorner, tmp_lat, tmp_lon); //HACK: replace by WGS84_to_local
			if(!IOUtils::checkEpsilonEquality(latitude, tmp_lat, 1.e-4) || !IOUtils::checkEpsilonEquality(longitude, tmp_lon, 1.e-4)) {
				throw InvalidArgumentException("Latitude/longitude and xllcorner/yllcorner don't match for Grid2DObject", AT);
			}
		}
	}

}

void Grid2DObject::set(const unsigned int ncols_in, const unsigned int nrows_in,
			const double xllcorner_in, const double yllcorner_in,
			const double latitude_in, const double longitude_in,
			const double cellsize_in, CArray2D<double>& grid2D_in)
{
	set(ncols_in, nrows_in, xllcorner_in, yllcorner_in, latitude_in, longitude_in, cellsize_in);

	//Test for equality in size: Only compatible CArray2D<double> grids are permitted
	int nx, ny;
	grid2D_in.GetSize(nx, ny);
	if (((int)ncols != nx) || ((int)nrows != ny)) {
		throw IOException("Mismatch in size of CArray2D<double> parameter grid2D_in and size of Grid2DObject", AT);
	}

	//Copy by value, after destroying the old grid
	grid2D.Create(ncols, nrows);
	for(unsigned int i=0; i<ncols; i++) {
		for(unsigned int j=0; j<nrows; j++) {
			grid2D(i,j) = grid2D_in(i,j);
		}
	}
}

Grid2DObject *Grid2DObject::sub(const unsigned int start_col, const unsigned int nb_cols)
{
	if(nb_cols==0) {
		throw InvalidArgumentException("requesting a subset of 0 columns for Grid2DObject", AT);
	}

	//allocate memory and get a pointer to it for the sub_grid
	Grid2DObject* sub_grid = new Grid2DObject(nb_cols, nrows, (xllcorner+start_col*cellsize), yllcorner,
					IOUtils::nodata, IOUtils::nodata, cellsize);

	//filling the grid
	for (unsigned int i=0; i<nb_cols; i++) {
		for (unsigned int j=0; j < nrows; j++){
			sub_grid->grid2D(i,j) = grid2D(i+start_col,j);
		}
	}

	//returning the pointer
	return sub_grid;

}

#ifdef _POPC_
#include "marshal_meteoio.h"
void Grid2DObject::Serialize(POPBuffer &buf, bool pack)
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
		marshal_TYPE_DOUBLE2D(buf, grid2D, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

