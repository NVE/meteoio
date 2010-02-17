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
#include "Grid2DObject.h"
#include "IOUtils.h"
#include "Coords.h"
#include <cmath>

/*
 * Default constructor.
 * grid2D attribute is initialized by Array2D default constructor.
 */
Grid2DObject::Grid2DObject() : grid2D()
{
	ncols = 0;
	nrows = 0;
	cellsize = 0.0;
}

Grid2DObject::Grid2DObject(const unsigned int& _ncols, const unsigned int& _nrows,
				const double& _cellsize, const Coords& _llcorner) : grid2D(_ncols, _nrows, IOUtils::nodata)
{
	//set metadata, grid2D already successfully created
	setValues(_ncols, _nrows, _cellsize, _llcorner);
}

Grid2DObject::Grid2DObject(const unsigned int& _ncols, const unsigned int& _nrows,
				const double& _cellsize, const Coords& _llcorner, const Array2D<double>& _grid2D) : grid2D()
{
	set(_ncols, _nrows, _cellsize, _llcorner, _grid2D);
}

Grid2DObject::Grid2DObject(const Grid2DObject& _grid2Dobj, const unsigned int& _nx, const unsigned int& _ny,
				const unsigned int& _ncols, const unsigned int& _nrows) 
	: grid2D(_grid2Dobj.grid2D, _nx,_ny, _ncols,_nrows)
{
	setValues(_ncols, _nrows, _grid2Dobj.cellsize);

	//we take the previous corner (so we use the same projection parameters)
	//and we shift it by the correct X and Y distance
	llcorner = _grid2Dobj.llcorner;
	if( (llcorner.getEasting()!=IOUtils::nodata) && (llcorner.getNorthing()!=IOUtils::nodata) ) {
		llcorner.setXY( llcorner.getEasting()+_nx*_grid2Dobj.cellsize,
				llcorner.getNorthing()+_ny*_grid2Dobj.cellsize);
	}
}

void Grid2DObject::grid_to_WGS84(const unsigned int& i, const unsigned int& j, Coords& point)
{
	const double easting = ((double)i) * cellsize; //The coordinate the ll corner of the cell
	const double northing = ((double)j) * cellsize;

	if(point.isSameProj(llcorner)==false) {
		point.copyProj(llcorner);
	}
	point.setXY(easting, northing);
}

int Grid2DObject::WGS84_to_grid(Coords point, unsigned int& i, unsigned int& j)
{
	if(point.isSameProj(llcorner)==false) {
		point.copyProj(llcorner);
	}
	double x = floor( (point.getEasting()-llcorner.getEasting()) / cellsize );
	double y = floor( (point.getNorthing()-llcorner.getNorthing()) / cellsize );

	int error_code=EXIT_SUCCESS;
	
	if(x<0.) {
		i=0;
		error_code=EXIT_FAILURE;
	} else if(x>(double)ncols) {
		i=ncols;
		error_code=EXIT_FAILURE;
	} else {
		i=(int)x;
	}
	
	if(y<0.) {
		j=0;
		error_code=EXIT_FAILURE;
	} else if(y>(double)nrows) {
		j=nrows;
		error_code=EXIT_FAILURE;
	} else {
		j=(int)y;
	}

	return error_code;
}

void Grid2DObject::set(const unsigned int& _ncols, const unsigned int& _nrows,
			const double& _cellsize, const Coords& _llcorner)
{
	setValues(_ncols, _nrows, _cellsize, _llcorner);
	grid2D.resize(ncols, nrows, IOUtils::nodata);
}

void Grid2DObject::set(const unsigned int& _ncols, const unsigned int& _nrows,
			const double& _cellsize, const Coords& _llcorner, const Array2D<double>& _grid2D)
{
	//Test for equality in size: Only compatible Array2D<double> grids are permitted
	unsigned int nx,ny;
	_grid2D.size(nx,ny);
	if ((_ncols != nx) || (_nrows != ny)) {
		throw IOException("Mismatch in size of Array2D<double> parameter _grid2D and size of Grid2DObject", AT);
	}

	setValues(_ncols, _nrows, _cellsize, _llcorner);

	//Copy by value, after destroying the old grid
	grid2D = _grid2D;
}

void Grid2DObject::setValues(const unsigned int& _ncols, const unsigned int& _nrows,
				const double& _cellsize)
{
	ncols = _ncols;
	nrows = _nrows;
	cellsize = _cellsize;
}

void Grid2DObject::setValues(const unsigned int& _ncols, const unsigned int& _nrows,
				const double& _cellsize, const Coords& _llcorner)
{
	setValues(_ncols, _nrows, _cellsize);
	llcorner = _llcorner;
}

bool Grid2DObject::isSameGeolocalization(const Grid2DObject& target)
{
	if( ncols==target.ncols && nrows==target.nrows &&
		llcorner==target.llcorner &&
		cellsize==target.cellsize) {
		return true;
	} else {
		return false;
	}
}

#ifdef _POPC_
#include "marshal_meteoio.h"
void Grid2DObject::Serialize(POPBuffer &buf, bool pack)
{
	if (pack)
	{
		buf.Pack(&ncols,1);
		buf.Pack(&nrows,1);
		buf.Pack(&cellsize,1);
		marshal_Coords(buf, llcorner, 0, FLAG_MARSHAL, NULL);
		marshal_TYPE_DOUBLE2D(buf, grid2D, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&ncols,1);
		buf.UnPack(&nrows,1);
		buf.UnPack(&cellsize,1);
		marshal_Coords(buf, llcorner, 0, !FLAG_MARSHAL, NULL);
		grid2D.clear();//if(grid2D!=NULL)delete(grid2D);
		marshal_TYPE_DOUBLE2D(buf, grid2D, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

