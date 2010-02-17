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
#include "Grid3DObject.h"

Grid3DObject::Grid3DObject() : grid3D() //using Array3D default constructor
{
	ncols = nrows = ndepth = 0;
	cellsize = 0.0;
}

Grid3DObject::Grid3DObject(const Grid3DObject& _grid3Dobj,
				const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
				const unsigned int& _nwidth, const unsigned int& _nheight, const unsigned int& _ndepth) 
	: grid3D(_grid3Dobj.grid3D, _nx,_ny,_nz, _nwidth,_nheight,_ndepth)
{
	setValues(_nwidth, _nheight, _ndepth, _grid3Dobj.cellsize);

	//we take the previous corner (so we use the same projection parameters)
	//and we shift it by the correct X and Y distance
	llcorner = _grid3Dobj.llcorner;
	if( (llcorner.getEasting()!=IOUtils::nodata) && (llcorner.getNorthing()!=IOUtils::nodata) ) {
		llcorner.setXY( llcorner.getEasting()+_nx*_grid3Dobj.cellsize,
				llcorner.getNorthing()+_ny*_grid3Dobj.cellsize);
	}
}

Grid3DObject::Grid3DObject(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _cellsize, const Coords& _llcorner) : grid3D(_ncols, _nrows, _ndepth, IOUtils::nodata)
{
	setValues(_ncols, _nrows, _ndepth, _cellsize, _llcorner);
}

Grid3DObject::Grid3DObject(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
					  const double& _cellsize, const Coords& _llcorner, const Array3D<double>& _grid3D) : grid3D()
{
	set(_ncols, _nrows, _ndepth, _cellsize, _llcorner, _grid3D);
}

void Grid3DObject::set(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _cellsize, const Coords& _llcorner)
{
	setValues(_ncols, _nrows, _ndepth, _cellsize, _llcorner);
	grid3D.resize(ncols, nrows, ndepth, IOUtils::nodata);	
}

void Grid3DObject::set(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _cellsize, const Coords& _llcorner, const Array3D<double>& _grid3D)
{
	//Test for equality in size: Only compatible Array3D<double> grids are permitted
	unsigned int nx, ny, nz;
	_grid3D.size(nx, ny, nz);
	if ((_ncols != nx) || (_nrows != ny) || (_ndepth != nz)) {
		throw IOException("Mismatch in size of Array3D<double> parameter grid3D and size of Grid3DObject", AT);
	}

	setValues(_ncols, _nrows, _ndepth, _cellsize, _llcorner);
	grid3D = _grid3D; //copy by value
}

void Grid3DObject::setValues(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _cellsize)
{
	ncols = _ncols;
	nrows = _nrows;
	ndepth = _ndepth;
	cellsize = _cellsize;
}

void Grid3DObject::setValues(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _cellsize, const Coords& _llcorner)
{
	setValues(_ncols, _nrows, _ndepth, _cellsize);
	llcorner = _llcorner;
}

bool Grid3DObject::isSameGeolocalization(const Grid3DObject& target)
{
	if( ncols==target.ncols && nrows==target.nrows && ndepth==target.ndepth &&
		cellsize==target.cellsize && llcorner==target.llcorner) {
		return true;
	} else {
		return false;
	}
}

#ifdef _POPC_
#include "marshal_meteoio.h"
void Grid3DObject::Serialize(POPBuffer &buf, bool pack)
{
	if (pack)
	{
		buf.Pack(&ncols,1);
		buf.Pack(&nrows,1);
		buf.Pack(&ndepth,1);
		buf.Pack(&cellsize,1);
		marshal_Coords(buf, llcorner, 0, FLAG_MARSHAL, NULL);
		//unsigned int x,y,z;
		//grid3D.size(x,y,z);
		marshal_TYPE_DOUBLE3D(buf, grid3D, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&ncols,1);
		buf.UnPack(&nrows,1);
		buf.UnPack(&ndepth,1);
		buf.UnPack(&cellsize,1);
		marshal_Coords(buf, llcorner, 0, !FLAG_MARSHAL, NULL);
		//grid3D.clear();//if(grid2D!=NULL)delete(grid2D);
		marshal_TYPE_DOUBLE3D(buf, grid3D, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

