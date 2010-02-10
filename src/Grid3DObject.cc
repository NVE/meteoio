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
#include "MapProj.h"

Grid3DObject::Grid3DObject() : grid3D() //using Array3D default constructor
{
	ncols = nrows = ndepth = 0;
	xllcorner = yllcorner = cellsize = 0.0;
	latitude = longitude = IOUtils::nodata;	
}

Grid3DObject::Grid3DObject(const Grid3DObject& _grid3Dobj,
				const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
				const unsigned int& _nwidth, const unsigned int& _nheight, const unsigned int& _ndepth) 
	: grid3D(_grid3Dobj.grid3D, _nx,_ny,_nz, _nwidth,_nheight,_ndepth)
{
	setValues(_nwidth, _nheight, _ndepth, 
			(_grid3Dobj.xllcorner+_nx*_grid3Dobj.cellsize), (_grid3Dobj.yllcorner+_ny*_grid3Dobj.cellsize), 
			_grid3Dobj.latitude, _grid3Dobj.longitude, _grid3Dobj.cellsize);
}

Grid3DObject::Grid3DObject(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude,
				const double& _cellsize) : grid3D(_ncols, _nrows, _ndepth, IOUtils::nodata)
{
	setValues(_ncols, _nrows, _ndepth, _xllcorner, _yllcorner, _latitude, _longitude, _cellsize);
}

Grid3DObject::Grid3DObject(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude,
					  const double& _cellsize, const Array3D<double>& _grid3D) : grid3D()
{
	set(_ncols, _nrows, _ndepth, _xllcorner, _yllcorner, _latitude, _longitude, _cellsize, _grid3D);
}

void Grid3DObject::set(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude,
				const double& _cellsize)
{
	setValues(_ncols, _nrows, _ndepth, _xllcorner, _yllcorner, _latitude, _longitude, _cellsize);
	grid3D.resize(ncols, nrows, ndepth, IOUtils::nodata);	
}

void Grid3DObject::set(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude,
				const double& _cellsize, const Array3D<double>& _grid3D)
{
	//Test for equality in size: Only compatible Array3D<double> grids are permitted
	unsigned int nx, ny, nz;
	_grid3D.size(nx, ny, nz);
	if ((_ncols != nx) || (_nrows != ny) || (_ndepth != nz)) {
		throw IOException("Mismatch in size of Array3D<double> parameter grid3D and size of Grid3DObject", AT);
	}

	setValues(_ncols, _nrows, _ndepth, _xllcorner, _yllcorner, _latitude, _longitude, _cellsize);
	grid3D = _grid3D; //copy by value
}


void Grid3DObject::setValues(const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth,
				const double& _xllcorner, const double& _yllcorner,
				const double& _latitude, const double& _longitude, const double& _cellsize)
{
	ncols = _ncols;
	nrows = _nrows;
	ndepth = _ndepth;
	cellsize = _cellsize;
	xllcorner = _xllcorner;
	yllcorner = _yllcorner;
	latitude = _latitude;
	longitude = _longitude;

	//const MapProj dummy;
	//if(proj!=dummy) {
	//	checkCoordinates(proj);
	//}
}

//TODO: add WGS84_to_local methods to Grid3DObject (with/without terrain following coords)

void Grid3DObject::checkCoordinates(const MapProj& proj)
{
	//HACK: altitude is missing!!
	//calculate/check coordinates if necessary
	if(latitude==IOUtils::nodata || longitude==IOUtils::nodata) {
		if(xllcorner==IOUtils::nodata || yllcorner==IOUtils::nodata) {
			throw InvalidArgumentException("missing positional parameters (xll,yll) or (lat,long) for Grid3DObject", AT);
		}
		proj.convert_to_WGS84(xllcorner, yllcorner, latitude, longitude);
	} else {
		if(xllcorner==IOUtils::nodata || yllcorner==IOUtils::nodata) {
			proj.convert_from_WGS84(latitude, longitude, xllcorner, yllcorner);
		} else {
			double tmp_lat, tmp_lon;
			proj.convert_to_WGS84(xllcorner, yllcorner, tmp_lat, tmp_lon);
			if(!IOUtils::checkEpsilonEquality(latitude, tmp_lat, 1.e-4) || !IOUtils::checkEpsilonEquality(longitude, tmp_lon, 1.e-4)) {
				throw InvalidArgumentException("Latitude/longitude and xllcorner/yllcorner don't match for Grid3DObject", AT);
			}
		}
	}
}

bool Grid3DObject::isSameGeolocalization(const Grid3DObject& target)
{//HACK: altitude is missing!!
	if( ncols==target.ncols && nrows==target.nrows && ndepth==target.ndepth &&
		IOUtils::checkEpsilonEquality(latitude, target.latitude, 1.e-4) && 
		IOUtils::checkEpsilonEquality(longitude, target.longitude, 1.e-4) &&
		cellsize==target.cellsize) {
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
		buf.Pack(&xllcorner,1);
		buf.Pack(&yllcorner,1);
		buf.Pack(&latitude,1);
		buf.Pack(&longitude,1);
		buf.Pack(&cellsize,1);
		//unsigned int x,y,z;
		//grid3D.size(x,y,z);
		//marshal_TYPE_DOUBLE2D(buf, grid3D, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&ncols,1);
		buf.UnPack(&nrows,1);
		buf.UnPack(&ndepth,1);
		buf.UnPack(&xllcorner,1);
		buf.UnPack(&yllcorner,1);
		buf.UnPack(&latitude,1);
		buf.UnPack(&longitude,1);
		buf.UnPack(&cellsize,1);
		//grid3D.clear();//if(grid2D!=NULL)delete(grid2D);
		//marshal_TYPE_DOUBLE2D(buf, grid3D, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

