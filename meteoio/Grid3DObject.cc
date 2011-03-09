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
#include <meteoio/Grid3DObject.h>
#include <cmath>

using namespace std;

namespace mio {

Grid3DObject& Grid3DObject::operator=(const Grid3DObject& source) {
	if(this != &source) {
		grid3D = source.grid3D;
		ncols = source.ncols;
		nrows = source.nrows;
		ndepth = source.ndepth;
		cellsize = source.cellsize;
		llcorner = source.llcorner;
	}
	return *this;
}

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
				llcorner.getNorthing()+_ny*_grid3Dobj.cellsize,
				llcorner.getAltitude()+_nz*_grid3Dobj.cellsize );
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

bool Grid3DObject::gridify(std::vector<Coords>& vec_points) const {
	bool status=true;

	std::vector<Coords>::iterator v_Itr = vec_points.begin();
	while ( v_Itr != vec_points.end() ) {
		if( gridify(*v_Itr)==false ) {
			v_Itr = vec_points.erase(v_Itr);
			status=false;
		} else {
			v_Itr++;
		}
	}

	return status;
}

bool Grid3DObject::gridify(Coords& point) const {
	std::string proj_type, proj_args;
	point.getProj(proj_type, proj_args);
	if(proj_type=="NULL") {
		//if the projection was "NULL", we set it to the grid's
		point.copyProj(llcorner);
	}

	if(point.getGridI()!=IOUtils::inodata && point.getGridJ()!=IOUtils::inodata && point.getGridK()!=IOUtils::inodata) {
		//we need to compute (easting,northing) and (lat,lon) and altitude
		return( grid_to_WGS84(point) );
	} else {
		//we need to compute (i,j,k)
		return( WGS84_to_grid(point) );
	}
}

bool Grid3DObject::grid_to_WGS84(Coords& point) const {
	int i=point.getGridI(), j=point.getGridJ(), k=point.getGridK();

	if(i==IOUtils::inodata || j==IOUtils::inodata || k==IOUtils::inodata) {
		//the point is invalid (outside the grid or contains nodata)
		return false;
	}

	if(i>(signed)ncols || i<0 || j>(signed)nrows || j<0 || k>(signed)ndepth || k<0) {
		//the point is outside the grid, we reset the indices to the closest values
		//still fitting in the grid and return an error
		if(i<0) i=0;
		if(j<0) j=0;
		if(k<0) k=0;
		if(i>(signed)ncols) i=(signed)ncols;
		if(j>(signed)nrows) j=(signed)nrows;
		if(k>(signed)ndepth) k=(signed)ndepth;
		point.setGridIndex(i, j, k, false);
		return false;
	}

	//easting and northing in the grid's projection
	const double easting = ((double)i) * cellsize + llcorner.getEasting();
	const double northing = ((double)j) * cellsize + llcorner.getNorthing();
	const double altitude = ((double)k) * cellsize + llcorner.getAltitude();

	if(point.isSameProj(llcorner)==true) {
		//same projection between the grid and the point -> precise, simple and efficient arithmetics
		point.setXY(easting, northing, altitude);
	} else {
		//projections are different, so we have to do an intermediate step...
		Coords tmp_proj;
		tmp_proj.copyProj(point); //making a copy of the original projection
		point.copyProj(llcorner); //taking the grid's projection
		point.setXY(easting, northing, altitude);
		point.copyProj(tmp_proj); //back to the original projection -> reproject the coordinates
	}
	return true;
}

bool Grid3DObject::WGS84_to_grid(Coords point) const {
	if(point.getLat()==IOUtils::nodata || point.getLon()==IOUtils::nodata || point.getAltitude()==IOUtils::nodata) {
			//if the point is invalid, there is nothing we can do
			return false;
	}

	bool error_code=true;
	int i,j,k;

	if(point.isSameProj(llcorner)==true) {
		//same projection between the grid and the point -> precise, simple and efficient arithmetics
		i = (int)floor( (point.getEasting()-llcorner.getEasting()) / cellsize );
		j = (int)floor( (point.getNorthing()-llcorner.getNorthing()) / cellsize );
		k = (int)floor( (point.getAltitude()-llcorner.getAltitude()) / cellsize );
	} else {
		//projections are different, so we have to do an intermediate step...
		Coords tmp_point(point);
		tmp_point.copyProj(llcorner); //getting the east/north coordinates in the grid's projection
		i = (int)floor( (tmp_point.getEasting()-llcorner.getEasting()) / cellsize );
		j = (int)floor( (tmp_point.getNorthing()-llcorner.getNorthing()) / cellsize );
		k = (int)floor( (point.getAltitude()-llcorner.getAltitude()) / cellsize );
	}

	//checking that the calculated indices fit in the grid2D
	//and giving them the closest value within the grid if not.
	if(i<0) {
		i=0;
		error_code=false;
	}
	if(i>(signed)ncols) {
		i=(signed)ncols;
		error_code=false;
	}
	if(j<0) {
		j=0;
		error_code=false;
	}
	if(j>(signed)nrows) {
		j=(signed)nrows;
		error_code=false;
	}
	if(k<0) {
		k=0;
		error_code=false;
	}
	if(k>(signed)ndepth) {
		k=(signed)ndepth;
		error_code=false;
	}

	point.setGridIndex(i, j, k, false);
	return error_code;
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

void Grid3DObject::extractLayer(const unsigned int& z, Grid2DObject& layer) 
{
	layer.set(ncols, nrows, cellsize, llcorner);
	for(unsigned int jj=0; jj<nrows; jj++) {
		for(unsigned int ii=0; ii<ncols; ii++) {
			layer.grid2D(ii,jj) = grid3D(ii,jj,z);
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Grid3DObject& grid)
{
	os << "<Grid3DObject>\n";
	os << grid.llcorner;
	os << grid.ncols << " x " << grid.nrows  << " x " << grid.ndepth << " @ " << grid.cellsize << "m\n";
	os << grid.grid3D;
	os << "</Grid3DObject>\n";
	return os;
}

} //end namespace

#ifdef _POPC_
#include "marshal_meteoio.h"
using namespace mio; //HACK for POPC
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
		marshal_DOUBLE3D(buf, grid3D, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&ncols,1);
		buf.UnPack(&nrows,1);
		buf.UnPack(&ndepth,1);
		buf.UnPack(&cellsize,1);
		marshal_Coords(buf, llcorner, 0, !FLAG_MARSHAL, NULL);
		//grid3D.clear();//if(grid2D!=NULL)delete(grid2D);
		marshal_DOUBLE3D(buf, grid3D, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

//} //namespace //HACK for POPC
