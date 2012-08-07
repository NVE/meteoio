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
		ndepths = source.ndepths;
		cellsize = source.cellsize;
		llcorner = source.llcorner;
		z = source.z;
		z_is_absolute = source.z_is_absolute;
	}
	return *this;
}

Grid3DObject::Grid3DObject() : grid3D() //using Array3D default constructor
{
	ncols = nrows = ndepths = 0;
	cellsize = 0.0;
	z_is_absolute = true;
}

Grid3DObject::Grid3DObject(const Grid3DObject& i_grid3Dobj,
                           const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
                           const unsigned int& i_nwidths, const unsigned int& i_nheights, const unsigned int& i_ndepths)
	: grid3D(i_grid3Dobj.grid3D, i_nx,i_ny,i_nz, i_nwidths,i_nheights,i_ndepths)
{
	setValues(i_nwidths, i_nheights, i_ndepths, i_grid3Dobj.cellsize);
	z_is_absolute = true;

	//we take the previous corner (so we use the same projection parameters)
	//and we shift it by the correct X and Y distance
	llcorner = i_grid3Dobj.llcorner;
	if( (llcorner.getEasting()!=IOUtils::nodata) && (llcorner.getNorthing()!=IOUtils::nodata) ) {
		llcorner.setXY( llcorner.getEasting()+i_nx*i_grid3Dobj.cellsize,
		                llcorner.getNorthing()+i_ny*i_grid3Dobj.cellsize,
		                llcorner.getAltitude()+i_nz*i_grid3Dobj.cellsize );
	}
}

Grid3DObject::Grid3DObject(const unsigned int& i_ncols, const unsigned int& i_nrows, const unsigned int& i_ndepths,
                           const double& i_cellsize, const Coords& i_llcorner) : grid3D(i_ncols, i_nrows, i_ndepths, IOUtils::nodata)
{
	setValues(i_ncols, i_nrows, i_ndepths, i_cellsize, i_llcorner);
	z_is_absolute = true;
}

Grid3DObject::Grid3DObject(const unsigned int& i_ncols, const unsigned int& i_nrows, const unsigned int& i_ndepths,
                           const double& i_cellsize, const Coords& i_llcorner, const Array3D<double>& i_grid3D) : grid3D()
{
	set(i_ncols, i_nrows, i_ndepths, i_cellsize, i_llcorner, i_grid3D);
	z_is_absolute = true;
}

bool Grid3DObject::gridify(std::vector<Coords>& vec_points) const
{
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

bool Grid3DObject::gridify(Coords& point) const
{
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

bool Grid3DObject::grid_to_WGS84(Coords& point) const
{
	int i=point.getGridI(), j=point.getGridJ(), k=point.getGridK();

	if(i==IOUtils::inodata || j==IOUtils::inodata || k==IOUtils::inodata) {
		//the point is invalid (outside the grid or contains nodata)
		return false;
	}

	if(i>(signed)ncols || i<0 || j>(signed)nrows || j<0 || k>(signed)ndepths || k<0) {
		//the point is outside the grid, we reset the indices to the closest values
		//still fitting in the grid and return an error
		if(i<0) i=0;
		if(j<0) j=0;
		if(k<0) k=0;
		if(i>(signed)ncols) i=(signed)ncols;
		if(j>(signed)nrows) j=(signed)nrows;
		if(k>(signed)ndepths) k=(signed)ndepths;
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

bool Grid3DObject::WGS84_to_grid(Coords point) const
{
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
	if(k>(signed)ndepths) {
		k=(signed)ndepths;
		error_code=false;
	}

	point.setGridIndex(i, j, k, false);
	return error_code;
}

void Grid3DObject::set(const unsigned int& i_ncols, const unsigned int& i_nrows, const unsigned int& i_ndepths,
                       const double& i_cellsize, const Coords& i_llcorner)
{
	setValues(i_ncols, i_nrows, i_ndepths, i_cellsize, i_llcorner);
	grid3D.resize(ncols, nrows, ndepths, IOUtils::nodata);
}

void Grid3DObject::set(const unsigned int& i_ncols, const unsigned int& i_nrows, const unsigned int& i_ndepths,
                       const double& i_cellsize, const Coords& i_llcorner, const Array3D<double>& i_grid3D)
{
	//Test for equality in size: Only compatible Array3D<double> grids are permitted
	unsigned int nx, ny, nz;
	i_grid3D.size(nx, ny, nz);
	if ((i_ncols != nx) || (i_nrows != ny) || (i_ndepths != nz)) {
		throw IOException("Mismatch in size of Array3D<double> parameter grid3D and size of Grid3DObject", AT);
	}

	setValues(i_ncols, i_nrows, i_ndepths, i_cellsize, i_llcorner);
	grid3D = i_grid3D; //copy by value
}

void Grid3DObject::size(unsigned int& o_ncols, unsigned int& o_nrows, unsigned int& o_ndepths) const
{
	o_ncols = ncols;
	o_nrows = nrows;
	o_ndepths = ndepths;
}

unsigned int Grid3DObject::getNx() const {
	return ncols;
}

unsigned int Grid3DObject::getNy() const {
	return nrows;
}

unsigned int Grid3DObject::getNz() const {
	return ndepths;
}

void Grid3DObject::clear() {
	grid3D.clear();
	ncols = nrows = ndepths = 0;
}

bool Grid3DObject::isEmpty() const {
	return (ncols==0 && nrows==0 && ndepths==0);
}

void Grid3DObject::setValues(const unsigned int& i_ncols, const unsigned int& i_nrows, const unsigned int& i_ndepths,
                             const double& i_cellsize)
{
	ncols = i_ncols;
	nrows = i_nrows;
	ndepths = i_ndepths;
	cellsize = i_cellsize;
}

void Grid3DObject::setValues(const unsigned int& i_ncols, const unsigned int& i_nrows, const unsigned int& i_ndepths,
                             const double& i_cellsize, const Coords& i_llcorner)
{
	setValues(i_ncols, i_nrows, i_ndepths, i_cellsize);
	llcorner = i_llcorner;
}

bool Grid3DObject::isSameGeolocalization(const Grid3DObject& target)
{
	if( ncols==target.ncols && nrows==target.nrows && ndepths==target.ndepths &&
		cellsize==target.cellsize && llcorner==target.llcorner) {
		return true;
	} else {
		return false;
	}
}

void Grid3DObject::extractLayer(const unsigned int& i_z, Grid2DObject& layer)
{
	layer.set(ncols, nrows, cellsize, llcorner);
	for(unsigned int jj=0; jj<nrows; jj++) {
		for(unsigned int ii=0; ii<ncols; ii++) {
			layer.grid2D(ii,jj) = grid3D(ii,jj,i_z);
		}
	}
}

double& Grid3DObject::operator()(const unsigned int& ix, const unsigned int& iy, const unsigned int& iz) {
	return grid3D(ix,iy,iz);
}

double Grid3DObject::operator()(const unsigned int& ix, const unsigned int& iy, const unsigned int& iz) const {
	return grid3D(ix,iy,iz);
}

double& Grid3DObject::operator()(const unsigned int& i) {
	return grid3D(i);
}

double Grid3DObject::operator()(const unsigned int& i) const {
	return grid3D(i);
}

std::ostream& operator<<(std::ostream& os, const Grid3DObject& grid)
{
	os << "<Grid3DObject>\n";
	os << grid.llcorner;
	os << grid.ncols << " x " << grid.nrows  << " x " << grid.ndepths << " @ " << grid.cellsize << "m\n";
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
		buf.Pack(&ndepths,1);
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
		buf.UnPack(&ndepths,1);
		buf.UnPack(&cellsize,1);
		marshal_Coords(buf, llcorner, 0, !FLAG_MARSHAL, NULL);
		//grid3D.clear();//if(grid2D!=NULL)delete(grid2D);
		marshal_DOUBLE3D(buf, grid3D, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

//} //namespace //HACK for POPC
