/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/IOUtils.h>
#include <meteoio/ResamplingAlgorithms2D.h>
#include <cmath>

using namespace std;

namespace mio {

/**
 * @brief Bilinear spatial data resampling
 */
const Grid2DObject ResamplingAlgorithms2D::BilinearResampling(const Grid2DObject &i_grid, const double &factor)
{
	const double cellsize = i_grid.cellsize/factor;
	const unsigned int ncols = (unsigned int)round( i_grid.ncols*factor );
	const unsigned int nrows = (unsigned int)round( i_grid.nrows*factor );
	Grid2DObject o_grid(ncols, nrows, cellsize, i_grid.llcorner);

	Bilinear_nodata(o_grid, i_grid); //GridObjects always keep nodata
	return o_grid;
}

const Grid2DObject ResamplingAlgorithms2D::NearestNeighbour(const Grid2DObject &i_grid, const double &factor)
{
	const double cellsize = i_grid.cellsize/factor;
	const unsigned int ncols = (unsigned int)round( i_grid.ncols*factor );
	const unsigned int nrows = (unsigned int)round( i_grid.nrows*factor );
	Grid2DObject o_grid(ncols, nrows, cellsize, i_grid.llcorner);

	NearestNeighbour(o_grid, i_grid); //GridObjects always keep nodata
	return o_grid;
}

///////////////////////////////////////////////////////////////////////
//Private Methods
///////////////////////////////////////////////////////////////////////
void ResamplingAlgorithms2D::NearestNeighbour(Grid2DObject &o_grid, const Grid2DObject &i_grid)
{
	const unsigned int org_ncols = i_grid.ncols;
	const unsigned int org_nrows = i_grid.nrows;
	const double scale_x = (double)o_grid.ncols / (double)org_ncols;
	const double scale_y = (double)o_grid.nrows / (double)org_nrows;

	for (unsigned int jj=0; jj<o_grid.nrows; jj++) {
		unsigned int org_jj = (unsigned int) round( (double)jj/scale_y );
		if(org_jj>=org_nrows) org_jj=org_nrows-1;

		for (unsigned int ii=0; ii<o_grid.ncols; ii++) {
			unsigned int org_ii = (unsigned int) round( (double)ii/scale_x );
			if(org_ii>=org_ncols) org_ii=org_ncols-1;
			o_grid(ii,jj) = i_grid(org_ii, org_jj);
		}
	}
}

void ResamplingAlgorithms2D::Bilinear_raw(Grid2DObject &o_grid, const Grid2DObject &i_grid)
{
	const unsigned int org_ncols = i_grid.ncols;
	const unsigned int org_nrows = i_grid.nrows;
	const double scale_x = (double)o_grid.ncols / (double)org_ncols;
	const double scale_y = (double)o_grid.nrows / (double)org_nrows;

	for (unsigned int jj=0; jj<o_grid.nrows; jj++) {
		const double org_y = (double)jj/scale_y;
		const unsigned int org_jj = (unsigned int) floor(org_y);
		const double y = org_y - (double)org_jj; //normalized y, between 0 and 1

		for (unsigned int ii=0; ii<o_grid.ncols; ii++) {
			const double org_x = (double)ii/scale_x;
			const unsigned int org_ii = (unsigned int) floor(org_x);
			const double x = org_x - (double)org_ii; //normalized x, between 0 and 1

			if(org_jj>=(org_nrows-1) || org_ii>=(org_ncols-1)) {
				o_grid(ii,jj) = i_grid(org_ii, org_jj);
				continue;
			}

			const double f_0_0 = i_grid(org_ii, org_jj);
			const double f_1_0 = i_grid(org_ii+1, org_jj);
			const double f_0_1 = i_grid(org_ii, org_jj+1);
			const double f_1_1 = i_grid(org_ii+1, org_jj+1);

			o_grid(ii,jj) = f_0_0 * (1.-x)*(1.-y) + f_1_0 * x*(1.-y) + f_0_1 * (1.-x)*y + f_1_1 *x*y;
		}
	}
}

void ResamplingAlgorithms2D::Bilinear_nodata(Grid2DObject &o_grid, const Grid2DObject &i_grid)
{
	const unsigned int org_ncols = i_grid.ncols;
	const unsigned int org_nrows = i_grid.nrows;
	const double scale_x = (double)o_grid.ncols / (double)org_ncols;
	const double scale_y = (double)o_grid.nrows / (double)org_nrows;

	for (unsigned int jj=0; jj<o_grid.nrows; jj++) {
		const double org_y = (double)jj/scale_y;
		const unsigned int org_jj = static_cast<unsigned int>( org_y );
		const double y = org_y - (double)org_jj; //normalized y, between 0 and 1

		for (unsigned int ii=0; ii<o_grid.ncols; ii++) {
			const double org_x = (double)ii/scale_x;
			const unsigned int org_ii = static_cast<unsigned int>( org_x );
			const double x = org_x - (double)org_ii; //normalized x, between 0 and 1

			if(org_jj>=(org_nrows-1) || org_ii>=(org_ncols-1)) {
				o_grid(ii,jj) = i_grid(org_ii, org_jj);
				continue;
			}

			const double f_0_0 = i_grid(org_ii, org_jj);
			const double f_1_0 = i_grid(org_ii+1, org_jj);
			const double f_0_1 = i_grid(org_ii, org_jj+1);
			const double f_1_1 = i_grid(org_ii+1, org_jj+1);

			double avg_value = 0.;
			unsigned int avg_count = 0;
			if(f_0_0!=IOUtils::nodata) {
				avg_value += f_0_0;
				avg_count++;
			}
			if(f_1_0!=IOUtils::nodata) {
				avg_value += f_1_0;
				avg_count++;
			}
			if(f_0_1!=IOUtils::nodata) {
				avg_value += f_0_1;
				avg_count++;
			}
			if(f_1_1!=IOUtils::nodata) {
				avg_value += f_1_1;
				avg_count++;
			}

			if(avg_count==4) {
				o_grid(ii,jj) = f_0_0 * (1.-x)*(1.-y) + f_1_0 * x*(1.-y) + f_0_1 * (1.-x)*y + f_1_1 *x*y;
				continue;
			}

			//special cases: less than two neighbours or three neighbours
			if(avg_count<=2) {
				o_grid(ii,jj) = IOUtils::nodata;
				continue;
			}

			double value = 0.;
			const double avg = avg_value/(double)avg_count;
			if(f_0_0!=IOUtils::nodata) value += f_0_0 * (1.-x)*(1.-y);
			else value += avg * (1.-x)*(1.-y);
			if(f_1_0!=IOUtils::nodata) value += f_1_0 * x*(1.-y);
			else value += avg * x*(1.-y);
			if(f_0_1!=IOUtils::nodata) value += f_0_1 * (1.-x)*y;
			else value += avg * (1.-x)*y;
			if(f_1_1!=IOUtils::nodata) value += f_1_1 *x*y;
			else value += avg *x*y;

			o_grid(ii,jj) = value;
		}
	}
}


} //namespace

