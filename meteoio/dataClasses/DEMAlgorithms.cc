// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2018 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <cmath>
#include <limits.h>
#include <algorithm>
#include <fstream>
#include <sstream>
//#include <cerrno>
#include <cstring>

#include <meteoio/dataClasses/DEMAlgorithms.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/MathOptim.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/meteoLaws/Meteoconst.h> //for math constants

/**
* @file DEMAlgorithms.cc
* @brief implementation of the static DEMAlgorithms class
*/

using namespace std;

namespace mio {

/**
* @brief Computes the hillshade for the dem
* This "fake illumination" method is used to better show the relief on maps.
* @param[in] dem DEM to work with
* @param[in] elev elevation (in degrees) of the source of light
* @param[in] azimuth azimuth (in degrees) of the source of light
* @return hillshade grid that containing the illumination
*
*/
Grid2DObject DEMAlgorithms::getHillshade(const DEMObject& dem, const double& elev, const double& azimuth)
{
	if (dem.slope.empty() || dem.azi.empty())
		throw InvalidArgumentException("Hillshade computation requires slope and azimuth!", AT);

	const double zenith_rad = (90.-elev)*Cst::to_rad;
	const double azimuth_rad = azimuth*Cst::to_rad;
	const size_t ncols = dem.getNx();
	const size_t nrows = dem.getNy();

	Grid2DObject hillshade(ncols, nrows, dem.cellsize, dem.llcorner);

	for (size_t j=0; j<nrows; j++) {
		for (size_t i=0; i<ncols; i++) {
			const double alt = dem.grid2D(i,j);
			const double sl = dem.slope(i,j);
			const double az = dem.azi(i,j);
			
			if (alt!=IOUtils::nodata && sl!=IOUtils::nodata && az!=IOUtils::nodata) {
				const double sl_rad = sl*Cst::to_rad;
				const double tmp = cos(zenith_rad) * cos(sl_rad) + sin(zenith_rad) * sin(sl_rad) * cos(azimuth_rad-az*Cst::to_rad);
				hillshade(i,j) = (tmp>=0.)? tmp : 0.;
			} else
				hillshade(i,j) = IOUtils::nodata;
		}
	}

	return hillshade;
}

double DEMAlgorithms::getSearchDistance(const DEMObject& dem)
{
	static const double sun_elev_thresh = 5.;
	if (dem.min_altitude!=IOUtils::nodata && dem.max_altitude!=IOUtils::nodata) {
		return max( (dem.max_altitude - dem.min_altitude) / tan(sun_elev_thresh*Cst::to_rad), dem.cellsize*1.5); //make sure we use at least 1 cell, even in the diagonal
	}
	
	return IOUtils::nodata;
}

/**
* @brief Returns the tangente of the horizon from a given point looking toward a given bearing
* @param[in] dem DEM to work with
* @param[in] ix1 x index of the origin point
* @param[in] iy1 y index of the origin point
* @param[in] bearing direction given by a compass bearing
* @return tangente of angle above the horizontal (in deg) or IOUtils::nodata if the point (ix1, iy1) does not fit within the provided DEM
*/
double DEMAlgorithms::getHorizon(const DEMObject& dem, const size_t& ix1, const size_t& iy1, const double& bearing)
{
	const double max_shade_distance = getSearchDistance(dem);
	if (max_shade_distance==IOUtils::nodata) 
		throw InvalidArgumentException("DEM not properly initialized or only filled with nodata", AT);
	
	const int dimx = (signed)dem.grid2D.getNx();
	const int dimy = (signed)dem.grid2D.getNy();
	if ((signed)ix1>dimx || (signed)iy1>dimy) return IOUtils::nodata; //in case the point does not fith within the provided DEM
	if (ix1==0 || (signed)ix1==dimx-1 || iy1==0 || (signed)iy1==dimy-1) return 0.; //a border cell is not shadded
	
	const double cell_alt = dem.grid2D(ix1, iy1);
	double horizon_tan_angle = 0.;
	const double sin_alpha = sin(bearing*Cst::to_rad);
	const double cos_alpha = cos(bearing*Cst::to_rad);
	
	size_t nb_cells = 1;
	bool horizon_found = false;
	while (!horizon_found) {
		nb_cells++;
		const int ix2 = (int)ix1 + (int)round( ((double)nb_cells)*sin_alpha ); //alpha is a bearing
		const int iy2 = (int)iy1 + (int)round( ((double)nb_cells)*cos_alpha ); //alpha is a bearing

		if (ix2<=0 || ix2>=dimx-1 || iy2<=0 || iy2>=dimy-1) break; //we are out of the dem

		const double new_altitude = dem.grid2D(ix2, iy2);
		if (new_altitude==mio::IOUtils::nodata) break; //we stop at nodata cells

		const double DeltaH = new_altitude - cell_alt;
		const double distance = sqrt( (double)( Optim::pow2(ix2-(signed)ix1) + Optim::pow2(iy2-(signed)iy1)) ) * dem.cellsize;
		const double tan_angle = DeltaH/distance;
		if (tan_angle>horizon_tan_angle) horizon_tan_angle = tan_angle;

		if (distance>max_shade_distance) horizon_found=true; //maximum lookup distance reached
	}

	return horizon_tan_angle;
}

/**
* @brief Returns the tangente of the horizon from a given point looking toward a given bearing
* @param[in] dem DEM to work with
* @param[in] point the origin point
* @param[in] bearing direction given by a compass bearing
* @return tangente of angle above the horizontal (in deg) or IOUtils::nodata if the provided point does not fir within the 
* provided DEM
*/
double DEMAlgorithms::getHorizon(const DEMObject& dem, Coords point, const double& bearing)
{
	if (!point.indexIsValid()) {
		if (!dem.gridify(point)) return IOUtils::nodata;
	}
	
	const size_t ix1 = static_cast<size_t>( point.getGridI() );
	const size_t iy1 = static_cast<size_t>( point.getGridJ() );
	return getHorizon(dem, ix1, iy1, bearing);
}

/**
* @brief Returns the horizon from a given point looking 360 degrees around by increments. If the provided point does not
* fit within the provided DEM, an empty result set is returned.
* @param[in] dem DEM to work with
* @param[in] point the origin point
* @param[in] increment to the bearing between two angles
* @return horizon vector of heights above a given angle
*/
std::vector< std::pair<double,double> > DEMAlgorithms::getHorizonScan(const DEMObject& dem, Coords point, const double& increment)
{
	std::vector< std::pair<double,double> > horizon;
	if (!point.indexIsValid()) {
		if (!dem.gridify(point)) return horizon;
	}
	
	const size_t ix1 = static_cast<size_t>( point.getGridI() );
	const size_t iy1 = static_cast<size_t>( point.getGridJ() );
	for (double bearing=0.0; bearing <360.; bearing += increment) {
		const double tan_alpha = getHorizon(dem, ix1, iy1, bearing);
		if (tan_alpha!=IOUtils::nodata)
			horizon.push_back( make_pair(bearing, atan(tan_alpha)*Cst::to_deg) );
	}
	
	return horizon;
}

//custom function for sorting the horizons
struct sort_horizons {
	bool operator()(const std::pair<double,double> &left, const std::pair<double,double> &right) {
		if (left.first < right.first) return true; else return false;
	}
};

std::vector< std::pair<double,double> > DEMAlgorithms::readHorizonScan(const std::string& where, const std::string& filename)
{
	std::vector< std::pair<double,double> > horizon;
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << where << ": " << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

	try {
		size_t lcount = 0;
		double azimuth, value;
		std::string line;
		do {
			lcount++;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss(line);
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);
			iss >> std::skipws >> azimuth;
			if ( !iss || azimuth<0. || azimuth>360.) {
				std::ostringstream ss;
				ss << "Invalid azimuth in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}
			iss >> std::skipws >> value;
			if ( !iss ){
				std::ostringstream ss;
				ss << "Invalid value in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}

			horizon.push_back( make_pair(azimuth, value) );
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
	
	if (horizon.empty()) throw InvalidArgumentException(where+", no valid mask found in file '"+filename+"'", AT);
	std::sort(horizon.begin(), horizon.end(), sort_horizons());
	return horizon;
}

//linear interpolation between the available points
double DEMAlgorithms::getMaskElevation(const std::vector< std::pair<double,double> > &mask, const double& azimuth)
{
	const std::vector< std::pair<double, double> >::const_iterator next = std::upper_bound(mask.begin(), mask.end(), make_pair(azimuth, 0.), sort_horizons()); //first element that is > azimuth
	
	double x1, y1, x2, y2;
	if (next!=mask.begin() && next!=mask.end()) { //normal case
		const size_t ii = next - mask.begin();
		x1 = mask[ii-1].first;
		y1 = mask[ii-1].second;
		x2 = mask[ii].first;
		y2 = mask[ii].second;
	} else {
		x1 = mask.back().first - 360.;
		y1 = mask.back().second;
		x2 = mask.front().first;
		y2 = mask.front().second;
	}
	
	const double a = (y2 - y1) / (x2 - x1);
	const double b = y2 - a * x2;
	
	return a*azimuth + b;
}

/**
 * @brief Compute the sky view factors for the terrain radiation based on the DEM.
 * This is inspired (ie with some changes) by Manners, J., S. B. Vosper, and N. Roberts, <i>"Radiative transfer over resolved
 * topographic features for high‐resolution weather prediction"</i>, Quarterly journal of the
 * royal meteorological society, <b>138.664</b>, pp720-733, 2012.
 *
 * @param[in] dem DEM to work with
 * @param[in] ii x coordinate of the cell whose view factor should be computed
 * @param[in] jj y coordinate of the cell whose view factor should be computed
 * @return sky view factor
 */

double DEMAlgorithms::getCellSkyViewFactor(const DEMObject& dem, const size_t& ii, const size_t& jj)
{
	const double tan_slope = tan( dem.slope(ii,jj)*Cst::to_rad );
	const double azi = dem.azi(ii,jj);
	const double max_shade_distance = getSearchDistance(dem);
	static const unsigned int nSectors = 32;

	double sum=0.;
	for (unsigned int sector=0; sector<nSectors; sector++) {
		const double bearing = 360. * (double)sector / (double)nSectors;
		const double cos_azi_diff = cos((bearing - azi)*Cst::to_rad);
		const double elev = atan( getTanMaxSlope(dem, max_shade_distance, bearing, ii, jj) );

		const double correction_horizon =  atan((tan_slope*cos_azi_diff));
		double new_horizon=elev +correction_horizon;
		if (new_horizon<0) new_horizon=0;

		const double sector_vf = Optim::pow2( sin(mio::Cst::PI2-new_horizon) );
		sum += sector_vf;
	}

	return sum / nSectors;
}

double DEMAlgorithms::getTanMaxSlope(const DEMObject& dem, const double& dmax, const double& bearing, const size_t& i, const size_t& j)
{
	const double inv_dmax = 1./dmax;
	const double sin_alpha = sin(bearing*Cst::to_rad);
	const double cos_alpha = cos(bearing*Cst::to_rad);
	const double ref_altitude = dem.grid2D(i, j);
	const double cellsize_sq = mio::Optim::pow2(dem.cellsize);
	const int ii = static_cast<int>(i), jj = static_cast<int>(j);
	const int ncols = static_cast<int>(dem.getNx()), nrows = static_cast<int>(dem.getNy());

	int ll=ii, mm=jj;

	double max_tan_slope = -99999.;
	size_t nb_cells = 0;
	while ( !(ll<0 || ll>ncols-1 || mm<0 || mm>nrows-1) ) {
		const double altitude = dem.grid2D(ll, mm);
		if ( (altitude!=mio::IOUtils::nodata) && !(ll==ii && mm==jj) ) {
			const double delta_elev = altitude - ref_altitude;
			const double inv_distance = Optim::invSqrt( cellsize_sq*(Optim::pow2(ll-ii) + Optim::pow2(mm-jj)) );
			if (inv_distance<inv_dmax) break; //stop if distance>dmax

			const double tan_slope = delta_elev*inv_distance;
			if ( tan_slope>max_tan_slope ) max_tan_slope = tan_slope;
		}

		//move to next cell
		nb_cells++;
		ll = ii + (int)round( ((double)nb_cells)*sin_alpha ); //alpha is a bearing
		mm = jj + (int)round( ((double)nb_cells)*cos_alpha ); //alpha is a bearing
	}

	return max_tan_slope;
}


} //end namespace
