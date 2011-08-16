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
#ifndef __SUN_H__
#define __SUN_H__

/*#include <ostream>
#include <iostream>
#include <iomanip>
#include <sstream>*/

#include <meteoio/IOUtils.h>
#include <meteoio/meteolaws/Suntrajectory.h>

namespace mio {

/**
 * @class SunObject
 * @brief A class to calculate Sun's characteristics
 *
 * @ingroup meteolaws
 * @author Mathias Bavay
 * @date   2010-06-10
 */
class SunObject {
	public:
		/// this enum provides the different position algorithms available
		typedef enum POSITION_ALGO {
			MEEUS ///<Jean Meeus' algorithm (Meeus, j. "Astronomical Algorithms", second edition, 1998, Willmann-Bell, Inc., Richmond, VA, USA)
		} position_algo;

		SunObject(const position_algo alg=MEEUS);
		SunObject(const double& _latitude, const double& _longitude, const double& _altitude);
		SunObject(const double& _latitude, const double& _longitude, const double& _altitude, const double& _julian, const double& TZ=0.);

		//local julian date and timezone
		void setDate(const double& _julian, const double& TZ=0.);
		void setLatLon(const double& _latitude, const double& _longitude, const double& _altitude);
		void setElevationThresh(const double& _elevation_threshold);

		void calculateRadiation(const double& ta, const double& rh, const double& pressure, const double& mean_albedo);
		void calculateRadiation(const double& ta, const double& rh, const double& mean_albedo);
		void getBeamRadiation(double& R_toa, double& R_direct, double& R_diffuse) const;
		void getHorizontalRadiation(double& R_toa, double& R_direct, double& R_diffuse) const;
		void getSlopeRadiation(const double& slope_azi, const double& slope_elev, double& R_toa, double& R_direct, double& R_diffuse) const;
		double getElevationThresh() const;

		double getSplitting(const double& iswr_modeled, const double& iswr_measured) const;

		//SunTrajectory position;
		SunMeeus position;

		double getJulian(const double& TZ);

		friend std::ostream& operator<<(std::ostream& os, const SunObject& data);
	private:
		void update();
		void getBeamPotential(const double& sun_elevation, const double& Eccentricity_corr,
		                  const double& ta, const double& rh, const double& pressure, const double& mean_albedo,
		                  double& R_toa, double& R_direct, double& R_diffuse);

	private:
		double julian_gmt;
		double latitude, longitude, altitude;
		double elevation_threshold;
		double beam_toa, beam_direct, beam_diffuse;
		static const double elevation_dftlThreshold;
};

} //end namespace

#endif
