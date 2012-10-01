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
#ifndef __ATMOSPHERE_H__
#define __ATMOSPHERE_H__

namespace mio {

/**
 * @class Atmosphere
 * @brief A class to calculate the atmosphere's parameters
 *
 * @ingroup meteolaws
 * @author Mathias Bavay
 * @date   2010-06-10
 */
class Atmosphere {
	public:
		static double stdAirPressure(const double& altitude);
		static double stdDryAirDensity(const double& altitude, const double& temperature);
		static double waterSaturationPressure(const double& T);
		static double virtualTemperatureFactor(const double& e, const double& p);
		static double waterVaporDensity(const double& Temperature, const double& VaporPressure);
		static double wetBulbTemperature(const double& T, const double& RH, const double& altitude);
		static double windLogProfile(const double& v_ref, const double& z_ref, const double& z, const double& z0=0.03);

		static double windChill(const double& TA, const double& VW);
		static double heatIndex(const double& TA, const double& RH);

		static double Omstedt_emissivity(const double& RH, const double& TA, const double& cloudiness);
		static double Omstedt_ilwr(const double& RH, const double& TA, const double& cloudiness);
		static double Brutsaert_emissivity(const double& RH, const double& TA);
		static double Brutsaert_ilwr(const double& RH, const double& TA);
		static double Crawford_ilwr(const double& RH, const double& TA, const double& iswr_meas, const double& iswr_clear_sky, const unsigned char& month);
		static double Crawford_ilwr(const double& lat, const double& lon, const double& altitude,
		                            const double& julian, const double& TZ,
		                            const double& RH, const double& TA, const double& ISWR);
		static double Unsworth_ilwr(const double& RH, const double& TA, const double& iswr_meas, const double& iswr_clear_sky);
		static double Unsworth_ilwr(const double& lat, const double& lon, const double& altitude,
		                            const double& julian, const double& TZ,
		                            const double& RH, const double& TA, const double& ISWR);
		static double Dilley_emissivity(const double& RH, const double& TA);
		static double Dilley_ilwr(const double& RH, const double& TA);
		static double Prata_emissivity(const double& RH, const double& TA);
		static double Prata_ilwr(const double& RH, const double& TA);
		static double Kasten_cloudiness(const double& solarIndex);
		static double ILWR_parametrized(const double& lat, const double& lon, const double& altitude,
		                                const double& julian, const double& TZ,
		                                const double& RH, const double& TA, const double& ISWR, const double& cloudiness);

		static double RhtoDewPoint(double RH, double TA, const bool& force_water);
		static double DewPointtoRh(double TD, double TA, const bool& force_water);
		static double specToRelHumidity(const double& altitude, const double& TA, const double& qi);
		static double relToSpecHumidity(const double& altitude, const double& TA, const double& RH);

		static double blkBody_Emissivity(const double& lwr, const double& T);
		static double blkBody_Radiation(const double& ea, const double& T);
};

} //end namespace

#endif
