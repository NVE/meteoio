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

#include <assert.h>
#include <cmath>
#include <meteoio/meteolaws/Atmosphere.h>
#include <meteoio/meteolaws/Meteoconst.h>
#include <meteoio/meteolaws/Sun.h>
#include <meteoio/IOUtils.h>

namespace mio {
/**
 * @brief Calculate the black body emissivity
 * @param lwr longwave radiation emitted by the body (W m-2)
 * @param T   surface temperature of the body (K)
 * @return black body emissivity (0-1)
 */
double Atmosphere::blkBody_Emissivity(const double& lwr, const double& T) {
	const double T2 = T*T;
	return ( lwr / (Cst::stefan_boltzmann * (T2*T2)) );
}

/**
 * @brief Calculates the black body long wave radiation knowing its emissivity
 * @param ea emissivity of the body (0-1)
 * @param T   surface temperature of the body (K)
 * @return black body radiation (W/m^2)
 */
double Atmosphere::blkBody_Radiation(const double& ea, const double& T) {
	const double T2 = T*T;
	return ( ea * (Cst::stefan_boltzmann * (T2*T2)) );
}

/**
* @brief Standard atmosphere pressure
* @param altitude altitude above sea level (m)
* @return standard pressure (Pa)
*/
double Atmosphere::stdAirPressure(const double& altitude) {
	const double p0 = Cst::std_press; // Air and standard pressure in Pa
	const double lapse_rate = 0.0065; // K m-1
	const double sea_level_temp = 288.15; // K
	const double expo = Cst::gravity / (lapse_rate * Cst::gaz_constant_dry_air);
	const double R0 = Cst::earth_R0; // Earth's radius in m

	const double p = p0 * pow( 1. - ( (lapse_rate * R0 * altitude) / (sea_level_temp * (R0 + altitude)) ), expo );

	return(p);
}

/**
* @brief Standard dry air pressure
* @param altitude altitude above sea level (m)
* @param temperature air temperature (K)
* @return standard pressure (Pa)
*/
double Atmosphere::stdDryAirDensity(const double& altitude, const double& temperature) {
	return stdAirPressure(altitude)/(Cst::gaz_constant_dry_air*temperature);
}

/**
* @brief Calculates the water vapor density, for a given temperature and vapor pressure
* @param Temperature air temperature (K)
* @param VaporPressure water vapor pressure (Pa)
* @return water vapor density (kg/m^3)
*/
double Atmosphere::waterVaporDensity(const double& Temperature, const double& VaporPressure) {
	return (Cst::water_molecular_mass*VaporPressure)/(Cst::gaz_constant*Temperature);
}

/**
* @brief Standard atmosphere wet bulb temperature.
* This gives the lowest temperature that could be reached by water evaporation. It is therefore linked to
* relative humidity. This implementation assumes a standard atmosphere for pressure and saturation pressure.
* @param T air temperature (K)
* @param RH relative humidity (between 0 and 1)
* @param altitude altitude above sea level (m)
* @return wet bulb temperature (K)
*/
double Atmosphere::wetBulbTemperature(const double& T, const double& RH, const double& altitude)
{
	const double L = Cst::l_water_vaporization; //latent heat of vaporisation
	const double mixing_ratio = Cst::gaz_constant_dry_air / Cst::gaz_constant_water_vapor;
	const double p = stdAirPressure(altitude);
	const double Vp = waterSaturationPressure(T);

	return ( T - (RH*Vp - Vp) * mixing_ratio * L / p / Cst::specific_heat_air );
}

/**
* @brief Wind log profile.
* This is used to compute the equivalent wind speed at a different height.
* @details It depends on the roughness length
* that can take a default value (0.03) or any other, for example:
* <center><table border="0">
* <tr><th>land class</th><th>roughness length</th></tr>
* <tr><td>sea</td><td>0.0002</td></tr>
* <tr><td>smooth</td><td>0.005</td></tr>
* <tr><td>open</td><td>0.03</td></tr>
* <tr><td>roughly open</td><td>0.10</td></tr>
* <tr><td>rough</td><td>0.25</td></tr>
* <tr><td>very rough</td><td>0.5</td></tr>
* <tr><td>closed</td><td>1.0</td></tr>
* <tr><td>chaotic</td><td>over 2.0</td></tr>
* </table></center>
*
* @param v_ref reference wind speed (m/s)
* @param z_ref height of the reference wind speed (m)
* @param z new height (m)
* @param z0 roughness length (m)
* @return wind speed at height z (m/s)
*/
double Atmosphere::windLogProfile(const double& v_ref, const double& z_ref, const double& z, const double& z0)
{
	return v_ref * ( log(z/z0) / log(z_ref/z0) );
}

/**
* @brief Wind chill temperature.
* This is an index aiming at expressing the human-percived feeling of air temperature on exposed skin due to wind.
* This is NOT a sceintific measurement, only an index to express a subjective feeling.
* It is inapplicable above 10 Celsius and a few m/s wind (5 m/s used here), therefore returning air temperature.
* @param TA air temperature (K)
* @param VW wind velocity (m/s)
* @return wind chill temperature (K)
*/
double Atmosphere::windChill(const double& TA, const double& VW)
{
	if(TA>(273.15+10.) || VW<5.) return TA; //not applicable in this case

	const double t = K_TO_C(TA); //in Celsius
	const double v = VW*3.6; //in km/h
	const double v016 = pow(v, 0.16);

	const double WCT = 13.12 + 0.6215*t - 11.37*v016 + 0.3965*t*v016;
	return C_TO_K(WCT);
}

/**
* @brief Heat index.
* This is an index aiming at expressing the human-perceived air temperature due to humidity.
* This is NOT a sceintific measurement, only an index to express a subjective feeling.
* It is inapplicable below 26.7 Celsius and below 40% humidity, therefore returning air temperature.
* @param TA air temperature (K)
* @param RH relative humidity (between 0 and 1)
* @return Heat index (K)
*/
double Atmosphere::heatIndex(const double& TA, const double& RH)
{
	if(TA<(273.15+26.7) || RH<0.4) return TA; //not applicable in this case

	const double t = K_TO_C(TA); //in Celsius
	const double t2 = t*t;
	const double rh = RH*100.; //in percent
	const double rh2 = rh*rh;

	const double c1=-8.784695, c2=1.61139411, c3=2.338549 , c4=-0.14611605;
	const double c5=-1.2308094e-2 , c6=-1.6424828e-2 , c7=2.211732e-3 , c8=7.2546e-4, c9=-3.582e-6;

	const double HI = c1 + c2*t + c3*rh + c4*t*rh + c5*t2 + c6*rh2 + c7*t2*rh + c8*t*rh2 + c9*t2*rh2;

	return C_TO_K(HI);
}

/**
* @brief Standard water vapor saturation
* @param T air temperature (K)
* @return standard water vapor saturation pressure (Pa)
*/
double Atmosphere::waterSaturationPressure(const double& T) {
	double c2, c3; // varying constants

	if ( T < 273.16 ) { // for a flat ice surface
		c2 = 21.88;
		c3 = 7.66;
	} else { // for a flat water surface
		c2 = 17.27;
		c3 = 35.86;
	}
	const double exp_p_sat = c2 *  (T - 273.16) / (T - c3); //exponent

	return( Cst::p_water_triple_pt * exp( exp_p_sat ) );
}

/**
* @brief Virtual temperature multiplying factor.
* In order to get a virtual temperature, multiply the original temperature by this factor. Note:
* since e/p is actually used, both pressures only need to have the same units.
* @param e vapor pressure (Pa)
* @param p air pressure (Pa)
* @return virtual temperature multiplying coefficient
*/
double Atmosphere::virtualTemperatureFactor(const double& e, const double& p) {
	const double epsilon = 0.622;
	return 1. / (1.-(1.-epsilon)*e/p);
}

/**
* @brief Evaluate the atmosphere emissivity from the water vapor pressure and cloudiness.
* This is according to A. Omstedt, "A coupled one-dimensional sea ice-ocean model applied to a semi-enclosed basin",
* Tellus, 42 A, 568-582, 1990, DOI:10.1034/j.1600-0870.1990.t01-3-00007.
* @param RH relative humidity (between 0 and 1)
* @param TA air temperature (K)
* @param cloudiness cloudiness (between 0 and 1, 0 being clear sky)
* @return emissivity (between 0 and 1)
*/
double Atmosphere::Omstedt_emissivity(const double& RH, const double& TA, const double& cloudiness) {
	const double e0 = RH * waterSaturationPressure(TA); //water vapor pressure
	const double eps_w = 0.97;
	const double a1 = 0.68;
	const double a2 = 0.0036;
	const double a3 = 0.18;

	const double ea = (eps_w * (a1 + a2 * sqrt(e0)) * (1. + a3 * cloudiness * cloudiness)); //emissivity
	if(ea > 1.0)
		return 1.0;

	return ea;
}

/**
* @brief Evaluate the long wave radiation from RH, TA and cloudiness.
* This is according to A. Omstedt, "A coupled one-dimensional sea ice-ocean model applied to a semi-enclosed basin",
* Tellus, 42 A, 568-582, 1990, DOI:10.1034/j.1600-0870.1990.t01-3-00007.
* @param RH relative humidity (between 0 and 1)
* @param TA air temperature (K)
* @param cloudiness cloudiness (between 0 and 1, 0 being clear sky)
* @return long wave radiation (W/m^2)
*/
double Atmosphere::Omstedt_ilwr(const double& RH, const double& TA, const double& cloudiness) {
	const double ea = Omstedt_emissivity(RH, TA, cloudiness);
	return blkBody_Radiation(ea, TA);
}

/**
 * @brief Evaluate the atmosphere emissivity for clear sky.
 * This uses the formula from Brutsaert -- "On a Derivable
 * Formula for Long-Wave Radiation From Clear Skies", Journal of Water Resources
 * Research, Vol. 11, No. 5, October 1975, pp 742-744.
 * Alternative: Satterlund (1979): Water Resources Research, 15, 1649-1650.
 * @param RH relative humidity (between 0 and 1)
 * @param TA Air temperature (K)
 * @return clear sky emissivity
 */
double Atmosphere::Brutsaert_emissivity(const double& RH, const double& TA) {
	const double e0 = RH * waterSaturationPressure(TA); //water vapor pressure
	const double e0_mBar = 0.01 * e0;
	const double exponent = 1./7.;

	return (1.24 * pow( (e0_mBar / TA), exponent) );
}

/**
 * @brief Evaluate the long wave radiation for clear sky.
 * This uses the formula from Brutsaert -- "On a Derivable
 * Formula for Long-Wave Radiation From Clear Skies", Journal of Water Resources
 * Research, Vol. 11, No. 5, October 1975, pp 742-744.
 * Alternative: Satterlund (1979): Water Resources Research, 15, 1649-1650.
 * @param RH relative humidity (between 0 and 1)
 * @param TA Air temperature (K)
 * @return long wave radiation (W/m^2)
*/
double Atmosphere::Brutsaert_ilwr(const double& RH, const double& TA) {
	const double ea = Brutsaert_emissivity(RH, TA);
	return blkBody_Radiation(ea, TA);
}

/**
 * @brief Evaluate the atmosphere emissivity for clear sky.
 * This uses the formula from Dilley and O'Brien -- "Estimating downward clear sky
 * long-wave irradiance at the surface from screen temperature and precipitable water",
 * Q. J. R. Meteorolo. Soc., Vol. 124, 1998, pp 1391-1401. The long wave is computed
 * and the ratio of this long wave to a black body emission gives an emissivity.
 * @param RH relative humidity (between 0 and 1)
 * @param TA near surface air temperature (K)
 * @return long wave radiation (W/m^2)
*/
double Atmosphere::Dilley_emissivity(const double& RH, const double& TA) {
	const double ilwr_dilley = Dilley_ilwr(RH, TA);
	const double ilwr_blkbody = blkBody_Radiation(1., TA);
	return ilwr_dilley/ilwr_blkbody;
}

/**
 * @brief Evaluate the long wave radiation for clear sky.
 * This uses the formula from Dilley and O'Brien -- "Estimating downward clear sky
 * long-wave irradiance at the surface from screen temperature and precipitable water",
 * Q. J. R. Meteorolo. Soc., Vol. 124, 1998, pp 1391-1401.
 * @param RH relative humidity (between 0 and 1)
 * @param TA near surface air temperature (K)
 * @return long wave radiation (W/m^2)
*/
double Atmosphere::Dilley_ilwr(const double& RH, const double& TA) {
	const double e0 = RH * waterSaturationPressure(TA) * 0.001; //water vapor pressure, kPa
	const double w = 4650.*e0/TA; //precipitable water, Prata 1996

	const double tmp = TA/Cst::t_water_triple_pt;
	const double pow_tmp6 = tmp*tmp*tmp*tmp*tmp*tmp;
	return 59.38 + 113.7*pow_tmp6 + 96.96*sqrt(w/25.);
}

/**
 * @brief Evaluate the atmosphere emissivity for clear sky.
 * This uses the formula from Prata -- "A new long-wave formula for estimating
 * downward clear-sky radiation at the surface", Q. J. R. Meteorolo. Soc., Vol. 122, 1996, pp 1127-1151.
 * @param RH relative humidity (between 0 and 1)
 * @param TA near surface air temperature (K)
 * @return long wave radiation (W/m^2)
*/
double Atmosphere::Prata_emissivity(const double& RH, const double& TA) {
	const double e0 = RH * waterSaturationPressure(TA) * 0.001; //water vapor pressure, kPa
	const double w = 4650.*e0/TA; //precipitable water, Prata 1996

	return 1. - (1.+w)*exp( -sqrt(1.2+3.*w) );
}

/**
 * @brief Evaluate the long wave radiation for clear sky.
 * This uses the formula from Prata -- "A new long-wave formula for estimating
 * downward clear-sky radiation at the surface", Q. J. R. Meteorolo. Soc., Vol. 122, 1996, pp 1127-1151.
 * @param RH relative humidity (between 0 and 1)
 * @param TA near surface air temperature (K)
 * @return long wave radiation (W/m^2)
*/
double Atmosphere::Prata_ilwr(const double& RH, const double& TA) {
	const double epsilon = Prata_emissivity(RH, TA);
	return blkBody_Radiation(epsilon, TA);
}

/**
 * @brief Evaluate the long wave radiation for clear or cloudy sky.
 * This uses the formula from Crawford and Duchon -- <i>"An Improved Parametrization
 * for Estimating Effective Atmospheric Emissivity for Use in Calculating Daytime
 * Downwelling Longwave Radiation"</i>, Journal of Applied Meteorology,
 * <b>38</b>, 1999, pp 474-480.
 * @param RH relative humidity (between 0 and 1)
 * @param TA Air temperature (K)
 * @param iswr_meas Measured Incoming Short Wave Radiation (W/m^2)
 * @param iswr_clear_sky Clear Sky Modelled Incoming Short Wave Radiation (W/m^2)
 * @param month current month (1-12, for a sinusoidal variation of the leading coefficients)
 * @return long wave radiation (W/m^2) or IOUtils::nodata at night time
*/
double Atmosphere::Crawford_ilwr(const double& RH, const double& TA, const double& iswr_meas, const double& iswr_clear_sky, const unsigned char& month) {
	if(iswr_meas<=0. || iswr_clear_sky<=0.)
		return IOUtils::nodata;

	double clf = 1. - iswr_meas/iswr_clear_sky; //cloud fraction estimate
	if(clf<0.) clf=0.;
	const double e = RH * waterSaturationPressure(TA); //near surface water vapor pressure
	const double e_mBar = 0.01 * e;

	const double epsilon = clf + (1.-clf) * (1.22 + 0.06*sin((month+2.)*Cst::PI/6.) ) * pow( (e_mBar/TA), 1./7.);
	return blkBody_Radiation(epsilon, TA);
}

/**
 * @brief Evaluate the long wave radiation for clear or cloudy sky.
 * This uses the formula from Crawford and Duchon -- <i>"An Improved Parametrization
 * for Estimating Effective Atmospheric Emissivity for Use in Calculating Daytime
 * Downwelling Longwave Radiation"</i>, Journal of Applied Meteorology,
 * <b>38</b>, 1999, pp 474-480.
 * @param lat latitude of the point of observation
 * @param lon longitude of the point of observation
 * @param altitude altitude of the point of observation
 * @param julian julian date at the point of observation
 * @param TZ time zone at the point of observation
 * @param RH relative humidity (between 0 and 1)
 * @param TA Air temperature (K)
 * @param ISWR Measured Incoming Short Wave Radiation (W/m^2)
 * @return long wave radiation (W/m^2) or IOUtils::nodata at night time
 * Please note that this call might NOT be efficient for processing large amounts of points,
 * since it internally builds complex objects at every call. You might want to copy/paste
 * its code in order to process data in bulk.
*/
double Atmosphere::Crawford_ilwr(const double& lat, const double& lon, const double& altitude,
                                 const double& julian, const double& TZ,
                                 const double& RH, const double& TA, const double& ISWR)
{
	if(TA==IOUtils::nodata || RH==IOUtils::nodata || ISWR==IOUtils::nodata) {
		return IOUtils::nodata;
	}
	SunObject Sun(lat, lon, altitude, julian, TZ);
	Sun.calculateRadiation(TA, RH, 0.5); //we force a terrain albedo of 0.5...
	double toa_h, direct_h, diffuse_h;
	Sun.getHorizontalRadiation(toa_h, direct_h, diffuse_h);

	Date date(julian, TZ, 0.);
	int year, month, day;
	date.getDate(year, month, day);

	return Atmosphere::Crawford_ilwr(RH, TA, ISWR, direct_h+diffuse_h, month);
}

/**
* @brief Convert a relative humidity to a dew point temperature.
* @param RH relative humidity between 0 and 1
* @param TA air temperature (K)
* @param force_water if set to true, compute over water. Otherwise, a smooth transition between over ice and over water is computed.
* @return dew point temperature (K)
*/
double Atmosphere::RhtoDewPoint(double RH, double TA, const bool& force_water)
{
	TA = K_TO_C(TA);
	double Es, E, Tdw, Tdi; //saturation and current water vapor pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97; //parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55; //parameters for ice
	const double Tfreeze = 0.;                          //freezing temperature
	const double Tnucl = -16.0;                         //nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);     //distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6); //distance to pure water

	//in order to avoid getting NaN if RH=0
	RH += 0.0001;
	assert(RH>0.);
	if (TA >= Tfreeze || force_water==true) {//above freezing point, water
		Es = Aw * exp( (Bw * TA) / (Cw + TA) );
		E = RH * Es;
		Tdw = ( Cw * log(E / Aw) ) / ( Bw - log(E / Aw) );
		return C_TO_K(Tdw);
	}
	if (TA < Tnucl) { //below nucleation, ice
		Es = Ai * exp( (Bi * TA) / (Ci + TA) );
		E = RH * Es;
		Tdi = ( Ci * log(E / Ai) ) / ( Bi - log(E / Ai) );
		return C_TO_K(Tdi);
	}

	//no clear state, we do a smooth interpolation between water and ice
	Es = Ai * exp( (Bi*TA) / (Ci + TA) );
	E = RH * Es;
	Tdi = ( Ci * log(E / Ai) ) / ( Bi - log(E / Ai) );

	Es = Aw * exp( (Bw * TA) / (Cw + TA) );
	E = RH * Es;
	Tdw = ( Cw * log(E / Aw) ) / ( Bw - log(E / Aw) );

	return C_TO_K( (di / (di + dw) * Tdi + dw / (di + dw) * Tdw) );
}

/**
* @brief Convert a dew point temperature to a relative humidity.
* @param TD dew point temperature (K)
* @param TA air temperature (K)
* @param force_water if set to true, compute over water. Otherwise, a smooth transition between over ice and over water is computed.
* @return relative humidity between 0 and 1
*/
double Atmosphere::DewPointtoRh(double TD, double TA, const bool& force_water)
{
	//Convert a dew point temperature into a Relative Humidity
	//TA, TD are in Kelvins, RH is returned between 0 and 1
	TA = K_TO_C(TA);
	TD = K_TO_C(TD);
	double Es, E, Rhi, Rhw, Rh;                         //saturation and current water vapro pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97; //parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55; //parameters for ice
	const double Tfreeze = 0.;                          //freezing temperature
	const double Tnucl = -16.0;                         //nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);     //distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6); //distance to pure water

	if (TA >= Tfreeze || force_water==true) {
		//above freezing point, water
		Es = Aw * exp( (Bw * TA) / (Cw + TA) );
		E  = Aw * exp( (Bw * TD) / (Cw + TD) );
		Rhw = (E / Es);
		if (Rhw > 1.) {
			return 1.;
		} else {
			return Rhw;
		}
	}
	if (TA < Tnucl) {
		//below nucleation, ice
		Es = Ai * exp( (Bi * TA) / (Ci + TA) );
		E  = Ai * exp( (Bi * TD) / (Ci + TD) );
		Rhi = (E / Es);
		if (Rhi > 1.) {
			return 1.;
		} else {
			return Rhi;
		}
	}

	//no clear state, we do a smooth interpolation between water and ice
	Es = Ai * exp( (Bi * TA) / (Ci + TA) );
	E  = Ai * exp( (Bi * TD) / (Ci + TD) );
	Rhi = E / Es;

	Es = Aw * exp( (Bw * TA) / (Cw + TA) );
	E  = Aw * exp( (Bw * TD) / (Cw + TD) );
	Rhw = E / Es;

	Rh = (di / (di + dw) * Rhi + dw / (di + dw) * Rhw);
	if(Rh > 1.) {
		return 1.;
	} else {
		return Rh;
	}
}

/**
* @brief Calculate the relative Humidity (RH) from specific humidity.
* @param altitude altitude over sea level (m)
* @param TA air temperature (K)
* @param qi specific humidity
* @return relative humidity between 0 and 1
*/
double Atmosphere::specToRelHumidity(const double& altitude, const double& TA, const double& qi)
{
	const double SatVaporDensity = waterVaporDensity(TA, waterSaturationPressure(TA));
	const double RH = (qi/(1.-qi))*stdDryAirDensity(altitude, TA)/SatVaporDensity;

	if(RH>1.) return 1.;
	else return RH;
}

/**
* @brief Calculate the specific Humidity from relative humidity (RH).
* @param altitude altitude over sea level (m)
* @param TA air temperature (K)
* @param RH relative humidity (between 0 and 1)
* @return specific humidity
*/
double Atmosphere::relToSpecHumidity(const double& altitude, const double& TA, const double& RH)
{
	const double dryAir_density = stdDryAirDensity(altitude, TA);
	const double wetAir_density = RH * waterVaporDensity(TA,waterSaturationPressure(TA));
	const double qi = 1./( dryAir_density/wetAir_density+1. );
	return qi;
}

} //namespace
