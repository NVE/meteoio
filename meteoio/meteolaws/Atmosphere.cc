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
#include <meteoio/IOUtils.h>

namespace mio {

/**
 * @brief Calculate the black body emissivity
 * @param lwr longwave radiation emitted by the body (W m-2)
 * @param T   surface temperature of the body (K)
 */
double Atmosphere::blkBody_Emissivity(const double& lwr, const double& T) {
	const double T2 = T*T;
	return ( lwr / (Cst::stefan_boltzmann * (T2*T2)) );
}

/**
 * @brief Calculates the black body long wave radiation knowing its emissivity
 * @param ea emissivity of the body (0-1)
 * @param T   surface temperature of the body (K)
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
{//HACK: should we check that RH in [0;1]?
	const double SatVaporDensity = waterVaporDensity(TA, waterSaturationPressure(TA));
	const double RH = (qi/(1.-qi))*stdDryAirDensity(altitude, TA)/SatVaporDensity;

	return RH;
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
