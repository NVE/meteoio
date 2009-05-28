/*
* PERRY BARTELT and MICHAEL LEHNING
* SLF-DAVOS/WSL-BIRMENSDORF
*/

/* **********************************************************************************************/
/*                                        VERSION 9.x                                          */
/*                               Derived from RESEARCH VERSION 9.0                             */
/* **********************************************************************************************/
/* **********************************************************************************/
/*  Copyright WSL Institute for Snow and Avalanche Research    SLF-DAVOS           */
/* **********************************************************************************/
/* This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/

//HACK: removed #include "Constants.h"
#include "Laws.h"

/**
 * @brief takes a long wave radiation and a temperature in K
 * @param lwr const double
 * @param T const double
 */
double lw_emissivity(const double lwr, const double T)
{
	
	return ( lwr / (STEFAN_BOLTZMANN * (T*T*T*T)) );
}

/**
* @brief Extrapolate air temperature
 * @param ta_ref const double
 * @param alti_ref const double
 * @param altitude const double
 * @return double
*/
double lw_TairLapseRate(const double ta_ref, const double alti_ref, const double altitude)
{
	double ta;
	const double lapse_rate = 0.0065; // K m-1
	
	ta = ta_ref - (altitude - alti_ref) * lapse_rate;
	
	return(ta);
}

/**
* @brief Equation for standard atmosphere; Stull "Meteorology for Scientists and Engineers", p.12-13
* standard pressure and standard temperature together with height as geopotential height
* with standard sea-level temperature 288.15 K and molecular scale temperature gradient
* -6.5 K/km
 * @param altitude const double
 * @return double
*/
double lw_AirPressure(const double altitude)
{
	double p;
	const double p0 = 101325.; 		// Air and standard pressure in Pa
	const double lapse_rate = 0.0065;	// K m-1
	const double sea_level_temp = 288.15;	// K
	const double expo = GRAVITY / (lapse_rate * GAS_CONSTANT_AIR);
	const double R0 = 6356766.0;		// Earth's radius in m
	
	p = p0 * pow( 1. - ( (lapse_rate * R0 * altitude) / (sea_level_temp * (R0 + altitude)) ), expo );
	
	return(p);
}

double lw_WetBulbTemperature(const double L, const double T, const double RH, const double altitude)
{
	// L:	latent heat of vaporisation
	const double mixing_ratio = 0.622;	//gas_constant(dry_air) / gas_constant(water_vapor)
	const double p = lw_AirPressure(altitude);
	const double Vp = lw_SaturationPressure(T);
	
	return ( T - (RH*Vp - Vp) * mixing_ratio * L / p / SPECIFIC_HEAT_AIR );

} // End of lw_WetBulbTemperature

/**
* @brief SATURATION PRESSURE
* This function finds the saturation pressure (in Pa). The saturation pressures are
* determined ( approximated ) using the formula:
*  e  = p0 * exp ( Qlh * ( T - T0 ) / ( R * T * T0 ) )
*  where
*  p0, T0:	triple point pressure and temperature of water
*  Qlh:	latent heat of evaporation, (water) or sublimation (ice)
*  @param T: absolute temperature (const double)
*/
double lw_SaturationPressure(const double T)
{
	double exp_p_sat; // exponent
	const double p0 = 610.78; // triple point pressure of water
	double c2, c3; // varying constants
	
	/*
	* computation of saturation vapor pressure in Pa with empirical formulas of "Beziehung von Murray"
	* for flat ice or water surfaces (others: Magnus formula, the precise formula after Golff-Gratch)
	* air temperature applied in K; note: the zero point of the Celsius scale of Murray is at 273.16 K!!!
	* from F.W. Murray (1967), J. of Applied Meteorology (transformed equation of Teten (1930) who gave the Magnus 
	* formula)
	*/
	
	if ( T < 273.16 ) { // for a flat ice surface
		c2 = 21.88;
		c3 = 7.66;
	} else { // for a flat water surface
		c2 = 17.27;
		c3 = 35.86;
	}
	exp_p_sat = c2 *  (T - 273.16) / (T - c3);
	
	return( p0 * exp( exp_p_sat ) );
} // End of lw_SaturationPressure


double RhtoDewPoint(double RH, double TA)
{
	//Convert a Relative Humidity into a dew point temperature
	//TA is in Kelvins, RH between 0 and 1, returns Td in Kelvins
	TA = K_TO_C(TA);
	double Es, E, Tdw, Tdi; //saturation and current water vapor pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97;	//parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55;	//parameters for ice
	const double Tfreeze = 0.;			//freezing temperature
	const double Tnucl = -16.0;			//nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);		//distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6);	//distance to pure water

	//in order to avoid getting NaN if RH=0
	RH += 0.0001;
	if (TA >= Tfreeze) {//above freezing point, water
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

double DewPointtoRh(double TD, double TA)
{
	//Convert a dew point temperature into a Relative Humidity
	//TA, TD are in Kelvins, RH is returned between 0 and 1
	TA = K_TO_C(TA);
	TD = K_TO_C(TD);
	double Es, E, Rhi, Rhw, Rh;			//saturation and current water vapro pressure
	const double Aw = 611.21, Bw = 17.502, Cw = 240.97;	//parameters for water
	const double Ai = 611.15, Bi = 22.452, Ci = 272.55;	//parameters for ice
	const double Tfreeze = 0.;			//freezing temperature
	const double Tnucl = -16.0;			//nucleation temperature
	const double di = 1. / ((TA - Tnucl) * (TA - Tnucl) + 1e-6);		//distance to pure ice
	const double dw = 1. / ((Tfreeze - TA) * (Tfreeze - TA) + 1e-6);	//distance to pure water

	if (TA >= Tfreeze) {
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

/*
* END of Laws.c
*/
