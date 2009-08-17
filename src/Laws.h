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

/*
* AUTHORS: MATHIAS BAVAY and CHARLES FIERZ
* SLF-DAVOS/WSL-BIRMENSDORF
*/

/*
* Compilation condition
*/
#ifndef LAWS_H
#define LAWS_H

#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

	//HACK: CONVERSIONS
	#define C_TO_K( T ) ( T + 273.15 )	  // Celsius to Kelvin
	#define K_TO_C( T ) ( T - 273.15 )	  // Kelvin to Celsius

	/*
	* CONSTANTS
	*/ //HACK: included some constants coming from Constants.h ...
	/*---------------------------------------------------------------------------------------------+ 
	 | Physical constants                                                                          |
	 +---------------------------------------------------------------------------------------------*/
	#define GRAVITY	9.80665		     // (m s-2)
	#define STEFAN_BOLTZMANN 5.67051e-8  // (W m-2 K-4)
	#define GAS_CONSTANT_AIR 287.	     // (J kg-1 K-1)   ( for air! )
	#define SPECIFIC_HEAT_AIR 1004.67    // see Stull "Meteorology for scientists and engineers" p44

	/*
	* STRUCTURES
	*/



// FUNCTION PROTOTYPES:
  double lw_emissivity(const double lwr, const double T);
  double lw_TairLapseRate(const double ta, const double ref_alti, const double altitude);
  double lw_AirPressure(const double altitude);
  double lw_WetBulbTemperature(const double L, const double T, const double RH, const double altitude);
  double lw_SaturationPressure(const double T);
  double RhtoDewPoint(double RH, double TA, const short int force_water);
  double DewPointtoRh(double TD, double TA, const short int force_water);
  double lw_LW_Brutsaert(const double e0, const double Ta);
  double lw_Omstedt(const double e0, const double cloud_frac);
  double lw_SnowResidualWaterContent(const double theta_i);
  double lw_AirEmissivity(const double input, const double Ta, const double Rh);

#ifdef __cplusplus
}
#endif

#endif
/*
* End of Laws.h
*/
