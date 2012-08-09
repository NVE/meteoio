/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __MATHOPTIM_H__
#define __MATHOPTIM_H__

//Quake3 fast 1/x² approximation
// For Magic Derivation see: Chris Lomont http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf
// Credited to Greg Walsh.
// 32  Bit float magic number - for 64 bits doubles: 0x5fe6ec85e7de30da
#define SQRT_MAGIC_D 0x5f3759df
#define SQRT_MAGIC_F 0x5f375a86

namespace mio {

namespace Optim {

	/**
	* @brief Optimized version of c++ round()
	* This version works with positive and negative numbers but does not
	* comply with IEEE handling of border cases (like infty, Nan, etc).
	* Please benchmark your code before deciding to use this!!
	* @param x number to round
	* @return rounded number cast as int
	*/
	inline long int round(const double& x) {
		if(x>=0.) return static_cast<long int>( x+.5 );
		else return static_cast<long int>( x-.5 );
	}

	/**
	* @brief Optimized version of c++ floor()
	* This version works with positive and negative numbers but does not
	* comply with IEEE handling of border cases (like infty, Nan, etc).
	* Please benchmark your code before deciding to use this!!
	* @param x number to floor
	* @return floored number cast as int
	*/
	inline long int floor(const double& x) {
		const long int xi = static_cast<long int>(x);
		if (x >= 0 || static_cast<double>(xi) == x) return xi ;
		else return xi - 1 ;
	}

	/**
	* @brief Optimized version of c++ ceil()
	* This version works with positive and negative numbers but does not
	* comply with IEEE handling of border cases (like infty, Nan, etc).
	* Please benchmark your code before deciding to use this!!
	* @param x number to ceil
	* @return ceiled number cast as int
	*/
	inline long int ceil(const double& x) {
		const long int xi = static_cast<long int>(x);
		if (x <= 0 || static_cast<double>(xi) == x) return xi ;
		else return xi + 1 ;
	}

	#ifdef _MSC_VER
	#pragma warning( push ) //for Visual C++
	#pragma warning(disable:4244) //Visual C++ righhtfully complains... but this behavior is what we want!
	#endif
	//maximum relative error is <1.7% while computation time for sqrt is <1/4. At 0, returns a large number
	//on a large scale interpolation test on TA, max relative error is 1e-6
	inline float invSqrt(const float x) {
		const float xhalf = 0.5f*x;

		union { // get bits for floating value
			float x;
			int i;
		} u;
		u.x = x;
		u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
		return u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
	}

	inline double invSqrt(const double x) {
		const double xhalf = 0.5f*x;

		union { // get bits for floating value
			float x;
			int i;
		} u;
		u.x = static_cast<float>(x);
		u.i = SQRT_MAGIC_D - (u.i >> 1);  // gives initial guess y0
		return u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
	}
	#ifdef _MSC_VER
	#pragma warning( pop ) //for Visual C++, restore previous warnings behavior
	#endif

	inline float fastSqrt_Q3(const float x) {
		return x * invSqrt(x);
	}

	inline double fastSqrt_Q3(const double x) {
		return x * invSqrt(x);
	}

	inline double pow2(const double& val) {return (val*val);}
	inline double pow3(const double& val) {return (val*val*val);}
	inline double pow4(const double& val) {return (val*val*val*val);}
}

} //end namespace

#endif
