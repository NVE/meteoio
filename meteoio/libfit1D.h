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
#ifndef __LIBFIT1D_H__
#define __LIBFIT1D_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Matrix.h>

//#include <cmath>
#include <vector>

namespace mio {

class Fit1D;

typedef double (Fit1D::*FitFctPtr)(const double& X);

/**
 * @class Fit1D
 * @brief A class to perform non-linear least square fitting.
 * It works on a time serie and uses matrix arithmetic to perform an arbitrary fit.
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-01-20
 */
class Fit1D {
 	public:
		Fit1D(const std::vector<double>& in_X, const std::vector<double>& in_Y);

		void setGuess(const std::vector<double> lambda_in);
		void leastSquareFit();

	private:
		static const double lambda_init; //initial default guess
		static const double delta_init_abs; //initial delta, absolute
		static const double delta_init_rel; //initial delta, relative
		static const double eps_conv; //convergence criteria
		FitFctPtr fitFct; //fit function pointer
		unsigned int nPts, nParam; //number of data points, number of parameters
		const std::vector<double>& X, Y; //X & Y of input data set to fit

		std::vector<double> Lambda; //parameters of the fit
		Matrix dLambda; //parameters variations;

		void initLambda();
		void initDLambda();
		double DDer(const double& X, const unsigned int& index);

		//various fit functions
		double LinFit(const double& X);
};

} //end namespace

#endif
