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
#ifndef __LIBINTERPOL1D_H__
#define __LIBINTERPOL1D_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <cmath>
#include <vector>
#include <algorithm>

namespace mio {

class Interpol1D {
 	public:
		static double linearInterpolation(const double& x1, const double& y1, 
                                            const double& x2, const double& y2, const double& x3);
  		static double linearInterpolation(const double& d1, const double& d2, const double& weight=1.0);
		static double arithmeticMean(const std::vector<double>& vecData);
		static double getMedian(const std::vector<double>& vecData);
		static double getMedianAverageDeviation(const std::vector<double>& vecData);
		static double variance(const std::vector<double>& X);
		static double std_dev(const std::vector<double>& X);
		static double covariance(const std::vector<double>& z1, const std::vector<double>& z2);

		static void LinRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg);
		static int NoisyLinRegression(const std::vector<double>& in_X, const std::vector<double>& in_Y, double& A, double& B, double& R, std::stringstream& mesg);
		static void LogRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg);
		static void ExpRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg);
};
} //end namespace

#endif
