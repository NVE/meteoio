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
#include "libinterpol1D.h"

using namespace std;
using namespace mio;

double Interpol1D::linearInterpolation(const double& d1, const double& d2, const double& weight)
{
	double tmp = abs(d1 - d2);
	
	if (d1 < d2) {
		return (d1 + tmp*weight);
	} else {
		return (d1 - tmp*weight); 
	}
	return 0;
}


double Interpol1D::arithmeticMean(const std::vector<double>& vecData)
{
	if (vecData.size() == 0)
		throw NoAvailableDataException("Trying to calculate an arithmetic mean with no data points", AT);

	double sum = 0.0;
	for (unsigned int ii=0; ii<vecData.size(); ii++){
		sum += vecData[ii];
	}

	return (sum/(double)vecData.size());
}
