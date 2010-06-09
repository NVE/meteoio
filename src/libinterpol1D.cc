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

namespace mio {

/**
 * @brief This function solves the equation y = kx + d for two given points and return the y for a given x
 * @param x1 x-coordinate of first point
 * @param y1 y-coordinate of first point
 * @param x2 x-coordinate of second point
 * @param y2 y-coordinate of second point
 * @param x3 x-coordinate of desired point
 * @return y-coordinate of desired point
 */
double Interpol1D::linearInterpolation(const double& x1, const double& y1, 
                                       const double& x2, const double& y2, const double& x3)
{
	if (x1 == x2)
		throw IOException("Attempted division by null", AT);

	//Solving y = kx +d
	double k = (y1 - y2) / (x1 - x2);
	double d = y2 - k*x2;

	return (k*x3 + d);
}

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

double Interpol1D::getMedian(const std::vector<double>& vecData)
{
	if (vecData.size() == 0)
		throw NoAvailableDataException("Trying to calculate a median with no data points", AT);

	vector<double> vecTemp(vecData); //copy by value of vecData
	
	double median = 0.0;
	unsigned int vecSize = vecTemp.size();
	unsigned int middle = (unsigned int)(vecSize/2);
	sort(vecTemp.begin(), vecTemp.end());

	if ((vecSize % 2) == 1){ //uneven
		median = vecTemp.at(middle);
	} else { //use arithmetic mean of element n/2 and n/2-1
		median = Interpol1D::linearInterpolation(vecTemp.at(middle-1), vecTemp.at(middle), 0.5);
	}

	return median;
}

double Interpol1D::getMedianAverageDeviation(const std::vector<double>& vecData)
{
	if (vecData.size() == 0)
		throw NoAvailableDataException("Trying to calculate MAD with no data points", AT);
	
	vector<double> vecWindow(vecData);

	double median = Interpol1D::getMedian(vecData);

	//Calculate vector of deviations and write each value back into the vecWindow
	for(unsigned int ii=0; ii<vecWindow.size(); ii++){
		vecWindow[ii] = std::abs(vecWindow[ii] - median);
	}

	//Calculate the median of the deviations
	double mad = Interpol1D::getMedian(vecWindow);

	return mad;
}

} //namespace
