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
#include <meteoio/libinterpol1D.h>

using namespace std;

namespace mio {

/**
 * @brief This function solves the equation y = kx + d for two given points and returns y for a given x
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

	unsigned int count=0;
	double sum = 0.0;
	for (unsigned int ii=0; ii<vecData.size(); ii++){
		const double value=vecData[ii];
		if(value!=IOUtils::nodata) {
			sum += vecData[ii];
			count++;
		}
	}

	if(count>0)
		return (sum/(double)count);
	else
		return IOUtils::nodata;
}

double Interpol1D::getMedian(const std::vector<double>& vecData)
{
	if (vecData.size() == 0)
		throw NoAvailableDataException("Trying to calculate a median with no data points", AT);
	
	vector<double> vecTemp;
	for(unsigned int i=0; i<vecData.size(); i++) {
		const double& value=vecData[i];
		if(value!=IOUtils::nodata)
			vecTemp.push_back(value);
	}

	if (vecTemp.size() == 0)
		return IOUtils::nodata;
	
	const unsigned int vecSize = vecTemp.size();
	const unsigned int middle = (unsigned int)(vecSize/2);
	sort(vecTemp.begin(), vecTemp.end());

	if ((vecSize % 2) == 1){ //uneven
		return vecTemp.at(middle);
	} else { //use arithmetic mean of element n/2 and n/2-1
		return Interpol1D::linearInterpolation(vecTemp.at(middle-1), vecTemp.at(middle), 0.5);
	}
}

double Interpol1D::getMedianAverageDeviation(const std::vector<double>& vecData)
{
	if (vecData.size() == 0)
		throw NoAvailableDataException("Trying to calculate MAD with no data points", AT);
	
	vector<double> vecWindow(vecData);

	double median = Interpol1D::getMedian(vecWindow);
	if(median==IOUtils::nodata)
		return IOUtils::nodata;

	//Calculate vector of deviations and write each value back into the vecWindow
	for(unsigned int ii=0; ii<vecWindow.size(); ii++){
		double& value=vecWindow[ii];
		if(value!=IOUtils::nodata)
			value = std::abs(value - median);
	}

	//Calculate the median of the deviations
	double mad = Interpol1D::getMedian(vecWindow);

	return mad;
}

double Interpol1D::variance(const std::vector<double>& X)
{
	const unsigned int n = X.size();

	unsigned int count=0;
	double sum=0., sum_sq=0.;

	for(unsigned int i=0; i<n; i++) {
		const double value=X[i];
		if(value!=IOUtils::nodata) {
			sum += value;
			sum_sq += value*value;
			count++;
		}
	}

	if(count<=1) return IOUtils::nodata;
	const double mean = sum/(double)count;
	return (sum_sq - sum*mean)/((double)count-1.);
}

double Interpol1D::std_dev(const std::vector<double>& X)
{
	return sqrt(variance(X));
}

double Interpol1D::covariance(const std::vector<double>& X, const std::vector<double>& Y)
{
	if(X.size()!=Y.size())
		throw IOException("Vectors should have the same size for covariance!", AT);
	const unsigned int n = X.size();
	if(n==0) return IOUtils::nodata;

	const double X_mean = Interpol1D::arithmeticMean(X);
	const double Y_mean = Interpol1D::arithmeticMean(Y);
	if(X_mean==IOUtils::nodata || Y_mean==IOUtils::nodata)
		return IOUtils::nodata;

	unsigned int count=0;
	double sum=0.;
	for(unsigned int i=0; i<n; i++) {
		if(X[i]!=IOUtils::nodata && Y[i]!=IOUtils::nodata) {
			sum += (X[i] - X_mean) * (Y[i] - Y_mean);
			count++;
		}
	}
	if(count<=1) return IOUtils::nodata;
	return sum/((double)count-1.);
}

/**
* @brief Computes the linear regression coefficients fitting the points given as X and Y in two vectors
* the linear regression has the form Y = aX + b with a regression coefficient r (it is nodata safe)
* @param X vector of X coordinates
* @param Y vector of Y coordinates (same order as X)
* @param a slope of the linear regression
* @param b origin of the linear regression
* @param r linear regression coefficient
* @param mesg information message if something fishy is detected
*/
void Interpol1D::LinRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg)
{	//check arguments
	const unsigned int n=X.size();
	if(n==0)
		throw NoAvailableDataException("Trying to calculate linear regression with no data points", AT);
	if(n!=Y.size())
		throw IOException("Vectors should have the same size for linear regression!", AT);

	//computing x_avg and y_avg
	double x_avg=0., y_avg=0.;
	unsigned int count=0;
	for (unsigned int i=0; i<n; i++) {
		if(X[i]!=IOUtils::nodata && Y[i]!=IOUtils::nodata) {
			x_avg += X[i];
			y_avg += Y[i];
			count++;
		}
	}
	if(count==0)
		throw NoAvailableDataException("Trying to calculate linear regression with no valid data points", AT);
	x_avg /= (double)count;
	y_avg /= (double)count;

	//computing sx, sy, sxy
	double sx=0., sy=0., sxy=0.;
	for (unsigned int i=0; i<n; i++) {
		if(X[i]!=IOUtils::nodata && Y[i]!=IOUtils::nodata) {
			sx += (X[i]-x_avg) * (X[i]-x_avg);
			sy += (Y[i]-y_avg) * (Y[i]-y_avg);
			sxy += (X[i]-x_avg) * (Y[i]-y_avg);
		}
	}

	//computing the regression line
	const double epsilon = 1e-6;
	if(sx <= abs(x_avg)*epsilon) { //sx and sy are always positive
		//all points have same X -> we return a constant value that is the average
		a = 0.;
		b = y_avg;
		r = 1.;
		mesg << "[W] Computing linear regression on data at identical X\n";
		return;
	}
	a = sxy / sx;
	b = y_avg - a*x_avg;
	if(sy==0) {
		//horizontal line: all y's are equals
		r = 1.;
	} else {
		//any other line
		r = sxy / sqrt(sx*sy);
	}
}

/**
* @brief Computes the linear regression coefficients fitting the points given as X and Y in two vectors
* the linear regression has the form Y = aX + b with a regression coefficient r. If the regression coefficient is not good enough, tries to remove bad points (up to 15% of the initial data set can be removed, keeping at least 4 points)
* @param in_X vector of X coordinates
* @param in_Y vector of Y coordinates (same order as X)
* @param A slope of the linear regression
* @param B origin of the linear regression
* @param R linear regression coefficient
* @param mesg information message if something fishy is detected
* @return EXIT_SUCCESS or EXIT_FAILURE
*/
int Interpol1D::NoisyLinRegression(const std::vector<double>& in_X, const std::vector<double>& in_Y, double& A, double& B, double& R, std::stringstream& mesg)
{
	//finds the linear regression for points (x,y,z,Value)
	const double r_thres=0.7;
	//we want at least 4 points AND 85% of the initial data set kept in the regression
	const unsigned int min_dataset=(unsigned int)floor(0.85*(double)in_X.size());
	const unsigned int min_pts=(min_dataset>4)?min_dataset:4;
	const unsigned int nb_pts = in_X.size();
	double a,b,r;

	if (nb_pts==2) {
		mesg << "[W] only two points for linear regression!\n";
	}
	if(nb_pts<2) { //this should not be needed, we should have refrained from calling LinRegression in such a case
		mesg << "[E] Not enough data point for linear regression!\n";
		A=0.;
		B=in_X[1];
		R=1.;
		return EXIT_FAILURE;
	}

	Interpol1D::LinRegression(in_X, in_Y, A, B, R, mesg);
	if(fabs(R)>=r_thres)
		return EXIT_SUCCESS;

	std::vector<double> X(in_X), Y(in_Y);
	unsigned int nb_valid_pts=nb_pts;

	while(fabs(R)<r_thres && nb_valid_pts>min_pts) {
		//we try to remove the one point in the data set that is the worst
		R=0.;
		unsigned int index_bad=0;
		for (unsigned int i=0; i<nb_pts; i++) {
			//invalidating alternatively each point
			const double Y_tmp=Y[i]; Y[i]=IOUtils::nodata;
			Interpol1D::LinRegression(X, Y, a, b, r, mesg);
			Y[i]=Y_tmp;

			if (fabs(r)>fabs(R)) {
				A=a;
				B=b;
				R=r;
				index_bad=i;
			}
		}
		//the worst point has been found, we overwrite it
		Y[index_bad]=IOUtils::nodata;
		nb_valid_pts--;
	}

	//check if r is reasonnable
	if (fabs(R)<r_thres) {
		mesg << "[W] Poor regression coefficient: " << std::setprecision(4) << R << "\n";
	}

	return EXIT_SUCCESS;
}

/**
* @brief Computes the Log regression coefficients fitting the points given as X and Y in two vectors
* the log regression has the form Y = a*ln(X) + b with a regression coefficient r (it is nodata safe)
* @param X vector of X coordinates
* @param Y vector of Y coordinates (same order as X)
* @param a slope of the regression
* @param b origin of the regression
* @param r regression coefficient
* @param mesg information message if something fishy is detected
*/
void Interpol1D::LogRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg)
{
	std::vector<double> x;

	for(unsigned int i=0; i<X.size(); i++) {
		const double val = X[i];
		if(val!=IOUtils::nodata) {
			x.push_back( log(val) );
		} else {
			x.push_back( IOUtils::nodata );
		}
	}

	LinRegression(x, Y, a, b, r, mesg); //HACK: how should we transform r?
}

/**
* @brief Computes the power regression coefficients fitting the points given as X and Y in two vectors
* the power regression has the form Y = b*X^a with a regression coefficient r (it is nodata safe)
* @param X vector of X coordinates
* @param Y vector of Y coordinates (same order as X)
* @param a slope of the regression
* @param b origin of the regression
* @param r regression coefficient
* @param mesg information message if something fishy is detected
*/
void Interpol1D::ExpRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg)
{
	std::vector<double> y;

	for(unsigned int i=0; i<Y.size(); i++) {
		const double val = Y[i];
		if(val!=IOUtils::nodata) {
			y.push_back( log(val) );
		} else {
			y.push_back( IOUtils::nodata );
		}
	}

	LinRegression(X, y, a, b, r, mesg); //HACK: how should we transform r?
}

} //namespace
