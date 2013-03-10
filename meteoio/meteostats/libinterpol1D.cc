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
#include <meteoio/meteostats/libinterpol1D.h>
#include <algorithm>
#include <cmath>

using namespace std;

namespace mio {

/**
 * @brief This function returns a vector of quantiles.
 * The vector does not have to be sorted. See https://secure.wikimedia.org/wikipedia/en/wiki/Quartile for more.
 * This code is heavily inspired by Ken Wilder, https://sites.google.com/site/jivsoft/Home/compute-ranks-of-elements-in-a-c---array-or-vector
 * (quantile method, replacing the nth-element call by direct access to a sorted vector).
 * @param X vector to classify
 * @param quartiles vector of quartiles, between 0 and 1
 * @return vector of ordinates of the quantiles
 */
std::vector<double> Interpol1D::quantiles(const std::vector<double>& X, const std::vector<double>& quartiles)
{
	const size_t Xsize = X.size();
	const size_t Qsize = quartiles.size();
	if (Xsize == 0)
		throw NoAvailableDataException("Trying to calculate quantiles with no data points", AT);
	if (Qsize == 0)
		throw NoAvailableDataException("No quantiles specified", AT);

	//in order to properly escape nodata points, we need to copy in a temporary vector
	vector<double> vecTemp;
	for(size_t i=0; i<Xsize; i++) {
		const double& value=X[i];
		if(value!=IOUtils::nodata)
			vecTemp.push_back(value);
	}
	std::sort( vecTemp.begin(), vecTemp.end()); //since we will process several values, we sort the vector

	//we will store results in a new vector
	std::vector<double> vecResults(Qsize, IOUtils::nodata);

	const size_t vecSize = vecTemp.size();
	if (vecSize == 0) {
		return vecResults; //ie: nodata values
	}
	if(vecSize == 1) {
		std::fill(vecResults.begin(), vecResults.end(), vecTemp[0]);
		return vecResults;
	}

	//compute quantiles
	for(size_t ii=0; ii<Qsize; ii++) {
		const double q = quartiles[ii];
		if(q<=0.) vecResults[ii] = vecTemp[0];
		else if(q>=1.) vecResults[ii] = vecTemp[vecSize-1];
		else {
			const double pos = static_cast<double>(vecSize - 1) * q;
			const unsigned int ind = static_cast<unsigned int>(pos);
			const double delta = pos - static_cast<double>(ind);
			const double i1 = vecTemp[ind];
			const double i2 = vecTemp[ind+1];
			vecResults[ii] = i1 * (1.0 - delta) + i2 * delta;
		}
	}
	return vecResults;
}

/**
 * @brief This function returns the vector of local derivatives, given a vector of abscissae and ordinates.
 * The vectors must be sorted by ascending x. The derivatives will be centered if possible, left or right otherwise or nodata
 * if nothing else can be computed.
 * @param X vector of abscissae
 * @param Y vector of ordinates
 * @return vector of local derivatives
 */
std::vector<double> Interpol1D::derivative(const std::vector<double>& X, const std::vector<double>& Y)
{
	if(X.size()!=Y.size()) {
		stringstream ss;
		ss << "X vector and Y vector don't match! " << X.size() << "!=" << Y.size() << "\n";
		throw InvalidArgumentException(ss.str(), AT);
	}

	std::vector<double> der;
	double right, centered, left, Dx_r, Dx_c, Dx_l;
	unsigned int i=0;
	//right hand derivative
	Dx_r=X[i+1]-X[i];
	if(Y[i]!=IOUtils::nodata && Y[i+1]!=IOUtils::nodata && Dx_r!=0.) der.push_back( (Y[i+1]-Y[i])/Dx_r );
		else der.push_back(IOUtils::nodata);

	//centered derivative if possible
	i++;
	for(; i<(X.size()-1); i++) {
		Dx_r=X[i-1]-X[i]; Dx_c=X[i+1]-X[i-1]; Dx_l=X[i]-X[i-1];

		if(Y[i]!=IOUtils::nodata && Y[i+1]!=IOUtils::nodata && Dx_r!=0.) right=(Y[i+1]-Y[i])/Dx_r; else right=IOUtils::nodata;
		if(Y[i]!=IOUtils::nodata && Y[i-1]!=IOUtils::nodata && Dx_l!=0.) left=(Y[i]-Y[i-1])/Dx_l; else left=IOUtils::nodata;
		if(Y[i-1]!=IOUtils::nodata && Y[i+1]!=IOUtils::nodata && Dx_c!=0.) centered=(Y[i+1]-Y[i-1])/Dx_c; else centered=IOUtils::nodata;

		if(centered!=IOUtils::nodata) der.push_back(centered);
		else if(right!=IOUtils::nodata) der.push_back(right);
		else if(left!=IOUtils::nodata) der.push_back(left);
		else der.push_back( IOUtils::nodata );
	}

	//left hand derivative
	Dx_l=X[i]-X[i-1];
	if(Y[i]!=IOUtils::nodata && Y[i-1]!=IOUtils::nodata && Dx_l!=0.) der.push_back( (Y[i]-Y[i-1])/Dx_l );
		else der.push_back(IOUtils::nodata);

	return der;
}

void Interpol1D::sort(std::vector<double>& X, std::vector<double>& Y)
{
	if(X.size()!=Y.size()) {
		stringstream ss;
		ss << "X vector and Y vector don't match! " << X.size() << "!=" << Y.size() << "\n";
		throw InvalidArgumentException(ss.str(), AT);
	}

	std::vector< std::pair<double,double> > new_vec( X.size() );
	for(unsigned int i=0; i<new_vec.size(); i++) {
		const std::pair<double,double> tmp(X[i],Y[i]);
		new_vec[i] = tmp;
	}

	std::sort( new_vec.begin(), new_vec.end(), pair_comparator );

	for(unsigned int i=0; i<new_vec.size(); i++) {
		X[i] = new_vec[i].first;
		Y[i] = new_vec[i].second;
	}
}

bool Interpol1D::pair_comparator(const std::pair<double, double>& l, const std::pair<double, double>& r) {
	return l.first < r.first;
}

/**
 * @brief This function returns the weighted aritmetic mean of two numbers.
 * A weight of 0 returns d1, a weight of 1 returns d2, a weight of 0.5 returns a centered mean.
 * See https://secure.wikimedia.org/wikipedia/en/wiki/Weighted_mean for more...
 * @param d1 first value
 * @param d2 second value
 * @param weight weight to apply to the mean
 * @return weighted aritmetic mean
 */
double Interpol1D::weightedMean(const double& d1, const double& d2, const double& weight)
{
	double tmp = abs(d1 - d2);

	if (d1 < d2) {
		return (d1 + tmp*weight);
	} else {
		return (d1 - tmp*weight);
	}
}

/**
 * @brief This function returns the weighted aritmetic mean of a vector.
 * See https://secure.wikimedia.org/wikipedia/en/wiki/Weighted_mean for more...
 * @param vecData vector of values
 * @param weight weights to apply to the mean
 * @return weighted aritmetic mean
 */
double Interpol1D::weightedMean(const std::vector<double>& vecData, const std::vector<double>& weight)
{
	const size_t nPts = vecData.size();
	if (nPts == 0)
		throw NoAvailableDataException("Trying to calculate an arithmetic mean with no data points", AT);
	if(nPts != weight.size()) {
		std::stringstream ss;
		ss << "Computing weighted mean of a vector of size " << nPts;
		ss << " with vector of weights of size " << weight.size();
		throw InvalidArgumentException(ss.str(), AT);
	}

	double sum = 0., count = 0.;
	for (size_t ii=0; ii<nPts; ii++){
		const double value=vecData[ii];
		if(value!=IOUtils::nodata) {
			const double w = weight[ii];
			sum += vecData[ii]*w;
			count += w;
		}
	}

	if(count>0.)
		return (sum/count);
	else
		return IOUtils::nodata;
}

double Interpol1D::arithmeticMean(const std::vector<double>& vecData)
{
	const size_t nPts = vecData.size();
	if (nPts == 0)
		throw NoAvailableDataException("Trying to calculate an arithmetic mean with no data points", AT);

	unsigned int count=0;
	double sum = 0.0;
	for (size_t ii=0; ii<nPts; ii++){
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
//This uses a sorting algorithm for getting middle element
// as much more efficient than full sorting (O(n) compared to O(n log(n))
	if (vecData.empty())
		throw NoAvailableDataException("Trying to calculate a median with no data points", AT);

	vector<double> vecTemp;
	for(size_t i=0; i<vecData.size(); i++) {
		const double& value=vecData[i];
		if(value!=IOUtils::nodata)
			vecTemp.push_back(value);
	}

	const size_t vecSize = vecTemp.size();
	if (vecSize == 0)
		return IOUtils::nodata;

	const int middle = (int)(vecSize/2);
	nth_element(vecTemp.begin(), vecTemp.begin()+middle, vecTemp.end());

	if ((vecSize % 2) == 1){ //uneven
		return *(vecTemp.begin()+middle);
	} else { //use arithmetic mean of element n/2 and n/2-1
		return weightedMean( *(vecTemp.begin()+middle), *(vecTemp.begin()+middle-1), 0.5);
	}
}

double Interpol1D::getMedianAverageDeviation(const std::vector<double>& vecData)
{
	if (vecData.empty())
		throw NoAvailableDataException("Trying to calculate MAD with no data points", AT);

	vector<double> vecWindow(vecData);

	const double median = Interpol1D::getMedian(vecWindow);
	if(median==IOUtils::nodata)
		return IOUtils::nodata;

	//Calculate vector of deviations and write each value back into the vecWindow
	for(size_t ii=0; ii<vecWindow.size(); ii++){
		double& value=vecWindow[ii];
		if(value!=IOUtils::nodata)
			value = std::abs(value - median);
	}

	//Calculate the median of the deviations
	const double mad = Interpol1D::getMedian(vecWindow);

	return mad;
}

double Interpol1D::variance(const std::vector<double>& X)
{//The variance is computed using a compensated variance algorithm,
//(see https://secure.wikimedia.org/wikipedia/en/wiki/Algorithms_for_calculating_variance)
//in order to be more robust to small variations around the mean.
	const size_t n = X.size();
	size_t count=0;
	double sum=0.;

	for(size_t i=0; i<n; i++) {
		const double value=X[i];
		if(value!=IOUtils::nodata) {
			sum += value;
			count++;
		}
	}

	if(count<=1) return IOUtils::nodata;

	const double mean = sum/(double)count;
	double sum2=0., sum3=0.;
	for(size_t i=0; i<n; i++) {
		const double value=X[i];
		if(value!=IOUtils::nodata) {
			sum2 = sum2 + (value - mean)*(value - mean);
			sum3 = sum3 + (value - mean);
		}
	}
	const double variance = (sum2 - sum3*sum3/static_cast<double>(count)) / static_cast<double>(count - 1);
	return variance;
}

double Interpol1D::std_dev(const std::vector<double>& X)
{
	return sqrt(variance(X));
}

double Interpol1D::covariance(const std::vector<double>& X, const std::vector<double>& Y)
{//this is a simple but still compensated covariance computation (see the notes on the variance)
	if(X.size()!=Y.size())
		throw IOException("Vectors should have the same size for covariance!", AT);
	const size_t n = X.size();
	if(n==0) return IOUtils::nodata;

	const double X_mean = Interpol1D::arithmeticMean(X);
	const double Y_mean = Interpol1D::arithmeticMean(Y);
	if(X_mean==IOUtils::nodata || Y_mean==IOUtils::nodata)
		return IOUtils::nodata;

	size_t count=0;
	double sum=0.;
	for(size_t i=0; i<n; i++) {
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
* @param r absolute value of linear regression coefficient
* @param mesg information message if something fishy is detected
*/
void Interpol1D::LinRegression(const std::vector<double>& X, const std::vector<double>& Y, double& a, double& b, double& r, std::stringstream& mesg)
{	//check arguments
	const size_t n=X.size();
	if(n==0)
		throw NoAvailableDataException("Trying to calculate linear regression with no data points", AT);
	if(n!=Y.size())
		throw IOException("Vectors should have the same size for linear regression!", AT);

	//computing x_avg and y_avg
	double x_avg=0., y_avg=0.;
	size_t count=0;
	for (size_t i=0; i<n; i++) {
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
	for (size_t i=0; i<n; i++) {
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
		r = fabs( sxy / sqrt(sx*sy) );
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
	const size_t nb_pts = in_X.size();
	double a,b,r;

	if (nb_pts==2) {
		mesg << "[W] only two points for linear regression!\n";
	}
	if (nb_pts<2) { //this should not be needed, we should have refrained from calling LinRegression in such a case
		mesg << "[E] Not enough data points for linear regression!\n";
		A=0.;
		if (!in_X.empty()) B=in_X[0]; else B=0.;
		R=1.;
		return EXIT_FAILURE;
	}

	LinRegression(in_X, in_Y, A, B, R, mesg);
	if(R>=r_thres)
		return EXIT_SUCCESS;

	std::vector<double> X(in_X), Y(in_Y);
	size_t nb_valid_pts = nb_pts;

	while(R<r_thres && nb_valid_pts>min_pts) {
		//we try to remove the one point in the data set that is the worst
		R=0.;
		size_t index_bad=0;
		for (size_t i=0; i<nb_pts; i++) {
			//invalidating alternatively each point
			const double Y_tmp=Y[i]; Y[i]=IOUtils::nodata;
			LinRegression(X, Y, a, b, r, mesg);
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
	if (R<r_thres) {
		mesg << "[W] Poor regression coefficient: " << std::setprecision(4) << R << "\n";
	}

	return EXIT_SUCCESS;
}

/**
* @brief Computes the bi-linear regression coefficients fitting the points given as X and Y in two vectors
* We consider that the regression can be made with 2 linear segments with a fixed inflection point. It relies on Interpol1D::NoisyLinRegression.
* @param in_X vector of X coordinates
* @param in_Y vector of Y coordinates (same order as X)
* @param bilin_inflection inflection point absissa
* @param coeffs a,b,r coefficients in a vector
* @return EXIT_SUCCESS or EXIT_FAILURE
*/
int Interpol1D::twoLinRegression(const std::vector<double>& in_X, const std::vector<double>& in_Y, const double& bilin_inflection, std::vector<double>& coeffs)
{
	//build segments
	std::vector<double> X1, Y1, X2, Y2;
	for(unsigned int ii=0; ii<in_X.size(); ii++) {
		if(in_X[ii]<bilin_inflection) { //first segment
			X1.push_back( in_X[ii] );
			Y1.push_back( in_Y.at(ii) );
		} else if(in_X[ii]>bilin_inflection) { //second segment
			X2.push_back( in_X[ii] );
			Y2.push_back( in_Y.at(ii) );
		} else { //point belongs to both segments
			X1.push_back( in_X[ii] );
			Y1.push_back( in_Y.at(ii) );
			X2.push_back( in_X[ii] );
			Y2.push_back( in_Y.at(ii) );
		}
	}

	double a1, b1, r1;
	double a2, b2, r2;
	//first segment
	std::stringstream mesg1;
	const int code1 = NoisyLinRegression(X1, Y1, a1, b1, r1, mesg1);

	//second segment
	std::stringstream mesg2;
	const int code2 = NoisyLinRegression(X2, Y2, a2, b2, r2, mesg2);

	if(code1==EXIT_FAILURE && code2==EXIT_FAILURE)
		return EXIT_FAILURE;

	coeffs.push_back(a1); coeffs.push_back(b1);
	coeffs.push_back(a2); coeffs.push_back(b2);
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
	std::vector<double> x(X.size());

	for(unsigned int i=0; i<X.size(); i++) {
		const double val = X[i];
		x[i] = (val!=IOUtils::nodata)? log(val) : IOUtils::nodata;
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
	std::vector<double> y( Y.size() );

	for(size_t i=0; i<Y.size(); i++) {
		const double val = Y[i];
		y[i] = (val!=IOUtils::nodata)? log(val) : IOUtils::nodata;
	}

	LinRegression(X, y, a, b, r, mesg); //HACK: how should we transform r?
}

} //namespace
