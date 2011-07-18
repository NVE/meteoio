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
#include <meteoio/libfit1D.h>
#include <meteoio/libinterpol1D.h>
#include <cmath>
#include <algorithm>

using namespace std;

namespace mio {

//default constructor
Fit1D::Fit1D(const std::string& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y) {
	if(regType=="SimpleLinear") fit=new SimpleLinear;
	if(regType=="NoisyLinear") fit=new NoisyLinear;
	if(regType=="SphericVario") fit=new SphericVario;
	if(regType=="LinVario") fit=new LinVario;
	if(regType=="LinearLS") fit=new LinearLS;
	if(regType=="Quadratic") fit=new Quadratic;

	fit->setData(in_X, in_Y);
}

Fit1D::~Fit1D() {
	delete fit;
}

//////////////////////////////////////////////////////////////
// regression models
void SimpleLinear::setData(const std::vector<double>& in_X, const std::vector<double>& in_Y) {
	X = in_X;
	Y = in_Y;

	//check input data consistency
	if( X.size()!=Y.size() ) {
		stringstream ss;
		ss << "X vector and Y vector don't match! " << X.size() << "!=" << Y.size() << "\n";
		throw InvalidArgumentException(ss.str(), AT);
	}

	const double nPts=X.size();
	if(nPts<min_nb_pts) {
		stringstream ss;
		ss << "Only " << nPts << " data points for " << regname << " regression model.";
		ss << " Expecting at least " << min_nb_pts << " for this model!\n";
		throw InvalidArgumentException(ss.str(), AT);
	}

	fit_ready = false;
}

double SimpleLinear::f(const double& x) {
	return Lambda.at(0)*x + Lambda.at(1);
}

bool SimpleLinear::initFit() {
	Lambda.clear();
	double a,b,r;
	std::stringstream mesg;
	Interpol1D::LinRegression(X, Y, a, b, r, mesg);
	Lambda.push_back(a);
	Lambda.push_back(b);
	mesg << "Computed regression with " << regname << " model - r=" << r;
	infoString = mesg.str();
	fit_ready = true;
	return true;
}

bool NoisyLinear::initFit() {
	Lambda.clear();
	double a,b,r;
	std::stringstream mesg;
	Interpol1D::NoisyLinRegression(X, Y, a, b, r, mesg);
	Lambda.push_back(a);
	Lambda.push_back(b);
	mesg << "Computed regression with " << regname << " model - r=" << r;
	infoString = mesg.str();
	fit_ready = true;
	return true;
}

//regression models using the standard least square algorithm
double SphericVario::f(const double& x) {
	const double c0 = Lambda.at(0);
	const double cs = Lambda.at(1);
	const double as = Lambda.at(2);

	if(x==0) return 0;

	const double abs_x = fabs(x);
	if(abs_x>0 && abs_x<=as) {
		const double val = abs_x/as;
		const double y = c0 + cs * ( 1.5*val - 0.5*val*val*val );
		return y;
	} else {
		return (c0+cs);
	}
}

void SphericVario::setDefaultGuess() {
	Lambda.push_back( *min_element(Y.begin(), Y.end()) );
	Lambda.push_back( *max_element(Y.begin(), Y.end()) );
	Lambda.push_back( *max_element(X.begin(), X.end()) );
}

double LinVario::f(const double& x) {
	const double c0 = Lambda.at(0);
	const double bl = Lambda.at(1);

	if(x==0) {
		return 0;
	} else {
		const double y = c0 + bl * abs(x);
		return y;
	}
}

void LinVario::setDefaultGuess() {
	double xzero=X[0];
	unsigned int xzero_idx=0;
	for(unsigned int i=1; i<X.size(); i++) {
		if(abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	const double slope = Interpol1D::arithmeticMean( Interpol1D::derivative(X, Y) );
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( slope );
}

double LinearLS::f(const double& x) {
	const double y = Lambda.at(0)*x + Lambda.at(1); //Lambda is a vector
	return y;
}

void LinearLS::setDefaultGuess() {
	double xzero=X[0];
	unsigned int xzero_idx=0;
	for(unsigned int i=1; i<X.size(); i++) {
		if(abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}

	const double slope = Interpol1D::arithmeticMean( Interpol1D::derivative(X, Y) );
	Lambda.push_back( slope );
	Lambda.push_back( Y[xzero_idx] );
}

double Quadratic::f(const double& x) {
	const double y = Lambda.at(0)*x*x + Lambda.at(1)*x + Lambda.at(2); //Lambda is a vector
	return y;
}

void Quadratic::setDefaultGuess() {
	std::vector<double> der = Interpol1D::derivative(X, Y);
	const double acc = 0.5 * Interpol1D::arithmeticMean( Interpol1D::derivative(X, der) );
	double xzero=der[0];
	unsigned int xzero_idx=0;
	for(unsigned int i=1; i<der.size(); i++) {
		if(abs(der[i])<xzero) { xzero=der[i]; xzero_idx=i;}
	}

	Lambda.push_back( acc ); //0
	Lambda.push_back( der[xzero_idx] ); //1
	if(acc>0.)
		Lambda.push_back( *( std::min_element( Y.begin(), Y.end() ) ) );
	else
		Lambda.push_back( *( std::max_element( Y.begin(), Y.end() ) ) );
}

} //namespace
