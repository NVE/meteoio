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
#include <meteoio/meteostats/libfit1D.h>
#include <meteoio/meteostats/libinterpol1D.h>
#include <cmath>
#include <algorithm>

using namespace std;

namespace mio {

Fit1D::Fit1D() : model(NULL)
{
}

//default constructor
Fit1D::Fit1D(const regression& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit) : model(NULL) {
	setModel(regType, in_X, in_Y, updatefit);
}

Fit1D::Fit1D(const std::string& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit) : model(NULL) {
	setModel(regType, in_X, in_Y, updatefit);
}

Fit1D::Fit1D(const Fit1D& i_fit) : model(NULL) { //HACK: teh pointer could not be valid anymore
	*this = i_fit;
}

Fit1D& Fit1D::operator=(const Fit1D& source) { //HACK: teh pointer could not be valid anymore
	if(this != &source) {
		model = new SimpleLinear; //this is only for memory allocation
		*model = *(source.model); //copy what is pointed to
	}
	return *this;
}

void Fit1D::setModel(const std::string& i_regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit) {
	regression regType;
	if(i_regType=="ZERO") regType=ZERO;
	else if(i_regType=="SIMPLE_LINEAR") regType=SIMPLE_LINEAR;
	else if(i_regType=="NOISYLINEAR") regType=NOISYLINEAR;
	else if(i_regType=="LINVARIO") regType=LINVARIO;
	else if(i_regType=="EXPVARIO") regType=EXPVARIO;
	else if(i_regType=="SPHERICVARIO") regType=SPHERICVARIO;
	else if(i_regType=="RATQUADVARIO") regType=RATQUADVARIO;
	else if(i_regType=="LINEARLS") regType=LINEARLS;
	else if(i_regType=="QUADRATIC") regType=QUADRATIC;
	else {
		throw IOException("The regression algorithm '"+i_regType+"' is not implemented" , AT);
	}

	setModel(regType, in_X, in_Y, updatefit);
}

void Fit1D::setModel(const regression& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y, const bool& updatefit) {
	if(model!=NULL) delete model;

	if(regType==ZERO) model=new Zero;
	if(regType==SIMPLE_LINEAR) model=new SimpleLinear;
	if(regType==NOISYLINEAR) model=new NoisyLinear;
	if(regType==LINVARIO) model=new LinVario;
	if(regType==EXPVARIO) model=new ExpVario;
	if(regType==SPHERICVARIO) model=new SphericVario;
	if(regType==RATQUADVARIO) model=new RatQuadVario;
	if(regType==LINEARLS) model=new LinearLS;
	if(regType==QUADRATIC) model=new Quadratic;

	model->setData(in_X, in_Y);
	if(updatefit) fit();
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

	nPts=X.size();
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

bool SimpleLinear::fit() {
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

bool NoisyLinear::fit() {
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
	//c0>=0, cs>=0, as>=0
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
	//c0>=0, b1>=0
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
	size_t xzero_idx=0;
	for(size_t i=1; i<X.size(); i++) {
		if(abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	const double slope = Interpol1D::arithmeticMean( Interpol1D::derivative(X, Y) );
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( slope );
}

double ExpVario::f(const double& x) {
	//c0>=0, ce>=0, ae>=0
	const double c0 = Lambda.at(0);
	const double ce = Lambda.at(1);
	const double ae = Lambda.at(2);

	if(x==0) {
		return 0;
	} else {
		const double y = c0 + ce * (1. - exp(-abs(x)/ae) );
		return y;
	}
}

void ExpVario::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for(size_t i=1; i<X.size(); i++) {
		if(abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( Y.back() - Y[xzero_idx] );
	Lambda.push_back( 1. );
}

double RatQuadVario::f(const double& x) {
	//c0>=0, cr>=0, ar>=0
	const double c0 = Lambda.at(0);
	const double cr = Lambda.at(1);
	const double ar = Lambda.at(2);

	if(x==0) {
		return 0;
	} else {
		const double y = c0 + cr*x*x / (1. + x*x/ar);
		return y;
	}
}

void RatQuadVario::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for(size_t i=1; i<X.size(); i++) {
		if(abs(X[i])<xzero) { xzero=X[i]; xzero_idx=i;}
	}
	Lambda.push_back( Y[xzero_idx] );
	Lambda.push_back( *( std::max_element( Y.begin(), Y.end() ) ) );
	Lambda.push_back( 1. );
}

double LinearLS::f(const double& x) {
	const double y = Lambda.at(0)*x + Lambda.at(1); //Lambda is a vector
	return y;
}

void LinearLS::setDefaultGuess() {
	double xzero=X[0];
	size_t xzero_idx=0;
	for(size_t i=1; i<X.size(); i++) {
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
	size_t xzero_idx=0;
	for(size_t i=1; i<der.size(); i++) {
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
