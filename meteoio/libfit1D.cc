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
//#include <algorithm>

using namespace std;

namespace mio {

const double Fit1D::lambda_init = 1.; //initial default guess
const double Fit1D::delta_init_abs = 1.; //initial delta, absolute
const double Fit1D::delta_init_rel = 0.2; //initial delta, relative
const double Fit1D::eps_conv = 1e-3; //convergence criteria
const unsigned int Fit1D::max_iter = 20; //maximum number of iterations

//static section
//HACK: do property map: reg_name, fct_ptr, nb_parameters, min_nb_pts
/*void Fit1D::registerInterpolations() {
	mapRegs["name"]  = Regmodel(fct_ptr, nb_param, min_nb_pts);
}*/

//default constructor
Fit1D::Fit1D(const std::vector<double>& in_X, const std::vector<double>& in_Y) : X(in_X), Y(in_Y) {
//TODO: set the interpolation type, establish nbPts, nbParam matching with fit type
	//fitFct = &Fit1D::LinFit;
	fitFct = &Fit1D::SqFit;
	if( X.size()!=Y.size() ) {
		stringstream ss;
		ss << "X vector and Y vector don't match! " << X.size() << "!=" << Y.size() << "\n";
		throw InvalidArgumentException(ss.str(), AT);
	}
	nPts=X.size();
	nParam=3;
	fit_ready = false;
}

/**
* @brief Calculate a value using the computed least square fit. 
* The fit has to be computed before.
* @param x abscissa
* @return f(x) using the computed least square fit
*/
double Fit1D::getF(const double& x) {
	if(fit_ready==true) {
		return (this->*(fitFct))(x);
	} else {
		throw InvalidArgumentException("The regression has not yet being computed!", AT);
	}
}

/**
* @brief Return a string of information about the fit. 
* The fit has to be computed before.
* @return info string
*/
std::string Fit1D::getInfo() {
	if(fit_ready==true) {
		return infoString;
	} else {
		throw InvalidArgumentException("The regression has not yet being computed!", AT);
	}
}

void Fit1D::setGuess(const std::vector<double> lambda_in) {
	const unsigned int nGuess = lambda_in.size();

	if(nGuess!=nParam) {
		stringstream ss;
		ss << "Provided " << nGuess << " guesses for " << nParam << " parameters!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for(unsigned int i=0; i<=nGuess; i++) {
		Lambda.push_back( lambda_in[i] );
	}
	fit_ready = false;
}

void Fit1D::initLambda() {
	if(Lambda.size()==0) //else, setGuess has been called
		Lambda.resize(nParam, lambda_init);
}

void Fit1D::initDLambda(Matrix& dLambda) const {
	dLambda.resize(nParam,1);
	for(unsigned int m=1; m<=nParam; m++) {
		const double var = Lambda[m-1]; //Lambda is a vector
		if(var==0) {
			dLambda(m,1) = delta_init_abs * 0.5;
		} else {
			dLambda(m,1) = delta_init_rel * var * 0.5;
		}
	}
}

double Fit1D::getDelta(const double& var) const {
//calculate a sensible delta for the partial derivative
	if(var==0) {
		return (delta_init_abs * 0.5);
	} else {
		return (delta_init_rel * var * 0.5);
	}
}

double Fit1D::DDer(const double& x, const unsigned int& index) {
	const double var = Lambda[index-1]; //Lambda is a vector
	const double delta = getDelta(var);
	const double v1 = var - delta;
	const double v2 = var + delta;

	Lambda[index-1] = v1;
	const double y1 = (this->*(fitFct))(x);
	Lambda[index-1] = v2;
	const double y2 = (this->*(fitFct))(x);
	Lambda[index-1] = var;

	return (y2-y1)/(v2-v1);
}

bool Fit1D::leastSquareFit(std::vector<double>& coefficients) {
	double max_delta = std::numeric_limits<double>::max();
	initLambda();
	Matrix dLambda; //parameters variations
	initDLambda(dLambda);

	Matrix A(nPts, nParam);
	Matrix dBeta(nPts,(unsigned)1);
	
	unsigned int iter = 0;
	do {
		iter++;
		//set A matrix
		for(unsigned int m=1; m<=nPts; m++) {
			for(unsigned int n=1; n<=nParam; n++) {
				const double value = DDer( X[m-1], n ); //X is a vector
				A(m,n) = value;
			}
		}

		//set dBeta matrix
		for(unsigned int m=1; m<=nPts; m++) {
			dBeta(m,1) = Y[m-1] - (this->*(fitFct))(X[m-1]); //X and Y are vectors
		}

		//calculate parameters deltas
		const Matrix a = A.getT() * A;
		const Matrix b = A.getT() * dBeta;
		Matrix::solve(a, b, dLambda);

		//apply the deltas to the parameters, record maximum delta
		max_delta = 0.;
		for(unsigned int m=1; m<=nParam; m++) {
			Lambda[m-1] += dLambda(m,1); //Lambda is a vector
			if( fabs(dLambda(m,1))>max_delta ) max_delta=fabs(dLambda(m,1));
		}

	} while (max_delta>eps_conv && iter<max_iter);

	//compute R2
	const double R2 = Matrix::dot(dBeta, dBeta);

	//building infoString
	stringstream ss;
	ss << "Computed regression with *** model - Sum of square residuals = " << R2 << " , max_delta = " << max_delta;
	infoString = ss.str();

	coefficients = Lambda;

	if(max_delta>eps_conv) { //it did not converge
		fit_ready = false;
		return false;
	} else {		 //it did converge
		fit_ready = true;
		return true;
	}
}

bool Fit1D::leastSquareFit() {
	vector<double> coefficients;
	return leastSquareFit(coefficients);
}

double Fit1D::LinFit(const double& x) {
	const double y = Lambda.at(0)*x + Lambda.at(1); //Lambda is a vector
	return y;
}

double Fit1D::SqFit(const double& x) {
	const double y = Lambda.at(0)*x*x + Lambda.at(1)*x + Lambda.at(2); //Lambda is a vector
	return y;
}

} //namespace
