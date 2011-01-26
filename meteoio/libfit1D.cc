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

//default constructor
Fit1D::Fit1D(const std::vector<double>& in_X, const std::vector<double>& in_Y) : X(in_X), Y(in_Y) {
//TODO: set the interpolation type, establish nbPts, nbParam matching with fit type
	fitFct = &Fit1D::LinFit;
	nPts=3;
	nParam=2;
	dLambda.resize(nParam,1);
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
}

void Fit1D::initLambda() {
	if(Lambda.size()==0) //else, setGuess has been called
		Lambda.resize(nParam, lambda_init);
}

void Fit1D::initDLambda() {
	for(unsigned int m=1; m<=nParam; m++) {
		const double var = Lambda[m];
		if(var==0) {
			dLambda(m,1) = delta_init_abs;
		} else {
			dLambda(m,1) = delta_init_rel * var;
		}
	}
}

double Fit1D::DDer(const double& X, const unsigned int& index) {
	const double var = Lambda.at(index);
	const double delta = dLambda(index,1)/2.;
	const double v1 = var - delta;
	const double v2 = var + delta;

	Lambda[index] = v1;
	const double Y1 = (this->*(fitFct))(X);
	Lambda[index] = v2;
	const double Y2 = (this->*(fitFct))(X);
	Lambda[index] = var;

	return (Y2-Y1)/(v2-v1);
}

void Fit1D::leastSquareFit() {
	double max_delta = std::numeric_limits<double>::max();
	initLambda();
	initDLambda();

	do {
		Matrix A(nPts, nParam);
		for(unsigned int m=1; m<=nPts; m++) {
			for(unsigned int n=1; n<=nParam; n++) {
				A(m,n) = DDer( X[m], n );
			}
		}

		//calculate parameters deltas
		const Matrix a = A.getT() * A;
		const Matrix dBeta = A*dLambda;
		const Matrix b = A.getT() * dBeta;
		dLambda = Matrix::solve(a,b);

		//apply the deltas to the parameters, record maximum delta
		max_delta = 0.;
		for(unsigned int m=1; m<=nParam; m++) {
			Lambda[m] += dLambda(m,1);
			if( fabs(dLambda(m,1))>max_delta ) max_delta=fabs(dLambda(m,1));
		}

		//compute R2
		//const double R2 = Matrix::dot(dBeta, dBeta);
	} while (max_delta>eps_conv);

	std::cout << "Coefficients:\n";
	for(unsigned int i=0; i<Lambda.size(); i++) {
		std::cout << Lambda[i] << "\n";
	}
	
}

double Fit1D::LinFit(const double& X) {
	const double Y = Lambda.at(0)*X + Lambda.at(1);
	return Y;
}

} //namespace
