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
#include <meteoio/meteostats/libfit1DCore.h>
#include <meteoio/meteostats/libinterpol1D.h>
#include <cmath>

using namespace std;

namespace mio {

void FitModel::getParams(std::vector<double>& o_coefficients) {
	if(fit_ready!=true) {
		throw InvalidArgumentException("The regression has not yet being computed!", AT);
	}
	o_coefficients = Lambda;
}

std::string FitModel::getInfo() {
	if(fit_ready==true) {
		return infoString;
	} else {
		throw InvalidArgumentException("The regression has not yet being computed!", AT);
	}
}

const double FitLeastSquare::lambda_init = 1.; //initial default guess
const double FitLeastSquare::delta_init_abs = 1.; //initial delta, absolute
const double FitLeastSquare::delta_init_rel = 0.2; //initial delta, relative
const double FitLeastSquare::eps_conv = 1e-6; //convergence criteria
const unsigned int FitLeastSquare::max_iter = 50; //maximum number of iterations

//default constructor
FitLeastSquare::FitLeastSquare() {
	fit_ready = false;
	nParam = 0;
	min_nb_pts = 0;
	regname = "none";
}

void FitLeastSquare::setData(const std::vector<double>& in_X, const std::vector<double>& in_Y) {
	X = in_X;
	Y = in_Y;
	checkInputs();
	//sort the data by increasing X
	Interpol1D::sort(X, Y);
	setDefaultGuess();
	fit_ready = false;
}

void FitLeastSquare::setDefaultGuess() {
	for(size_t i=0; i<nParam; i++) {
		Lambda.push_back( lambda_init );
	}
}

void FitLeastSquare::setGuess(const std::vector<double> lambda_in) {
	const size_t nGuess = lambda_in.size();

	if(nGuess!=nParam) {
		stringstream ss;
		ss << "Provided " << nGuess << " guesses for " << nParam << " parameters!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for(size_t i=0; i<nGuess; i++) {
		Lambda.push_back( lambda_in[i] );
	}
	fit_ready = false;
}

bool FitLeastSquare::initFit() {
	return computeFit();
}

////////////////////////////////////////////////////////////
//// End of public methods

void FitLeastSquare::checkInputs() {
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
}

bool FitLeastSquare::computeFit() {
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
			dBeta(m,1) = Y[m-1] - f(X[m-1]); //X and Y are vectors
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
	ss << "Computed regression with " << regname << " model ";
	ss << "- Sum of square residuals = " << R2 << " , max_delta = " << max_delta << " ";
	ss << "- " << iter << " iterations";
	infoString = ss.str();

	if(max_delta>eps_conv) { //it did not converge
		fit_ready = false;
		return false;
	} else {		 //it did converge
		fit_ready = true;
		return true;
	}
}

void FitLeastSquare::initLambda() {
	if(Lambda.size()==0) //else, setGuess has been called
		Lambda.resize(nParam, lambda_init);
}

void FitLeastSquare::initDLambda(Matrix& dLambda) const {
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

double FitLeastSquare::getDelta(const double& var) const {
//calculate a sensible delta for the partial derivative
	if(var==0) {
		return (delta_init_abs * 0.5);
	} else {
		return (delta_init_rel * var * 0.5);
	}
}

double FitLeastSquare::DDer(const double& x, const unsigned int& index) {
	const double var = Lambda[index-1]; //Lambda is a vector
	const double delta = getDelta(var);
	const double v1 = var - delta;
	const double v2 = var + delta;

	Lambda[index-1] = v1;
	const double y1 = f(x);
	Lambda[index-1] = v2;
	const double y2 = f(x);
	Lambda[index-1] = var;

	return (y2-y1)/(v2-v1);
}

} //namespace
