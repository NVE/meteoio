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
#ifndef __LIBFIT1D_H__
#define __LIBFIT1D_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/meteostats/libfit1DCore.h>

#include <vector>

namespace mio {

class SimpleLinear : public FitModel {
	public:
		SimpleLinear() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "SimpleLinear";};
		void setData(const std::vector<double>& in_X, const std::vector<double>& in_Y);
		void setGuess(const std::vector<double> /*lambda_in*/) {};
		bool initFit();
		double f(const double& x);
};

class NoisyLinear : public SimpleLinear {
	public:
		NoisyLinear() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "NoisyLinear";};
		bool initFit();
};

class SphericVario : public FitLeastSquare {
	public:
		SphericVario() {fit_ready = false; nParam = 3; min_nb_pts = 4; regname = "SphericVario";};
		void setDefaultGuess();
		double f(const double& x);
};

class LinVario : public FitLeastSquare {
	public:
		LinVario() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "LinVario";};
		void setDefaultGuess();
		double f(const double& x);
};

class LinearLS : public FitLeastSquare {
	public:
		LinearLS() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "LinearLS";};
		void setDefaultGuess();
		double f(const double& x);
};

class Quadratic : public FitLeastSquare {
	public:
		Quadratic() {fit_ready = false; nParam = 2; min_nb_pts = 3; regname = "Quadratic";};
		void setDefaultGuess();
		double f(const double& x);
};


/**
 * @class Fit1D
 * @brief A class to perform 1D regressions
 * It works on a time serie and uses matrix arithmetic to perform an arbitrary fit.
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-01-20
 */
class Fit1D {
 	public:
		/**
		* @brief Constructor.
		* @param regType regression model to use
		* @param in_X vector of data points abscissae
		* @param in_Y vector of data points ordinates
		*/
		Fit1D(const std::string& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y);

		~Fit1D();

		/**
		* @brief Provide a set of initial values for the model parameters.
		* @param lambda_in one initial value per model parameter
		*/
		void setGuess(const std::vector<double> lambda_in) {fit->setGuess(lambda_in);};

		/**
		* @brief Compute the regression parameters
		* @return false if could not find the parameters
		*/
		bool initFit() {return fit->initFit();};

		/**
		* @brief Calculate a value using the computed least square fit.
		* The fit has to be computed before.
		* @param x abscissa
		* @return f(x) using the computed least square fit
		*/
		double f(const double& x) {return fit->f(x);};

		/**
		* @brief Calculate the parameters of the fit.
		* The fit has to be computed before.
		* @param coefficients vector containing the coefficients
		*/
		void getParams(std::vector<double>& coefficients) {fit->getParams(coefficients);};

		/**
		* @brief Return the name of the fit model.
		* @return model name
		*/
		std::string getModel() {return fit->getModel();};

		/**
		* @brief Return a string of information about the fit.
		* The fit has to be computed before.
		* @return info string
		*/
		std::string getInfo() {return fit->getInfo();};

	private:
		FitModel *fit;
};

} //end namespace

#endif
