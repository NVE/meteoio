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
#ifndef __LIBFIT1DCORE_H__
#define __LIBFIT1DCORE_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Matrix.h>

#include <vector>

namespace mio {

//purely virtual class to use as an interface
class FitModel {
	public:
		FitModel() {};
		virtual void setData(const std::vector<double>& in_X, const std::vector<double>& in_Y) = 0;
		virtual void setGuess(const std::vector<double> lambda_in) = 0;
		virtual bool initFit() = 0;
		virtual double f(const double& x) = 0;
		void getParams(std::vector<double>& o_coefficients);
		std::string getModel() {return regname;};
		std::string getInfo();
	protected:
		std::string regname; //model name
		unsigned int nParam; //number of parameters
		unsigned int min_nb_pts; //minimum number of data points
		bool fit_ready;
		std::string infoString;
		std::vector<double> Lambda; //parameters of the fit
		std::vector<double> X; //X of input data set to fit
		std::vector<double> Y; //Y of input data set to fit
};

/**
 * @class FitLeastSquare
 * @brief A class to perform non-linear least square fitting.
 * It works on a time serie and uses matrix arithmetic to perform an arbitrary fit.
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-01-20
 */
class FitLeastSquare : public FitModel {
 	public:
		FitLeastSquare();
		void setData(const std::vector<double>& in_X, const std::vector<double>& in_Y);
		void setGuess(const std::vector<double> lambda_in);
		bool initFit();
		virtual double f(const double& x) = 0;

	protected:
		virtual void setDefaultGuess(); //set defaults guess values. Called by setData

	private:
		unsigned int nPts; //number of data points
		void checkInputs();
		void initLambda();
		void initDLambda(Matrix& dLambda) const;
		double getDelta(const double& var) const;
		double DDer(const double& x, const unsigned int& index);
		bool computeFit();

		static const double lambda_init; //initial default guess
		static const double delta_init_abs; //initial delta, absolute
		static const double delta_init_rel; //initial delta, relative
		static const double eps_conv; //convergence criteria
		static const unsigned int max_iter; //maximum number of iterations
};

} //end namespace

#endif
