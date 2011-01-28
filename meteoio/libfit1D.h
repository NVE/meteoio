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
#include <meteoio/IOUtils.h>
#include <meteoio/Matrix.h>

//#include <cmath>
#include <vector>

namespace mio {

class Fit1D;

typedef double (Fit1D::*FitFctPtr)(const double& X);

class RegModel {
	public:
		std::string name; ///< A string representing the regression model
		FitFctPtr fitFct; ///< The pointer to the regression model
		unsigned int nb_param; ///< Number of parameters of the model
		unsigned int min_nb_pts; ///< Minimum number of data points required to use the model
		
		/**
		 * @brief The main constructor for the RegModel class
		 *
		 * @param _s1 A std::string representing the file to be opened (or "" if plugin is statically linked)
		 * @param _s2 A std::string that is the classname of the object to be loaded (e.g. "A3DIO", "GSNIO")
		 * @param p1  A pointer to the loaded object of type IOInterface (or NULL)
		 * @param p2  A pointer to the loaded dynamic library (or NULL)
		 */
		RegModel(const std::string in_name, FitFctPtr in_fitFct, const unsigned int& in_nb_param, const unsigned int& in_min_nb_pts) : name(in_name), fitFct(in_fitFct), nb_param(in_nb_param), min_nb_pts(in_min_nb_pts){}
		RegModel() : name(""), fitFct(NULL), nb_param(0), min_nb_pts(0){}

		friend std::ostream& operator<<(std::ostream& os, const RegModel& data);
		static const std::string header; //to contain a helpful header for understanding the output of <<

};

/**
 * @class Fit1D
 * @brief A class to perform non-linear least square fitting.
 * It works on a time serie and uses matrix arithmetic to perform an arbitrary fit.
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-01-20
 */
class Fit1D {
 	public:
		Fit1D(const std::string& regType, const std::vector<double>& in_X, const std::vector<double>& in_Y);

		void setGuess(const std::vector<double> lambda_in);
		bool leastSquareFit(std::vector<double>& coefficients);
		bool leastSquareFit();
		double getF(const double& x);
		std::string getInfo();

	private:
		static const double lambda_init; //initial default guess
		static const double delta_init_abs; //initial delta, absolute
		static const double delta_init_rel; //initial delta, relative
		static const double eps_conv; //convergence criteria
		static const unsigned int max_iter; //maximum number of iterations

	private:
		FitFctPtr fitFct; //fit function pointer
		unsigned int nPts, nParam; //number of data points, number of parameters
		std::string regname; //human readable regression model name
		bool fit_ready;
		std::string infoString;
		const std::vector<double>& X; //X of input data set to fit
		const std::vector<double>& Y; //Y of input data set to fit

		std::vector<double> Lambda; //parameters of the fit
		std::map<std::string, RegModel::RegModel> mapRegs;

		void registerRegressions();
		void initLambda();
		void initDLambda(Matrix& dLambda) const;
		double getDelta(const double& var) const;
		double DDer(const double& x, const unsigned int& index);

		//various fit models
		double LinFit(const double& x);
		double SqFit(const double& x);

		//variogram fit models
		double LinVario(const double& x);
		double SphericVario(const double& x);
};

} //end namespace

#endif
