/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef ILWREPSALGORITHM_H
#define ILWREPSALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class ILWREpsAlgorithm
 * @ingroup spatialization
 * @brief Incoming Long Wave Radiation interpolation algorithm.
 * @details
 * Each ILWR is converted to an emissivity (using the local air temperature), interpolated using IDW_LAPSE and reconverted to ILWR. As
 * a side effect, the user must have defined algorithms to be used for air temperature (since this is needed for emissivity to ILWR conversion).
 * The lapse rate definition arguments as parsed by InterpolationAlgorithm::setTrendParams are supported as well as the following arguments:
 *  - SCALE: this is a scaling parameter to smooth the IDW distribution. In effect, this is added to the distance in order
 * to move into the tail of the 1/d distribution (default: 1000m);
 *  - ALPHA: this is an exponent to the 1/d distribution (default: 1);
 *
 * @code
 * ILWR::algorithms = ILWR_EPS
 * ILWR::ilwr_eps::soft = true
 * ILWR::ilwr_eps::rate = -0.03125
 * @endcode
 */
class ILWREpsAlgorithm : public InterpolationAlgorithm {
	public:
		ILWREpsAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
		                                 Meteo2DInterpolator& i_mi);
		virtual double getQualityRating(const Date& i_date);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		Trend trend;
		Meteo2DInterpolator& mi;
		std::vector<double> vecDataEA; ///<vectors of extracted emissivities
		double scale, alpha; ///<a scale parameter to smooth out the 1/dist and an exponent
};

} //end namespace mio

#endif
