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
#ifndef RYANALGORITHM_H
#define RYANALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class RyanAlgorithm
 * @ingroup spatialization
 * @brief DEM-based wind direction interpolation algorithm.
 * @details
 * This is an implementation of a method that alters a wind field (that is first computed with IDW)) based on the DEM, as described in Ryan,
 * <i>"a mathematical model for diagnosis and prediction of surface winds in mountainous terrain"</i>,
 * 1977, journal of applied meteorology, <b>16</b>, 6. It takes the following arguments:
 *  - SCALE: this is a scaling parameter to smooth the IDW distribution. In effect, this is added to the distance in order
 * to move into the tail of the 1/d distribution (default: 1000m);
 *  - ALPHA: this is an exponent to the 1/d distribution (default: 1);
 *
 * @code
 * DW::algorithms    = RYAN
 * @endcode
 */
class RyanAlgorithm : public InterpolationAlgorithm {
	public:
		RyanAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector< std::pair<std::string, std::string> >& vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataVW, vecDataDW; ///<vectors of extracted VW and DW
		double scale, alpha; ///<a scale parameter to smooth out the 1/dist and an exponent
		bool inputIsAllZeroes;
};

} //end namespace mio

#endif
