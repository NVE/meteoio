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
#ifndef IDWALGORITHM_H
#define IDWALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class IDWAlgorithm
 * @ingroup spatialization
 * @brief Inverse Distance Weighting interpolation algorithm.
 * @details
 * Each cell receives the weighted average of the whole data set with weights being 1/r²
 * (r being the distance of the current cell to the contributing station) and renormalized
 * (so that the sum of the weights is equal to 1.0).
 */
class IDWAlgorithm : public InterpolationAlgorithm {
	public:
		IDWAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector< std::pair<std::string, std::string> >& vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

} //end namespace mio

#endif
