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
#ifndef SNOWPSUMINTERPOLATION_H
#define SNOWPSUMINTERPOLATION_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class SnowPSUMInterpolation
 * @brief Precipitation distribution according to the local slope and curvature.
 * The precipitation distribution is initialized using a specified algorithm (IDW_LAPSE by default, see IDWLapseAlgorithm).
 * An optional parameter can be given to specify which algorithm has to be used for initializing the grid.
 * Please do not forget to provide the arguments of the chosen algorithm itself if necessary!
 *
 * After this initialization, the pixels whose air temperatures are below or at freezing are modified according
 * to the method described in <i>"Quantitative evaluation of different hydrological modelling approaches
 * in a partly glacierized Swiss watershed"</i>, Magnusson et Al., Hydrological Processes, <b>25</b>, 2071-2084, 2011 and
 * <i>"Modelling runoff from highly glacierized alpine catchments in a changing climate"</i>, Huss et All., Hydrological Processes, <b>22</b>, 3888-3902, 2008.
 *
 * An example using this algorithm, initializing the grid with a constant lapse rate fill using +0.05% precipitation increase per meter of elevation, is given below:
 * @code
 * PSUM::algorithms = PSUM_SNOW
 * PSUM::psum_snow = avg_lapse
 * PSUM::avg_lapse = 0.0005 frac
 * @endcode
 *
 * @author Florian Kobierska, Jan Magnusson and Mathias Bavay
 */
class SnowPSUMInterpolation : public InterpolationAlgorithm {
	public:
		SnowPSUMInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
  			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

} //end namespace mio

#endif
