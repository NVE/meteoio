// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2021 MobyGIS Srl, Trento, Italy                                      */
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

#include <meteoio/gridResampling/GridResamplingAlgorithms.h>
#include <meteoio/gridResampling/GridLinearResampling.h>
#include <meteoio/gridResampling/GridTimeseriesResampling.h>

namespace mio {

/**
 * @brief Facade constructor for a generic grid resampling algorithm.
 */
GridResamplingAlgorithm::GridResamplingAlgorithm(const std::string& algorithm, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& /*vecArgs*/)
	: algo(algorithm), parname(i_parname), grid_window_size(dflt_window_size)
{
	//do nothing
}

/**
 * @brief Set this algorithm's window size to something other than the default value.
 * @param[in] window_size Desired window size in seconds.
 */
void GridResamplingAlgorithm::setWindowSize(const double& window_size)
{
	if (window_size <= 0.)
		throw InvalidArgumentException("Invalid WINDOW_SIZE for grid resampling algorithm", AT);
	grid_window_size = window_size / 86400.; //end user enters seconds, Julian internally
}

/**
 * @brief Object factory for temporal grid resampling algorithms.
 * @param[in] i_algorithm Semantic name of algorithm (as given in the INI file) to build.
 * @param[in] parname Meteo parameter to build the algorithm for.
 * @param[in] grid_window_size Standard window size for temporal grid resampling.
 * @param[in] vecArgs The algorithm's parameters as parsed from the user setings.
 * @return A resampling algorithm object for the desired parameter.
 */
GridResamplingAlgorithm* GridResamplingAlgorithmsFactory::getAlgorithm(const std::string& i_algorithm, const std::string& parname,
	const double& grid_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string algorithm( IOUtils::strToUpper(i_algorithm) );
	if (algorithm == "LINEAR")
		return new GridLinearResampling(algorithm, parname, grid_window_size, vecArgs);
	else if (algorithm == "TIMESERIES")
		return new GridTimeseriesResampling(algorithm, parname, grid_window_size, vecArgs);
	else
		throw IOException("The temporal grid resampling algorithm '" + algorithm + "' is not implemented", AT);

	return nullptr;
}

} //namespace

