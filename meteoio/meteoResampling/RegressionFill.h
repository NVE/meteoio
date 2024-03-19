// SPDX-License-Identifier: LGPL-3.0-or-later
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
#ifndef RESAMPLINGREGRESSIONFILL_H
#define RESAMPLINGREGRESSIONFILL_H

#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <unordered_map>

namespace mio {

/**
 * @brief Brief description
 * @details
 * Longer description of the algorithm as well as example of use
 */
class RegressionFill : public ResamplingAlgorithms {
	public:
		RegressionFill(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs);

		void resample(const std::string& stationHash, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
        void resample(const std::string& stationHash, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                const std::vector<MeteoData>& vecM, MeteoData& md, const std::vector<METEO_SET>& additional_stations);
		std::string toString() const;

        double linear(double julian_date, const std::vector<double>& coefficients);

    private: 
        std::unordered_map<int, std::vector<double>> regression_coefficients;
};

} //end namespace mio

#endif
