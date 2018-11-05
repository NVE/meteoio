/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PROCQUANTILEMAPPING_H
#define PROCQUANTILEMAPPING_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcQuantileMapping
 * @ingroup processing
 * @brief Quantile Mapping correction
 * @details The statistical distribution of the chosen parameter is computed and the multiplying factors provided as arguments
 * are used to correct each provided quantiles (see https://rcmes.jpl.nasa.gov/content/statistical-downscaling#download).
 * 
 * The quantiles must be provided as increasing decimal numbers between 0 (0%) and 1 (100%) and the full range must be covered.
 *
 * @code
 * TA::filter1    = QM
 * TA::arg1::corrections = correctionsTA.dat
 * @endcode
 *
 * Example of correction file (quantiles in the first column and correction factors as second column):
 * @code
 * 0 1.2
 * 0.25 1.2
 * 0.5 1.05
 * 0.8 0.95
 * 0.9 0.89
 * 1 0.89
 * @endcode
 */

class ProcQuantileMapping : public ProcessingBlock {
	public:
		ProcQuantileMapping(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& i_root_path);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	protected:
		void correctPeriod(const unsigned int& param, const size_t& idx_start, const size_t& idx_end, std::vector<MeteoData>& ovec) const;
		double getCorrection(const std::vector<double>& thresholds, const double& value) const;
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);

		std::vector<double> quantiles, corrections;
		std::string root_path;
		char period;
};

} //end namespace

#endif
