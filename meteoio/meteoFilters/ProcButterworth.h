/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PROCBUTTERWORTH_H
#define PROCBUTTERWORTH_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcButterworth
 * @ingroup processing
 * @brief Simple 2 poles Butterworth low pass filter.
 * The cutoff <b>period</b> (defined as the frequency at a -3dB gain) is given in seconds as argument. The phase is removed by running
 * the filter twice, first backward and then forward (this also squares the amplitude response, ie doubles the order of the filter, see
 * http://www.dspguide.com/ch19/4.htm or http://unicorn.us.com/trading/allpolefilters.html)
 *
 * The original paper is <i>On the Theory of Filters Amplifiers</i>, S. Butterworth, Experimental wireless & the wireless engineer,
 * <b>7</b>, pp 536-541, 1930.
 * @code
 * VW::filter1	= BUTTERWORTH
 * VW::arg1	= 10800 ;3 hours
 * @endcode
 *
 */

class ProcButterworth : public ProcessingBlock {
	public:
		ProcButterworth(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void computeCoefficients(const double& samplerate, const double& f_cutoff, double A[3], double B[3]) const;
		void parse_args(std::vector<std::string> vec_args);

		double cutoff;
		bool bidirectional;

		static const double n, p, g, c; ///< filter definition: number of passes, polynomial coefficients, 3dB cutoff correction
};

} //end namespace

#endif
