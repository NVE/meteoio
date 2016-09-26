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
#ifndef PROCIIRLOWPASS_H
#define PROCIIRLOWPASS_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcIIRLowPass
 * @ingroup processing
 * @brief Two poles Critically Damped low pass filter.
 * The cutoff <b>period</b> (defined as the frequency at a -3dB gain) is given in seconds as argument. The phase is removed by
 * bidirectional filtering, ie. running the filter twice, first backward and then forward (this also squares the amplitude response, see
 * http://www.dspguide.com/ch19/4.htm or http://unicorn.us.com/trading/allpolefilters.html) with a mechanism to disable filtering for
 * the points that generate overshooting (basically, it runs forward/backward as well as backward/forward and removes that points that
 * are too different between the two since these are indicators of overshooting).
 *
 * @code
 * VW::filter1	= Low_Pass
 * VW::arg1	= 10800 ;3 hours
 * @endcode
 *
 * It is possible to disable the bidirectional filtering by adding the "single_pass" argument.
 */

class ProcIIRLowPass : public ProcessingBlock {
	public:
		ProcIIRLowPass(const std::vector<std::string>& vec_args, const std::string& name);

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
