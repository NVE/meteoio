/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef NOISE_H
#define NOISE_H

#include <meteoio/meteoFilters/FilterBlock.h>

#include <vector>
#include <string>

namespace mio {

/**
 * @class ProcNoise
 * @ingroup processing
 * @brief Generate a noise signal to modify the input. 
 * @details
 * The noise signal is either added ("add")  to the input or used as a fraction and multiplied by the input signal ("mult").
 * This filter always takes three arguments:
 *  - RANGE: the scaling factor to apply;
 *  - TYPE: to specify if the noise is added (ADD) or multiplied (MULT);
 *  - DISTRIBUTION: to specify the random numbers distribution as either
 *      + uniform: the range represents the maximum amplitude of the noise;
 *      + normal: the range represents the standard deviation of the noise.
 *
 * @note When multiplying the input signal by the noise, the range is interpreted as a fraction. For example, to modify the input by
 * +/-10%, select "mult" with "uniform" and "0.1" as a range.
 * @code
 *  #add a normally distributed noise (mean=0) of standard deviation 5 to TA
 * TA::filter1            = Noise
 * TA::arg1::type         = add
 * TA::arg1::distribution = normal
 * TA::arg1::range        = 5
 * 
 * #add +/-10% of variation, uniformly distributed to HS
 * HS::filter1            = Noise
 * HS::arg1::type         = mult
 * HS::arg1::distribution = uniform
 * HS::arg1::range        = 0.1
 * @endcode
 */

class ProcNoise : public ProcessingBlock {
	public:
		ProcNoise(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void uniform_add(const unsigned int& param, std::vector<MeteoData>& ovec) const;
		void uniform_mult(const unsigned int& param, std::vector<MeteoData>& ovec) const;
		void normal_add(const unsigned int& param, std::vector<MeteoData>& ovec) const;
		void normal_mult(const unsigned int& param, std::vector<MeteoData>& ovec) const;
		double getBoxMuller() const;
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		
		double range;
		char distribution, type;
};

} //end namespace

#endif
