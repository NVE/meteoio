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
#ifndef __PROCPASSIVET_H__
#define __PROCPASSIVET_H__

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcPassiveT
 * @ingroup processing
 * @brief Filters and corrects temperatures from unventilated sensor.
 * This implements the correction described in <i>"Air Temperature Measurement Errors in Naturally Ventilated Radiation Shields"</i>, Reina Nakamura, L. Mahrt, J. Atmos. Oceanic Technol., <b>22</b>, 2005, pp 1046–1058
 * with an albedo dependency as introduced in <i>"Albedo effect on radiative errors in air temperature measurements"</i>, H. Huwald, C. W. Higgins, M.-O. Boldi, E. Bou-Zeid, M. Lehning, and M. B. Parlange, Water Resour. Res., <b>45</b>, W08431, 2009.
 *
 * If the "soft" option is given, the albedo has a value different according to snow (or no snow) on the ground. The default albedo (that is also used if snow height is nodata) can be given as a second option.
 *
 * @note This filter can ONLY be applied to air temperatures.
 * @code
 * TA::filter2	= passive_T
 * TA::arg2	= soft 0.23
 * @endcode
 *
 * @author Mathias Bavay
 * @date   2012-05-03
 */

class ProcPassiveT : public ProcessingBlock {
	public:
		ProcPassiveT(const std::vector<std::string>& vec_args);

		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);

		bool is_soft;
		double albedo;
		static const double dflt_albedo, soil_albedo, snow_albedo;
		static const double snow_thresh;
};

} //end namespace

#endif
