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
#ifndef __UNVENTILATED_T_H__
#define __UNVENTILATED_T_H__

#include <meteoio/meteofilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcUnventilatedT
 * @ingroup processing
 * @brief Filters and corrects temperatures from unventilated sensor.
 * This implements the correction described in <i>"Air Temperature Measurement Errors in Naturally Ventilated Radiation Shields"</i>, Reina Nakamura, L. Mahrt, J. Atmos. Oceanic Technol., <b>22</b>, 2005, pp 1046–1058
 * with an albedo dependency as introduced in <i>"Albedo effect on radiative errors in air temperature measurements"</i>, H. Huwald, C. W. Higgins, M.-O. Boldi, E. Bou-Zeid, M. Lehning, and M. B. Parlange, Water Resour. Res., <b>45</b>, W08431, 2009.
 *
 * If the "soft" option is given, the albedo has a value different according to snow (or no snow) on the ground. The default albedo (that is also used if snow height is nodata) can be given as a second option.
 *
 * @note This filter can ONLY be applied to air temperatures. Moreover, since it <i>might</i> depend on the radiation shield, it is highly recommended to do some tests (ie. comparisons between ventillated and unventillated sensors) before using it on a new sensor type. Hopefully a new paper would come and clarify its usability per sensor types...
 *
 * @code
 * TA::filter2	= unventilated_T
 * TA::arg2	= soft 0.23
 * @endcode
 *
 * @author Mathias Bavay
 * @date   2012-05-03
 */

class ProcUnventilatedT : public ProcessingBlock {
	public:
		ProcUnventilatedT(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);

		double usr_albedo;
		static const double dflt_albedo, soil_albedo, snow_albedo;
		static const double snow_thresh, vw_thresh;
		bool is_soft;
		bool nakamura; //use Nakamura or Huwald model

};

} //end namespace

#endif
