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
#ifndef PROCSHADE_H
#define PROCSHADE_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <meteoio/meteoLaws/Sun.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcShade
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2016-01-22
 * @brief Apply a shading mask to the Incoming or Reflected Short Wave Radiation
 * A shading mask (that will be linearly interpolated between the provided points) must be provided in a separate file,
 * simply containing the horizon elevation (in deg.) as a function of azimuth (in deg.):
 * @code
 * 0	5
 * 15	25
 * 45	12
 * 180	30
 * 270	20
 * @endcode
 * 
 * Then the filter is declared with the file name containing the horizon mask as argument:
 * @code
 * ISWR::filter1 = SHADE
 * ISWR::arg1    = ../input/iswr_mask.dat
 * @endcode
 *
 */

class ProcShade : public ProcessingBlock {
	public:
		ProcShade(const std::vector<std::string>& vec_args, const std::string& name, const std::string& i_root_path);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		static void readMask(const std::string& filter, const std::string& filename, std::vector< std::pair<double,double> > &o_mask);
		
		void parse_args(const std::vector<std::string>& vec_args);
		double getMaskElevation(const double& azimuth) const;

		std::map<std::string, SunObject> Suns;
		std::vector< std::pair<double,double> > mask;
		std::string root_path;
		
		static const double soil_albedo, snow_albedo, snow_thresh; ///< parametrize the albedo from HS
};

} //end namespace

#endif