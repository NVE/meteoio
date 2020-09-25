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
#include <meteoio/meteoFilters/FilterMinMax.h>

using namespace std;

namespace mio {

FilterMinMax::FilterMinMax(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
             : ProcessingBlock(vecArgs, name, cfg), min_val(0.), max_val(0.), min_soft(0.), max_soft(0.), is_soft(false)

{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::both; //for the rest: default values
}

void FilterMinMax::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                           std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	for (size_t ii=0; ii<ovec.size(); ii++){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		if (tmp < min_val){
			if (is_soft){
				tmp = min_soft;
			} else {
				tmp = IOUtils::nodata;
			}
		} else if (tmp > max_val){
			if (is_soft){
				tmp = max_soft;
			} else {
				tmp = IOUtils::nodata;
			}
		}
	}
}

void FilterMinMax::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	bool has_min=false, has_min_reset=false;
	bool has_max=false, has_max_reset=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="SOFT") {
			IOUtils::parseArg(vecArgs[ii], where, is_soft);
		} else if (vecArgs[ii].first=="MAX") {
			IOUtils::parseArg(vecArgs[ii], where, max_val);
			has_max = true;
		} else if (vecArgs[ii].first=="MAX_RESET") {
			IOUtils::parseArg(vecArgs[ii], where, max_soft);
			has_max_reset = true;
		} else if (vecArgs[ii].first=="MIN") {
			IOUtils::parseArg(vecArgs[ii], where, min_val);
			has_min = true;
		} else if (vecArgs[ii].first=="MIN_RESET") {
			IOUtils::parseArg(vecArgs[ii], where, min_soft);
			has_min_reset = true;
		}
	}

	if (!has_max) throw InvalidArgumentException("Please provide a MAX value for "+where, AT);
	if (!has_min) throw InvalidArgumentException("Please provide a MIN value for "+where, AT);
	if (is_soft && !has_max_reset) max_soft = max_val;
	if (is_soft && !has_min_reset) min_soft = min_val;
}

} //end namespace
