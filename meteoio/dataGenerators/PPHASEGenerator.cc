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

#include <meteoio/dataGenerators/PPHASEGenerator.h>

namespace mio {

void PPhaseGenerator::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	bool has_type=false, has_thresh=false, has_thresh1=false, has_thresh2=false;
	double range_thresh1, range_thresh2;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string user_algo( IOUtils::strToUpper(vecArgs[ii].second) );

			if (user_algo=="THRESH") model = THRESH;
			else if (user_algo=="RANGE") model = RANGE;
			else
				throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for the "+algo+" generator", AT);

			has_type = true;
		} else if(vecArgs[ii].first=="THRESH") {
			parseArg(vecArgs[ii], fixed_thresh);
			has_thresh = true;
		} else if(vecArgs[ii].first=="SNOW") {
			parseArg(vecArgs[ii], range_thresh1);
			has_thresh1 = true;
		} else if(vecArgs[ii].first=="RAIN") {
			parseArg(vecArgs[ii], range_thresh2);
			has_thresh2 = true;
		}
	}

	if (!has_type) throw InvalidArgumentException("Please provide a TYPE for algorithm "+algo, AT);
	if (model == THRESH && !has_thresh) throw InvalidArgumentException("Please provide a threshold for algorithm "+algo, AT);
	if (model == RANGE) {
		if (!has_thresh1 || !has_thresh2) throw InvalidArgumentException("Please provide a a snow and a rain threshold for algorithm "+algo, AT);
		range_start = range_thresh1;
		range_norm = 1. / (range_thresh2-range_thresh1);
	}
}

bool PPhaseGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA);
		if (TA==IOUtils::nodata) return false;
		
		if (model==THRESH) {
			value = (TA>=fixed_thresh)? 1. : 0.;
		} else if (model==RANGE) {
			const double tmp_rainfraction = range_norm * (TA - range_start);
			value = (tmp_rainfraction>1)? 1. : (tmp_rainfraction<0.)? 0. : tmp_rainfraction;
		}
	}

	return true; //all missing values could be filled
}

bool PPhaseGenerator::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;
	
	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii]))
			all_filled = false;
	}

	return all_filled;
}

} //namespace
