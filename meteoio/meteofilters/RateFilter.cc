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
#include <meteoio/meteofilters/RateFilter.h>
#include <cmath>

using namespace std;

namespace mio {

RateFilter::RateFilter(const std::vector<std::string>& vec_args) : FilterBlock("RATE")
{
	parse_args(vec_args);
	properties.for_second_pass = true; //for the rest: default values
}

void RateFilter::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
                           std::vector<MeteoData>& ovec)
{
	ovec.clear();
	size_t last_good = IOUtils::npos;

	//Find first point that is not IOUtils::nodata
	for (size_t ii=0; ii<ivec.size(); ii++){
		ovec.push_back(ivec[ii]);
		if (ovec[ii].param(index) != IOUtils::nodata){
			last_good = ii;
			break;
		}
	}

	if (last_good == IOUtils::npos) //can not find a good point to start
		return;

	for (size_t ii=(last_good+1); ii<ivec.size(); ii++) {
		ovec.push_back(ivec[ii]);

		double& curr_value       = ovec[ii].param(index);
		const double& prev_value = ovec[last_good].param(index);
		const double curr_time   = ovec[ii].date.getJulianDate();
		const double prev_time   = ovec[last_good].date.getJulianDate();

		if (curr_value == IOUtils::nodata)
			continue;

		const double local_rate = (curr_value-prev_value) / ((curr_time-prev_time)*24.*3600.); //per seconds

		if( local_rate>max_rate_of_change || local_rate<min_rate_of_change ) {
			curr_value = IOUtils::nodata;
		} else {
			last_good = ii;
		}
	}
}

void RateFilter::parse_args(std::vector<std::string> vec_args) {
	vector<double> filter_args;
	FilterBlock::convert_args(1, 2, vec_args, filter_args);

	const size_t nb_args = filter_args.size();
	if (nb_args == 2) {
		min_rate_of_change = filter_args[0];
		max_rate_of_change = filter_args[1];
	} else if(nb_args == 1) {
		min_rate_of_change = -filter_args[0];
		max_rate_of_change = filter_args[0];
	} else
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName() + " - Please provide 1 or 2 arguments!", AT);

}

}
