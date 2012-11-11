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
#include <meteoio/meteofilters/ProcMult.h>

using namespace std;

namespace mio {

ProcMult::ProcMult(const std::vector<std::string>& vec_args) : ProcessingBlock("MULT"), factor(1.)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcMult::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec.clear();
	ovec.reserve(ivec.size());

	for (size_t ii=0; ii<ivec.size(); ii++){
		ovec.push_back(ivec[ii]);

		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		tmp *= factor;
	}
}


void ProcMult::parse_args(const std::vector<std::string>& vec_args) {
	vector<double> filter_args;

	convert_args(1, 1, vec_args, filter_args);

	if (filter_args.size() > 1)
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);

	factor = filter_args[0];
}

} //end namespace
