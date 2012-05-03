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
#include <meteoio/meteofilters/FilterPassiveT.h>

using namespace std;

namespace mio {

FilterPassiveT::FilterPassiveT(const std::vector<std::string>& vec_args) : FilterBlock("PASSIVE_T") {
	parse_args(vec_args);
	properties.stage = ProcessingProperties::second; //for the rest: default values
}

void FilterPassiveT::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	if(index!=MeteoData::TA) {
		stringstream ss;
		ss << "Can not use " << getName() << " processing on " << MeteoData::getParameterName(index);
		throw InvalidArgumentException(ss.str(), AT);
	}
	ovec.clear();
	ovec.reserve(ivec.size());

	for (unsigned int ii=0; ii<ivec.size(); ii++){
		ovec.push_back(ivec[ii]);

		double& tmp = ovec[ii](index);
		if(tmp == IOUtils::nodata) continue; //preserve nodata values

		if(type==nakamura) {
			const double iswr = ivec[ii](MeteoData::ISWR);
			const double ta = ivec[ii](MeteoData::TA);
			const double vw = ivec[ii](MeteoData::VW);

			if(iswr==IOUtils::nodata || ta==IOUtils::nodata || vw==IOUtils::nodata)
				continue;

			const double rho = 1200.;
			const double Cp = 1004.;
			const double C0 = 0.13;
			const double C1 = 373.40;
			const double X = iswr / (rho*Cp*ta*vw);

			if(X<1e-4) continue; //the correction does not work well for small X values
			const double RE = C0 + C1*X;
			tmp += RE;
		}
	}
}


void FilterPassiveT::parse_args(std::vector<std::string> filter_args) {

	for(size_t ii=0; ii<filter_args.size(); ii++) {
		IOUtils::toLower(filter_args[ii]);
	}

	if(filter_args[0]=="nakamura") {
		type=nakamura;
	} else {
		throw InvalidArgumentException("Air temperature correction \""+ filter_args[0] +"\" unknown for filter "+getName(), AT);
	}

}

} //end namespace
