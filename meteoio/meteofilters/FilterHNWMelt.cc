/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteofilters/FilterHNWMelt.h>

using namespace std;

namespace mio {

FilterHNWMelt::FilterHNWMelt(const std::vector<std::string>& vec_args) : FilterBlock("HNW_MELT") {
	parse_args(vec_args);
	properties.stage = ProcessingProperties::both; //for the rest: default values
}

void FilterHNWMelt::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	if(index!=MeteoData::HNW) {
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

		if(tmp>0.) {
			const double rh = ivec[ii](MeteoData::RH);
			const double ta = ivec[ii](MeteoData::TA);
			const double tss = ivec[ii](MeteoData::TSS);

			if (rh!=IOUtils::nodata &&  rh<thresh_rh) //not enough humidity for precipitation
				tmp = 0.;
			if (ta!=IOUtils::nodata && tss!=IOUtils::nodata && (ta-tss)>thresh_Dt ) //clear sky condition
				tmp = 0.;
                        if ( rh==IOUtils::nodata && ((ta==IOUtils::nodata) || (tss==IOUtils::nodata)) )
                                tmp = IOUtils::nodata; //we could not even try to validate the data point -> we delete it for safety
		}
	}
}


void FilterHNWMelt::parse_args(const std::vector<std::string>& vec_args) {
	vector<double> filter_args;

	convert_args(0, 2, vec_args, filter_args);

	const size_t nb_args = filter_args.size();
	if (nb_args == 0) {
		thresh_rh = 0.5;
		thresh_Dt = 3.0;
	} else if(nb_args == 2) {
		thresh_rh = filter_args[0];
		thresh_Dt = filter_args[1];
	} else
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName() + " - Please provide 0 or 2 arguments!", AT);


}

} //end namespace
