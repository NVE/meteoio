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
#include <meteoio/meteofilters/ProcPassiveT.h>
#include <cmath>

using namespace std;

namespace mio {

const double ProcPassiveT::dflt_albedo = .23;
const double ProcPassiveT::soil_albedo = .23; //grass
const double ProcPassiveT::snow_albedo = .56; //white surface
const double ProcPassiveT::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo
const double ProcPassiveT::vw_thresh = 0.1; //wind speed threshold
const bool ProcPassiveT::nakamura = false; //use Nakamura or Huwald model

ProcPassiveT::ProcPassiveT(const std::vector<std::string>& vec_args) : ProcessingBlock("PASSIVE_T") {
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcPassiveT::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
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

		double albedo = usr_albedo;
		double iswr = ivec[ii](MeteoData::ISWR);
		const double rswr = ivec[ii](MeteoData::RSWR);
		const double ta = ivec[ii](MeteoData::TA);
		double vw = ivec[ii](MeteoData::VW);
		double hs = ivec[ii](MeteoData::HS);

		if(is_soft && iswr!=IOUtils::nodata && rswr!=IOUtils::nodata && rswr>5. && iswr>5.) {
			albedo = iswr / rswr;
			hs = IOUtils::nodata; //to make sure we would not try to recompute a pseudo albedo later
		}

		if(is_soft && hs!=IOUtils::nodata && iswr==IOUtils::nodata && rswr!=IOUtils::nodata) {
			if(hs>snow_thresh) iswr = snow_albedo*rswr;
			else iswr = soil_albedo*rswr;
		}

		if(iswr==IOUtils::nodata || ta==IOUtils::nodata || vw==IOUtils::nodata)
			continue;

		if(is_soft && hs!=IOUtils::nodata) { //try to get snow height in order to adjust the albedo
			if(hs>snow_thresh) albedo = snow_albedo;
			else albedo = soil_albedo;
		}

		if(vw<vw_thresh) vw=vw_thresh; //this should be around the minimum measurable wind speed on regular instruments
		const double rho = 1.2; // in kg/m3
		const double Cp = 1004.;
		const double X = iswr / (rho*Cp*ta*vw);
		if(X<1e-4) continue; //the correction does not work well for small X values
		if(nakamura) {
			const double C0 = 0.13;
			const double C1 = 373.40 * albedo / dflt_albedo; //in order to introduce the albedo as a scaling factor
			const double RE = C0 + C1*X;
			tmp -= RE; //substracting the radiative error
		} else {
			const double RE = 3.1 * sqrt(X);
			tmp -= RE; //substracting the radiative error
		}
	}
}

void ProcPassiveT::parse_args(std::vector<std::string> vec_args) {
	vector<double> filter_args;

	is_soft = false;
	if (vec_args.size() >= 1){
		is_soft = FilterBlock::is_soft(vec_args);
	}

	convert_args(0, 1, vec_args, filter_args);

	if (filter_args.size() > 1)
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);

	if (filter_args.size() == 0){
		usr_albedo = dflt_albedo;
	} else {
		usr_albedo = filter_args[0];
	}
}

} //end namespace
