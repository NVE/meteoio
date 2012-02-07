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
#include <meteoio/meteofilters/ProcUndercatch.h>
#include <cmath>

using namespace std;

namespace mio {

const double ProcUndercatch::Tsnow=-2., ProcUndercatch::Train=2.; //WMO values from Yan et al (2001)

ProcUndercatch::ProcUndercatch(const std::vector<std::string>& vec_args) : ProcessingBlock("UNDERCATCH") {
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcUndercatch::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	if(index!=MeteoData::HNW)
		throw InvalidArgumentException("Trying to use UNDERCATCH filter on " + MeteoData::getParameterName(index) + " but it can only be applied to precipitation!!" + getName(), AT);
	ovec.clear();
	ovec.reserve(ivec.size());

	for (unsigned int ii=0; ii<ivec.size(); ii++){
		ovec.push_back(ivec[ii]);

		double& tmp = ovec[ii](index);
		const double VW = ovec[ii](MeteoData::VW);
		double t = ovec[ii](MeteoData::TA);
		if(t!=IOUtils::nodata) t=K_TO_C(t); //t in celsius
		precip_type precip=mixed;
		if(t<=Tsnow) precip=snow;
		if(t>=Train) precip=rain;

		//HACK:: we don't use Tmax, Tmin, Tmean but only the current temperature instead
		if (tmp == IOUtils::nodata || tmp==0. || precip==rain) {
			continue; //preserve nodata values and no precip or purely liquid precip
		} else if(type==cst) {
			if(precip==snow) tmp *= factor_snow;
			if(precip==mixed) tmp *= factor_mixed;
		} else if(type==nipher) {
			if(VW==IOUtils::nodata) continue;
			double k;
			if(precip==snow) k=100.-0.44*VW*VW-1.98*VW;
			if(precip==mixed) {
				if(t==IOUtils::nodata) continue;
				k=97.29-3.18*VW+0.58*t-0.67*t; //Tmax, Tmin
			}
			tmp *= 100./k;
		} else if(type==tretyakov) {
			if(VW==IOUtils::nodata || t==IOUtils::nodata) continue;
			double k;
			if(precip==snow) k=103.11-8.67*VW+0.30*t; //Tmax
			if(precip==mixed) k=96.99-4.46*VW+0.88*t+0.22*t; //Tmax, Tmin
			tmp *= 100./k;
		} else if(type==us8sh) {
			if(VW==IOUtils::nodata) continue;
			double k;
			if(precip==snow) k=exp(4.61-0.04*pow(VW, 1.75));
			if(precip==mixed) k=101.04-5.62*VW;
			tmp *= 100./k;
		} else if(type==us8unsh) {
			if(VW==IOUtils::nodata) continue;
			double k;
			if(precip==snow) k=exp(4.61-0.16*pow(VW, 1.28));
			if(precip==mixed) k=100.77-8.34*VW;
			tmp *= 100./k;
		} else if(type==hellmann) {
			if(VW==IOUtils::nodata) continue;
			double k;
			if(precip==snow) k=100.+1.13*VW*VW-19.45*VW;
			if(precip==mixed) {
				if(t==IOUtils::nodata) continue;
				k=96.63+0.41*VW*VW-9.84*VW+5.95*t; //Tmean
			}
			tmp *= 100./k;
		}
	}
}

void ProcUndercatch::parse_args(std::vector<std::string> filter_args)
{
	if (filter_args.size() < 1)
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);

	for(size_t ii=0; ii<filter_args.size(); ii++) {
		IOUtils::toLower(filter_args[ii]);
	}

	if(filter_args[0]=="cst") {
		type=cst;
		if (filter_args.size() < 3)
			throw InvalidArgumentException("Wrong number of arguments for filter "+getName()+" with rain gauge type \"cst\"", AT);
		IOUtils::convertString(factor_snow, filter_args[1]);
		IOUtils::convertString(factor_mixed, filter_args[2]);
	} else if(filter_args[0]=="nipher") {
		type=nipher;
	} else if(filter_args[0]=="tretyakov") {
		type=tretyakov;
	} else if(filter_args[0]=="us8sh") {
		type=us8sh;
	} else if(filter_args[0]=="us8unsh") {
		type=us8unsh;
	} else if(filter_args[0]=="hellmann") {
		type=hellmann;
	} else {
		throw InvalidArgumentException("Rain gauge type \""+ filter_args[0] +"\" unknown for filter "+getName(), AT);
	}
}

} //end namespace
