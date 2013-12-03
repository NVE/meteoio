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
#include <meteoio/meteofilters/ProcUnshade.h>
#include <cmath>

using namespace std;

namespace mio {

ProcUnshade::ProcUnshade(const std::vector<std::string>& vec_args, const std::string& name) : WindowedFilter(name), max_gap(3.5)
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = Duration(.5, 0.);
	properties.time_after  = Duration(.5, 0.);
	properties.points_before = 1;
	properties.points_after = 1;
}

void ProcUnshade::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                           std::vector<MeteoData>& ovec)
{
	if(param!=MeteoData::ISWR && param!=MeteoData::RSWR)
		throw InvalidArgumentException("Trying to use "+getName()+" filter on " + MeteoData::getParameterName(param) + " but it can only be applied to iswr/rswr!!" + getName(), AT);
	
	ovec = ivec;
	vector<double> vecAlbedo, vecJulian;
	
	//build the vector of albedos
	for (size_t ii=0; ii<ovec.size(); ii++) {
		const double iswr = ovec[ii](MeteoData::ISWR);
		const double rswr = ovec[ii](MeteoData::RSWR);
		
		if(iswr<=0. && rswr<=0.) continue;
		
		const double albedo = (iswr!=IOUtils::nodata && rswr!=IOUtils::nodata)? rswr / (iswr+1e-6) : IOUtils::nodata;
		if(albedo>0. && albedo<1.)
			vecAlbedo.push_back( albedo );
		else
			vecAlbedo.push_back( IOUtils::nodata );
		
		vecJulian.push_back( ovec[ii].date.getJulian(true) );
	}

	//filter the albedo
	if( !MADFilter(vecAlbedo) ) return; //HACK: do it piecewise on a window!
	
	//reinterpolate the albedos
	if( !linInterpolate(vecJulian, vecAlbedo) ) return;
	
	//apply the albedo correction
	size_t alb_idx = 0;
	for (size_t ii=0; ii<ovec.size(); ii++) {
		double& iswr = ovec[ii](MeteoData::ISWR);
		double& rswr = ovec[ii](MeteoData::RSWR);
		
		if(iswr<=0. && rswr<=0.) continue; //skip nights
		
		const double albedo = vecAlbedo[alb_idx];
		
		if(param==MeteoData::ISWR) {
			if(rswr!=IOUtils::nodata && albedo!=IOUtils::nodata) {
				if(rswr>0. && albedo>0.)
					iswr = rswr / albedo;
				else 
					iswr = 0.;
			}
		} else { //on RSWR
			if(iswr!=IOUtils::nodata && albedo!=IOUtils::nodata) {
				if(iswr>0. && albedo>0.)
					rswr = iswr * albedo;
				else 
					rswr = 0.;
			}
		}
		
		alb_idx++;
	}
}

bool ProcUnshade::MADFilter(std::vector<double> &data)
{
	const double K = 1. / 0.6745;
	double mad     = IOUtils::nodata;
	double median  = IOUtils::nodata;
	
	try {
		median = Interpol1D::getMedian(data);
		mad    = Interpol1D::getMedianAverageDeviation(data);
	} catch(const exception&){
		return false;
	}

	if( median==IOUtils::nodata || mad==IOUtils::nodata ) return false;

	const double sigma = mad * K;
	const double lower_lim = median - 3.*sigma;
	const double upper_lim = median + 3.*sigma;
	
	for(size_t ii=0; ii<data.size(); ii++) {
		double& value = data[ii];
		if( value!=IOUtils::nodata && ((value>upper_lim) || (value<lower_lim)) ) {
			value = IOUtils::nodata;
		}
	}
	
	return true;
}

void ProcUnshade::interpolFill(const size_t& start_idx, const size_t& end_idx, const std::vector<double>& julian, std::vector<double> &data)
{
	const double x1 = julian[start_idx];
	const double x2 = julian[end_idx];
	const double y1 = data[start_idx];
	const double y2 = data[end_idx];
	
	if (x1 == x2)
		throw IOException("Attempted division by zero", AT);

	//Solving y = ax + b
	const double a = (y2 - y1) / (x2 - x1);
	const double b = y2 - a*x2;
	
	for(size_t ii = start_idx+1; ii<end_idx; ii++)
		data[ii] = a * julian[ii] + b;
}

bool ProcUnshade::linInterpolate(const std::vector<double>& julian, std::vector<double> &data)
{
	//get starting point
	size_t start_idx = IOUtils::npos;
	for(size_t ii=0; ii<data.size(); ii++) {
		if(data[ii]!=IOUtils::nodata) {
			start_idx = ii;
			break;
		}
	}
	
	if(start_idx==IOUtils::npos) return false;
	
	for(size_t ii=start_idx+1; ii<data.size(); ii++) {
		if(data[ii]!=IOUtils::nodata) {
			if( ii!=start_idx+1 && (julian[ii]-julian[start_idx])<=max_gap )
				interpolFill(start_idx, ii, julian, data);
			start_idx = ii;
		}
	}
	
	return true;
}

void ProcUnshade::parse_args(std::vector<std::string> vec_args)
{
	if (vec_args.size() > 1){
		throw InvalidArgumentException("Invalid arguments for filter " + getName(), AT);
	}
	
	vector<double> filter_args;
	
	if(vec_args.size()==1) {
		IOUtils::convertString( max_gap, vec_args[0]);
		max_gap /= 86400.; //convert back to days
	}
}

}
