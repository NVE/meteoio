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
#include <meteoio/meteofilters/FilterMAD.h>
#include <meteoio/meteostats/libinterpol1D.h>
#include <cmath>

using namespace std;

namespace mio {

FilterMAD::FilterMAD(const std::vector<std::string>& vec_args) : WindowedFilter("MAD")
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise: //HACK: does not account for centering!
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void FilterMAD::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec.clear();
	ovec.reserve(ivec.size());

	for (unsigned int ii=0; ii<ivec.size(); ii++){ //for every element in ivec, get a window
		ovec.push_back(ivec[ii]);
		double& value = ovec[ii](index);

		const vector<const MeteoData*>& vec_window = get_window(ii, ivec);
		if(value==IOUtils::nodata) continue; //because get_window needs to move the index 1 by 1

		if (is_soft){
			if (vec_window.size() > 0){
				MAD_filter_point(vec_window, index, value);
			}
		} else {
			if (vec_window.size() >= min_data_points){
				MAD_filter_point(vec_window, index, value);
			} else {
				value = IOUtils::nodata;
			}
		}

	}
}

void FilterMAD::MAD_filter_point(const std::vector<const MeteoData*>& vec_window, const unsigned int& index, double& value)
{
	const double K = 1. / 0.6745;
	double mad     = IOUtils::nodata;
	double median  = IOUtils::nodata;

	std::vector<double> data;
	for(unsigned int ii=0; ii<vec_window.size(); ii++) data.push_back( (*vec_window[ii])(index) );

	//Calculate MAD
	try {
		median = Interpol1D::getMedian(data);
		mad    = Interpol1D::getMedianAverageDeviation(data);
	} catch(const exception&){
		return;
	}

	if( median==IOUtils::nodata || mad==IOUtils::nodata ) return;

	const double sigma = mad * K;
	const double upper_lim = median + 3.*sigma;
	const double lower_lim = median - 3.*sigma;

	//cout << lower_lim << " " << value << " " << upper_lim << endl;
	if( (value>upper_lim) || (value<lower_lim) ) {
		value = IOUtils::nodata;
	}
}

void FilterMAD::parse_args(std::vector<std::string> vec_args)
{
	vector<double> filter_args;

	if (vec_args.size() > 2){
		is_soft = FilterBlock::is_soft(vec_args);
	}

	if (vec_args.size() > 2)
		centering = (WindowedFilter::Centering)WindowedFilter::get_centering(vec_args);

	convert_args(2, 2, vec_args, filter_args);

	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT);
	}

	min_data_points = (unsigned int)floor(filter_args[0]);
	min_time_span = Duration(filter_args[1] / 86400.0, 0.);
}

}
