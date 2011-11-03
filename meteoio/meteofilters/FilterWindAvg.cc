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
#include <cmath>
#include <meteoio/meteofilters/FilterWindAvg.h>

using namespace std;

namespace mio {

FilterWindAvg::FilterWindAvg(const std::vector<std::string>& vec_args) : WindowedFilter("WIND_AVG")
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void FilterWindAvg::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
                            std::vector<MeteoData>& ovec)
{
	if(index!=MeteoData::VW && index!=MeteoData::DW) {
		stringstream ss;
		ss << "Can not use WIND_AVG processing on " << MeteoData::getParameterName(index);
		throw InvalidArgumentException(ss.str(), AT);
	}
	ovec.clear();
	ovec.reserve(ivec.size());

	for (size_t ii=0; ii<ivec.size(); ii++){ //for every element in ivec, get a window
		ovec.push_back(ivec[ii]);
		double& value = ovec[ii](index);

		const vector<const MeteoData*>& vec_window = get_window(ii, ivec);

		if (is_soft){
			if (vec_window.size() > 0){
				value = calc_avg(index, vec_window);
			} else {
				value = IOUtils::nodata;
			}
		} else {
			if (vec_window.size() >= min_data_points){
				value = calc_avg(index, vec_window);
			} else {
				value = IOUtils::nodata;
			}
		}
	}
}

/**
 * @brief Actual algorithm to calculate the average value for all values in vec_window(index)
 * @param index The MeteoData parameter to be averaged (e.g. MeteoData::TA, etc)
 * @param vec_window A vector of pointers to MeteoData that shall be used for the averaging
 * @return A double either representing the average or IOUtils::nodata if averaging fails
 */
double FilterWindAvg::calc_avg(const unsigned int& index, const std::vector<const MeteoData*>& vec_window)
{
		const size_t vecSize = vec_window.size();
		double meanspeed     = IOUtils::nodata;
		double meandirection = IOUtils::nodata;

		if (vecSize == 0){
			return IOUtils::nodata;
		} else {
			//calculate ve and vn
			double ve=0.0, vn=0.0;
			for (size_t jj=0; jj<vecSize; jj++){
				ve += vec_window[jj]->operator()(MeteoData::VW) * sin(vec_window[jj]->operator()(MeteoData::DW) * M_PI / 180.); //turn into radians
				vn += vec_window[jj]->operator()(MeteoData::VW) * cos(vec_window[jj]->operator()(MeteoData::DW) * M_PI / 180.); //turn into radians
			}
			ve /= vecSize;
			vn /= vecSize;

			meanspeed = sqrt(ve*ve + vn*vn);
			meandirection = fmod( atan2(ve,vn) * 180. / M_PI + 360. , 360.); // turn into degrees [0;360)
		}

		if(index==MeteoData::VW)
			return meanspeed;
		else
			return meandirection;
}

void FilterWindAvg::parse_args(std::vector<std::string> vec_args)
{
	vector<double> filter_args;

	if (vec_args.size() > 2){
		is_soft = FilterBlock::is_soft(vec_args);
	}

	if (vec_args.size() > 2)
		centering = (WindowedFilter::Centering)WindowedFilter::get_centering(vec_args);

	FilterBlock::convert_args(2, 2, vec_args, filter_args);

	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT);
	}

	min_data_points = (unsigned int)floor(filter_args[0]);
	min_time_span = Duration(filter_args[1] / 86400.0, 0.);
}

}
