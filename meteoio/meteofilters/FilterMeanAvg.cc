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
#include <meteoio/meteofilters/FilterMeanAvg.h>
#include <cmath>

using namespace std;

namespace mio {

FilterMeanAvg::FilterMeanAvg(const std::vector<std::string>& vec_args) : WindowedFilter("MEAN_AVG")
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void FilterMeanAvg::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
					   std::vector<MeteoData>& ovec)
{
	ovec.clear();

	for (unsigned int ii=0; ii<ivec.size(); ii++){ //for every element in ivec, get a window
		ovec.push_back(ivec[ii]);
		double& value = ovec[ii].param(index);

		const vector<const MeteoData*>& vec_window = get_window(ii, ivec);

		//cout << "left : "<<elements_left << endl << "right: "<<elements_right << endl;
		//cout << "size : "<<vec_window.size() << endl;

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

		//cout << endl << "Final value: " << value << " (" << ovec[ii].date.toString(Date::ISO) << ")" << endl;
		//cout << "============================" << endl;
	}
}

/**
 * @brief Actual algorithm to calculate the average value for all values in vec_window.param(index)
 * @param index The MeteoData parameter to be averaged (e.g. MeteoData::TA, etc)
 * @param vec_window A vector of pointers to MeteoData that shall be used for the averaging
 * @return A double either representing the average or IOUtils::nodata if averaging fails 
 */
double FilterMeanAvg::calc_avg(const unsigned int& index, const std::vector<const MeteoData*>& vec_window)
{
	if (vec_window.size() == 0)
		return IOUtils::nodata;

	double sum = 0;
	unsigned int counter = 0;
	for (unsigned int ii=0; ii<vec_window.size(); ii++){
		const double& value = (*vec_window[ii]).param(index);
		if (value != IOUtils::nodata){
			sum += value;
			counter++;
			//cout << "Adding value: " << value << "(" << (*vec_window[ii]).date.toString(Date::ISO) << ")" << endl;
		}
	}

	if (counter == 0){
		return IOUtils::nodata;
	} else if (counter == 1){
		return sum;
	}

	return (sum / (double)counter);
}	

void FilterMeanAvg::parse_args(std::vector<std::string> vec_args)
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
	min_time_span = Date(filter_args[1] / 86400.0);
} 

}
