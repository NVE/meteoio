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
#include <meteoio/FilterMeanAvg.h>

using namespace std;

namespace mio {

FilterMeanAvg::FilterMeanAvg(const std::vector<std::string>& vec_args) : FilterBlock("MEAN_AVG")
{
	parse_args(vec_args);
}

void FilterMeanAvg::process(const unsigned int& index, const std::vector<MeteoData>& ivec,
					   std::vector<MeteoData>& ovec)
{
	ovec.clear();

	vector<const MeteoData*> vec_window;
	unsigned int elements_left=0, elements_right=0;
	unsigned int kk = 0;

	for (kk=0; kk<min_data_points; kk++){
		if (ivec.size() > kk) {
			elements_right++;
			vec_window.push_back(&ivec[kk]);
		}
	}
	
	if (elements_right > 0){
		elements_left = 1; //just as an initial setup

		if (centering == FilterBlock::left){
			vec_window.clear();
			vec_window.push_back(&ivec[0]);
			elements_right = 1;
		}
	}


	/*
	if (vec_window.size() > 0){
		while ((((*vec_window[vec_window.size()-1]).date - (*vec_window[0]).date) < min_time_span)
			  && ((ivec.size() > kk))){
			elements_right++;
			vec_window.push_back(&ivec[kk]);
		}
	}
	*/
	for (unsigned int ii=0; ii<ivec.size(); ii++){
		ovec.push_back(ivec[ii]);
		double& value = ovec[ii].param(index);

		//check whether a window is available, calculate point and output it
		if (is_soft){
			//wait a minute
		} else {
			if (centering == FilterBlock::right){
				if (elements_right >= min_data_points){
					value = calc_avg(index, vec_window);
					if (vec_window.size() > 0){
						vec_window.erase(vec_window.begin());
						elements_right--;
					}
					if (ivec.size() > (ii+min_data_points)) { //shift window one point to the right
						vec_window.push_back(&ivec[ii+min_data_points]);
						elements_right++; //elements_left will stay at a constant 1
					}
				} else {
					value = IOUtils::nodata;
				}
			} else if (centering == FilterBlock::left){
				if (elements_left >= min_data_points){
					value = calc_avg(index, vec_window);
					if (vec_window.size() > 0){
						vec_window.erase(vec_window.begin());
						elements_left--;
					}

					if (ivec.size() > (ii+1)) { //shift window one point to the right
						vec_window.push_back(&ivec[ii+1]);
						elements_left++; //elements_left will stay at a constant 1
					}
					
				} else {
					value = IOUtils::nodata;
					if (ivec.size() > (ii+1)) { //broaden window
						vec_window.push_back(&ivec[ii+1]);
						elements_left++; 
					}
				}
			} else if (centering == FilterBlock::center){
				throw IOException("Centered MeanAvgFilter currently not implemented!", AT);
			}
		}
	}
	
}

double FilterMeanAvg::calc_avg(const unsigned int& index, std::vector<const MeteoData*>& vec_window)
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

	if (vec_args.size() > 2){
		centering = (FilterBlock::WindowOrientation)FilterBlock::get_orientation(vec_args);
	}
	
	FilterBlock::convert_args(2, 2, vec_args, filter_args);

	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT); 
	}

	min_data_points = filter_args[0];
	min_time_span = Date(filter_args[1] / 86400.0);
} 

}
