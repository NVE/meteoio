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
#include <meteoio/meteofilters/WindowedFilter.h>

using namespace std;

namespace mio {

WindowedFilter::WindowedFilter(const std::string& name)
	: FilterBlock(name), is_soft(false), min_data_points(1), min_time_span(0.0, 0.),
	  centering(WindowedFilter::center), elements_left(0), elements_right(0), last_index(IOUtils::npos)
{}

unsigned int WindowedFilter::get_centering(std::vector<std::string>& vec_args)
{
	if (vec_args.size() > 0){
		if (vec_args[0] == "left"){
			vec_args.erase(vec_args.begin());
			return WindowedFilter::left;
		} else if (vec_args[0] == "right"){
			vec_args.erase(vec_args.begin());
			return WindowedFilter::right;
		} else if (vec_args[0] == "center"){
			vec_args.erase(vec_args.begin());
			return WindowedFilter::center;
		}
	}

	return WindowedFilter::center; //the default
}

/**
 * @brief A function that cuts out the desired window for the 'index' element within
 *        ivec, the window elements are stored into vec_window
 *        Calls to this function have to start with index 0, then 1, 2, 3, ...
 *        vec_window is not allowed to be changed between two calls
 * @param index The index of the element in ivec that requires a window
 * @param ivec The original sequence of data points
 */
const std::vector<const MeteoData*>& WindowedFilter::get_window(const size_t& index,
                                                                const std::vector<MeteoData>& ivec)
{
	if ((index == 0) || (last_index > index)){ //reset global variables
		vec_window.clear();
		elements_left = elements_right = 0;

		if ((centering == WindowedFilter::right) || (is_soft)){
			const size_t end_kk = MIN(min_data_points, ivec.size());
			//vec_window.reserve(end_kk+1);
			for (size_t kk=0; kk<end_kk; kk++){
				elements_right++;
				vec_window.push_back(&ivec[kk]);
			}

			if (elements_right > 0) elements_left = 1;
		} else if ((centering == WindowedFilter::left) || (centering == WindowedFilter::center)){
			if (ivec.size() > 0) {
				elements_left = elements_right = 1;
				vec_window.push_back(&ivec[0]); //HACK ?? this leads to 1 element in window!
			}
		}

		last_index = 0;
		return vec_window;
	}

	if (index != (last_index+1))
		throw IOException("get_window function only to be used with increments of 1 for the index", AT);

	//check whether a window is available
	if (is_soft){
		if (vec_window.size() > 0){
			if (centering == WindowedFilter::right){
				//Try to move right, if it doesn't work, don't change anything
				if (ivec.size() > (index+vec_window.size()-1)) { //shift window one point to the right
					vec_window.push_back(&ivec[index+vec_window.size()-1]);
					vec_window.erase(vec_window.begin());
				}
			} else if (centering == WindowedFilter::left){
				if (index >= (vec_window.size())){
					if (ivec.size() > index){ //otherwise don't touch the whole thing
						vec_window.erase(vec_window.begin());
						vec_window.push_back(&ivec[index]);
					}
				} else {
					elements_left++;
					elements_right--;
				}
			} else if (centering == WindowedFilter::center){
				if (elements_right <= elements_left){
					if (ivec.size() > (index+elements_right-1)){ //otherwise don't touch the whole thing
						vec_window.push_back(&ivec[index+elements_right-1]);
						vec_window.erase(vec_window.begin());
					} else {
						elements_right--;
						elements_left++;
					}
				} else {
					elements_right--;
					elements_left++;
				}
			}
		} //HACK what happens if vec_window is empty?
	} else { //!is_soft
		if (centering == WindowedFilter::right){
			if (elements_right >= min_data_points){
				vec_window.erase(vec_window.begin());
				elements_right--;

				if (ivec.size() > (index+min_data_points-1)) { //shift window one point to the right
					vec_window.push_back(&ivec[index+min_data_points-1]);
					elements_right++; //elements_left will stay at a constant 1
				}
			}
		} else if (centering == WindowedFilter::left){
			if (elements_left >= min_data_points){
				vec_window.erase(vec_window.begin());
				elements_left--;

				if (ivec.size() > index) { //shift window one point to the right
					vec_window.push_back(&ivec[index]);
					elements_left++; //elements_left will stay at a constant 1
				}
			} else {
				if (ivec.size() > index) { //broaden window
					vec_window.push_back(&ivec[index]);
					elements_left++;
				}
			}
		} else if (centering == WindowedFilter::center){
			if ((elements_left + elements_right - 1) >= min_data_points){
				vec_window.erase(vec_window.begin());
				if (elements_right > 0) elements_right--;

				if (ivec.size() > (index+elements_left-1)) { //shift window one point to the right
					vec_window.push_back(&ivec[index+elements_left-1]);
					elements_right++; //elements_left will stay at a constant
				}
			} else {
				if (ivec.size() > (index+elements_left-1)) { //shift window one point to the right
					vec_window.push_back(&ivec[index+elements_left-1]);
					elements_left++;
				}

				if ((elements_left + elements_right - 1) < min_data_points){
					if (ivec.size() > (index+elements_left-1)){ //shift window one point to the right
						vec_window.push_back(&ivec[index+elements_left-1]);
						elements_right++;
					}
				}
			}
		}
	}

	last_index = index;

	return vec_window;
}

/**
 * @brief A function that computes the start and end for a window for the 'index' element from ivec
 * @param index The index of the element in ivec that requires a window
 * @param ivec The original sequence of data points
 * @param start the start index of the window
 * @param end the end index of the window
 * @return true if success, false if a window could not be computed
 */
bool WindowedFilter::get_window_specs(const size_t& index, const std::vector<MeteoData>& ivec, size_t& start, size_t& end)
{	/*
	The principle is too compute the first index that matches the minimum number of points criteria,
	the the one that matches the minimum time window,
	then combine them (with the equivalent of OR: we take the MIN index).
	Afterward, we compute the last index [...] for number of points
	and the last index [...] for the time window
	and combine them (with the equivalent of OR: we take the MIN index).
	(or vice versa for right centering)
	*/
	const Date date = ivec[index].date;
	start = end = index; //for proper initial value, so we can bail out without worries

	if(centering == WindowedFilter::left) {
		//get start of window
		size_t start_elements = min_data_points - 1; //start elements criteria
		if(start_elements>index) {
			if(!is_soft) return false;
			start_elements = index; //as many as possible
		}
		const Date start_date = date - min_time_span;
		size_t start_time_idx = IOUtils::seek(start_date, ivec, false); //start time criteria
		if(start_time_idx==IOUtils::npos) {
			if(!is_soft) return false;
			start_time_idx=0; //first possible element
		}
		const size_t elements_left = MAX(index - start_time_idx, start_elements);
		start = index - elements_left;

		//get end of window
		if(!is_soft) return true; //with end=index
		size_t end_elements = (min_data_points>(elements_left+1))?min_data_points - (elements_left + 1):0;
		const Date end_date = ivec[start].date+min_time_span;
		size_t end_time_idx = (end_date>date)?IOUtils::seek(end_date, ivec, false):index; //end time criteria
		if(end_time_idx==IOUtils::npos) {
			if(!is_soft) return false;
			end_time_idx=ivec.size()-1; //last possible element
		}
		const size_t elements_right = MAX(end_time_idx - index, end_elements);
		end = index + elements_right;
	}

	if(centering == WindowedFilter::right) {
		//get end of window
		size_t end_elements = min_data_points - 1; //end elements criteria
		if(end_elements>(ivec.size()-1-index)) {
			if(!is_soft) return false;
			end_elements = (ivec.size()-1-index); //as many as possible
		}
		const Date end_date = date+min_time_span;
		size_t end_time_idx = IOUtils::seek(end_date, ivec, false); //end time criteria
		if(end_time_idx==IOUtils::npos) {
			if(!is_soft) return false;
			end_time_idx=ivec.size()-1; //last possible element
		}
		const size_t elements_right = MAX(end_time_idx - index, end_elements);
		end = index + elements_right;

		//get start of window
		if(!is_soft) return true; //with start=index
		size_t start_elements = (min_data_points>(elements_right+1))?min_data_points - (elements_right + 1):0;
		const Date start_date = ivec[end].date-min_time_span;
		size_t start_time_idx = (start_date<date)?IOUtils::seek(start_date, ivec, false):index; //start time criteria
		if(start_time_idx==IOUtils::npos) {
			if(!is_soft) return false;
			end_time_idx=0; //first possible element
		}
		const size_t elements_left = MAX(index - start_time_idx, start_elements);
		start = index - elements_left;
	}

	if(centering == WindowedFilter::center) {
		//get start of ideal window
		size_t start_elements = min_data_points/2; //start elements criteria
		if(start_elements>index) {
			if(!is_soft) return false;
			start_elements = index; //as many as possible
		}
		const Date start_date = date - min_time_span/2;
		size_t start_time_idx = IOUtils::seek(start_date, ivec, false); //start time criteria
		if(start_time_idx==IOUtils::npos) {
			if(!is_soft) return false;
			start_time_idx=0; //first possible element
		}
		const size_t elements_left = MAX(index - start_time_idx, start_elements);
		start = index - elements_left;

		//get end of ideal window
		size_t end_elements = min_data_points/2; //end elements criteria
		if(end_elements>(ivec.size()-1-index)) {
			if(!is_soft) return false;
			end_elements = (ivec.size()-1-index); //as many as possible
		}
		const Date end_date = date+min_time_span/2;
		size_t end_time_idx = IOUtils::seek(end_date, ivec, false); //end time criteria
		if(end_time_idx==IOUtils::npos) {
			if(!is_soft) return false;
			end_time_idx=ivec.size()-1; //last possible element
		}
		const size_t elements_right = MAX(end_time_idx - index, end_elements);
		end = index + elements_right;

		//now, check (and modify) if the window could not be centered
		if(elements_left==elements_right) return true;
		if(!is_soft) return false;
		if(elements_left<elements_right) { //we hit the left border
			//get again the end of window
			size_t end_elements = (min_data_points>(elements_left+1))?min_data_points - (elements_left + 1):0;
			const Date end_date = ivec[start].date+min_time_span;
			size_t end_time_idx = (end_date>date)?IOUtils::seek(end_date, ivec, false):index; //end time criteria
			if(end_time_idx==IOUtils::npos) {
				end_time_idx=ivec.size()-1; //last possible element
			}
			const size_t elements_right = MAX(end_time_idx - index, end_elements);
			end = index + elements_right;
		} else { //we hit the right border
			//get again the start of window
			size_t start_elements = (min_data_points>(elements_right+1))?min_data_points - (elements_right + 1):0;
			const Date start_date = ivec[end].date-min_time_span;
			size_t start_time_idx = (start_date<date)?IOUtils::seek(start_date, ivec, false):index; //start time criteria
			if(start_time_idx==IOUtils::npos) {
				end_time_idx=0; //first possible element
			}
			const size_t elements_left = MAX(index - start_time_idx, start_elements);
			start = index - elements_left;
		}
	}


	return true;
}

} //namespace
