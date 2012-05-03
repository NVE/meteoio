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
			for (size_t kk=0; kk<min_data_points; kk++){
				if (ivec.size() > kk) {
					elements_right++;
					vec_window.push_back(&ivec[kk]);
				}
			}

			if (elements_right > 0) elements_left = 1;
		} else if ((centering == WindowedFilter::left) || (centering == WindowedFilter::center)){
			if (ivec.size() > 0) {
				elements_left = elements_right = 1;
				vec_window.push_back(&ivec[0]);
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
		}
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

void WindowedFilter::get_window_fast(const unsigned int& index, const unsigned int& ivec_size,
                                     unsigned int& index_start, unsigned int& index_end)
{
	if ((ivec_size == 0) || (ivec_size <= index)){
		index_end = 0;
		index_start = 1;
	}

	if ((index == 0) || (last_index > index)){ //reset global variables
		elements_left = elements_right = 0;

		if ((centering == WindowedFilter::right) || (is_soft)){
			for (unsigned int kk=0; kk<min_data_points; kk++){
				if (ivec_size > kk) elements_right++;
			}

			if (elements_right > 0) elements_left = 1;
		} else if ((centering == WindowedFilter::left) || (centering == WindowedFilter::center)){
			if (ivec_size > 0) {
				elements_left = elements_right = 1;
			}
		}

		last_index = 0;

		startIndex = index_start = index + 1 - elements_left;
		endIndex = index_end = index + elements_right - 1;

		return;
	}

	if (index != (last_index+1))
		throw IOException("get_window function only to be used with increments of 1 for the index", AT);

	unsigned int window_size = endIndex - startIndex + 1;

	//check whether a window is available
	if (is_soft){
		if (startIndex <= endIndex){
			if (centering == WindowedFilter::right){
				//Try to move right, if it doesn't work, don't change anything
				if (ivec_size > (index + window_size - 1)) { //shift window one point to the right
					endIndex++;
					startIndex++;
				}
			} else if (centering == WindowedFilter::left){
				if (index >= window_size){
					if (ivec_size > index){ //otherwise don't touch the whole thing
						startIndex++;
						endIndex++;
					}
				} else {
					elements_left++;
					elements_right--;
				}
			} else if (centering == WindowedFilter::center){
				if (elements_right <= elements_left){
					if (ivec_size > (index+elements_right-1)){ //otherwise don't touch the whole thing
						endIndex = index+elements_right-1;
						startIndex++;
					} else {
						elements_right--;
						elements_left++;
					}
				} else {
					elements_right--;
					elements_left++;
				}
			}
		}
		index_start = startIndex;
		index_end = endIndex;
	} else { //!is_soft
		if (centering == WindowedFilter::right){
			if (elements_right >= min_data_points){
				startIndex++;
				elements_right--;

				if (ivec_size > (index+min_data_points-1)) { //shift window one point to the right
					endIndex++;
					elements_right++; //elements_left will stay at a constant 1
				}
			}
		} else if (centering == WindowedFilter::left){
			if (elements_left >= min_data_points){
				startIndex++;
				elements_left--;

				if (ivec_size > index) { //shift window one point to the right
					endIndex++;
					elements_left++; //elements_left will stay at a constant 1
				}
			} else {
				if (ivec_size > index) { //broaden window
					endIndex++;
					elements_left++;
				}
			}
		} else if (centering == WindowedFilter::center){
			if ((elements_left + elements_right - 1) >= min_data_points){
				startIndex++;
				if (elements_right > 0) elements_right--;

				if (ivec_size > (index+elements_left-1)) { //shift window one point to the right
					endIndex++;
					elements_right++; //elements_left will stay at a constant
				}
			} else {
				if (ivec_size > (index+elements_left-1)) { //shift window one point to the right
					endIndex++;
					elements_left++;
				}

				if ((elements_left + elements_right - 1) < min_data_points){
					if (ivec_size > (index+elements_left-1)){ //shift window one point to the right
						endIndex++;
						elements_right++;
					}
				}
			}
		}

		if ((elements_left + elements_right - 1) >= min_data_points){
			index_start = startIndex;
			index_end = endIndex;
		} else {
			index_end = 0;
			index_start = 1;
		}
	}

	last_index = index;
}

void WindowedFilter::get_window(const unsigned int& index, const unsigned int& ivec_size,
                                unsigned int& index_start, unsigned int& index_end)
{
	if ((centering == WindowedFilter::right)){
		index_end   = index + min_data_points - 1;
		index_start = index;

		if (index_end >= ivec_size){
			if (is_soft){
				index_end = ivec_size - 1;
				if (ivec_size >= min_data_points){
					index_start = ivec_size + 1 - min_data_points;
				} else {
					index_start = 0;
				}
			} else {
				index_start = index_end + 1;
			}
		}

		return;
	}

	if ((centering == WindowedFilter::left)){
		index_end = index;

		if (index <= min_data_points){
			if (is_soft){
				index_start = 0;
				index_end = MIN(min_data_points - 1, ivec_size - 1);
			} else {
				index_start = index_end + 1;
			}
		} else {
			index_start = index + 1 - min_data_points;
		}

		return;
	}

	if ((centering == WindowedFilter::center)){
		unsigned int el_left = min_data_points / 2;
		index_start = index_end = index;

		//first calc index_start
		if (el_left > index){
			if (is_soft){
				index_start = 0;
			} else {
				index_start = index_end + 1;
				return;
			}
		} else {
			index_start -= el_left;
		}

		if ((index_end - index_start + 1) == min_data_points) return; //Nothing more to do

		index_end = index_start + min_data_points - 1;
		if (index_end >= ivec_size){
			if (is_soft){
				index_end = ivec_size - 1;
				while ((index_start != 0) && ((index_end - index_start + 1) < min_data_points)){
					index_start--;
				}
			} else {
				index_start = index_end + 1;
				return;
			}
		}
	}
}

} //namespace
