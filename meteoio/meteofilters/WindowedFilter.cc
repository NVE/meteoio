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
	: FilterBlock(name), is_soft(false), min_data_points(1), min_time_span(0.0), 
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
 * @param vec_window A vector of pointers to MeteoData
 */
const std::vector<const MeteoData*>& WindowedFilter::get_window(const unsigned int& index, 
                                                          const std::vector<MeteoData>& ivec)
{
	//cout << "Requesting index " << index << endl;

	if ((index == 0) || (last_index > index)){ //reset global variables
		vec_window.clear();
		elements_left = elements_right = 0;

		if ((centering == WindowedFilter::right) || (is_soft)){	
			for (unsigned int kk=0; kk<min_data_points; kk++){
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

		//cout << "Init: " << elements_left << "/" << elements_right << endl;
		
		last_index = 0;

		return vec_window;
	}

	if (index != (last_index+1))
		throw IOException("get_window function only to be used with increments of 1 for the index", AT);

	//check whether a window is available
	if (is_soft){
		if (vec_window.size() > 0){
			if (centering == WindowedFilter::right){
				if (vec_window.size() > 0){
					//Try to move right, if it doesn't work, don't change anything
					if (ivec.size() > (index+vec_window.size()-1)) { //shift window one point to the right
						vec_window.push_back(&ivec[index+vec_window.size()-1]);
						vec_window.erase(vec_window.begin());
					}
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

}
