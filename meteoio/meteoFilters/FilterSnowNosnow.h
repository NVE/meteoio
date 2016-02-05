/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __FilterSnowNosnow_H__
#define __FilterSnowNosnow_H__

//#include <meteoio/meteoFilters/WindowedFilter.h> //use this one for filters relying on a data window, for example std_dev
#include <meteoio/meteoFilters/FilterBlock.h> //use this one for all others

#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterSnowNosnow
 * @ingroup processing
 * @author Anna-Maria Tilg and Mathias Bavay
 * @date   2015-12-16
 * @brief This filter is used to distinguish if snow (HS) is on the ground or not, because the 
 * ultrasonic sensor cannot distinguish between snow or vegetation/grass on the ground. 
 * The filter is based on total snow depth (HS), snow surface temperature (TSS), ground surface temperature (TSG)
 * and reflected shortwave radiation (RSWR). 
 * Different steps to do: 
 * #step 1#: calculate possible offset of TSS (raison: at some stations in some springs the TSS increases 
 * although snow is still on the ground) 
 * #step 2#: calculate correlation of TSS and TSG in spring (normally both temperatures increase at the same 
 * time in spring which results in a high correlation; low correlation if TSS and TSG increase not parallel 
 * which leads in connection with a high offset of TSS to the assumption, that TSS is false) 
 * #step 3a#: if TSS has a low offset and the correlation between TSS and TSG is high, the algorithm 
 * analyses based on the daily Max/Min/Mean of TSS if snow is on the ground or not 
 * #step 3b#: if no TSS is available or the offset of TSS is high/correlation of TSS and TSG is low, the 
 * algorithm analyses based on the variance of TSG and the value of TSG if snow is on the ground or 
 * not. 
 * References/Literature: Tilg, A.-M., Marty C. and G. Klein, 2015: An automatic algorithm for validating snow 
 * depth measurements of IMIS stations. Swiss Geoscience Meeting 2015    
 * Remarks:
 * - one of the used criteria is currently only valid for mid-northern latitudes;
 * - nodata values are excluded ?????? => checken  
 * - Two arguments expected (both have to be fulfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - only window position "center" possible
 * - keyword "soft" not allowed 
 *
 * @code
 * HS::filter1	= FilterSnowNosnow
 * @endcode
 */

class FilterSnowNosnow : public FilterBlock { //use this one for simple filter that only look at one data point at a time, for example min_max
//class FilterSnowNosnow : public ProcessingBlock { //use this one for data corrections
//class FilterSnowNosnow : public WindowedFilter { //use this one for filters relying on a data window, for example std_dev
	public:
		FilterSnowNosnow(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void filterOnTsg(const unsigned int& param, const size_t& ii, std::vector<MeteoData>& ovec);
		void filterOnTss(const unsigned int& param, const size_t& ii, const double& tss_offset, std::vector<MeteoData>& ovec);
		void parse_args(std::vector<std::string> vec_args);
		
		static double getTssTsgCorrelation(const std::vector<MeteoData>& ovec, const size_t& firstWarmDay_idx);
		static void findFirstWarmDay(const std::vector<MeteoData>& ovec, size_t &tssWarmDay_idx, size_t &tsgWarmDay_idx);
		static double getTSSOffset(const unsigned int& param, const std::vector<MeteoData>& ivec);
		static void getDailyParameters(const std::vector<MeteoData>& ivec, const Date day_start, double &HS_daily_median, double &TSS_daily_median, double &RSWR_daily_10pc);
		static void getTSSDailyPpt(const std::vector<MeteoData>& ivec, const Date day_start, double &TSS_daily_min, double &TSS_daily_max, double &TSS_daily_mean);
		static double getDailyTSGVariance(const std::vector<MeteoData>& ivec, const Date day_start);
		static Date getDailyStart(const Date& resampling_date);
		
		Date prev_day;
		double TSS_daily_max, TSS_daily_min, TSS_daily_mean, TSG_daily_var;
		int month;
};

} //end namespace

#endif
