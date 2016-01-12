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
#include <meteoio/meteoFilters/FilterSnowNosnow.h>
#include <meteoio/meteoStats/libinterpol1D.h>

using namespace std;

namespace mio {

FilterSnowNosnow::FilterSnowNosnow(const std::vector<std::string>& vec_args, const std::string& name)
          : FilterBlock(name)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first;
}

void FilterSnowNosnow::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	const double TSS_offset = getTSSOffset(param, ivec);
	
	if (TSS_offset<IOUtils::C_TO_K(1.5)) {
		//TODO: TSS-TSG correlation and decide if TSS or TSG is used below
		Date prev_day;
		double TSS_daily_max, TSS_daily_min, TSS_daily_mean;
		for (size_t ii=0; ii<ovec.size(); ii++){ //now really perform the comparison
			double& value = ovec[ii](param);
			if (value==IOUtils::nodata) continue;
			
			int year, month, day;
			const Date day_start = getDailyStart(ivec[ii].date);
			if (day_start!=prev_day) {
				ovec[ii].date.getDate(year, month, day);
				getTSSDailyPpt(ivec, day_start, TSS_daily_min, TSS_daily_max, TSS_daily_mean);
				prev_day = day_start;
			}
			
			if (month>=7 && month<=12) { //early snow season
				if (TSS_daily_min>=IOUtils::C_TO_K(2.)) 
					value = 0.;
				else if (TSS_daily_min>=IOUtils::C_TO_K(0.)) {
					if (value<0.4) value = 0.;
				} else {
					if (TSS_daily_max>=IOUtils::C_TO_K(2.)) value = 0.;
				}
			} else { //late snow season
				if (TSS_daily_mean>IOUtils::C_TO_K(0.) && TSS_daily_max>=IOUtils::C_TO_K(5.)) value = 0.;
			}
		}
	} else {
		Date prev_day;
		double TSG_daily_var;
		for (size_t ii=0; ii<ovec.size(); ii++){ //now really perform the comparison
			double& value = ovec[ii](param);
			if (value==IOUtils::nodata) continue;
			
			const double TSG = ovec[ii](MeteoData::TSG);
			if (TSG==IOUtils::nodata) continue;
			
			//compute TSG daily variance
			const Date day_start = getDailyStart(ivec[ii].date);
			if (day_start!=prev_day) {
				TSG_daily_var = getDailyTSGVariance(ivec, day_start);
				prev_day = day_start;
			}
			
			if (TSG>IOUtils::C_TO_K(7.) || TSG_daily_var>1.) value = 0.;
		}
	}
}

//compute the TSS offset/correction
double FilterSnowNosnow::getTSSOffset(const unsigned int& param, const std::vector<MeteoData>& ivec) const 
{
	Date prev_day;
	double HS_daily_median, TSS_daily_median, RSWR_daily_10pc;
	bool high_tss_day;
	std::vector<double> tss_dat;
	
	for (size_t ii=0; ii<ivec.size(); ii++){
		if (ivec[ii](param)==IOUtils::nodata) continue;

		const Date day_start = getDailyStart(ivec[ii].date);
		if (day_start!=prev_day) {
			getDailyParameters(ivec, day_start, HS_daily_median, TSS_daily_median, RSWR_daily_10pc);
			high_tss_day = (HS_daily_median>0.3) && (TSS_daily_median>IOUtils::C_TO_K(1.)) && (RSWR_daily_10pc>350.);
			prev_day = day_start;
		}
		
		if (high_tss_day) tss_dat.push_back( ivec[ii](MeteoData::TSS) );
	}
	
	return Interpol1D::getMedian(tss_dat);
}

//daily values for TSS offset calc 
void FilterSnowNosnow::getDailyParameters(const std::vector<MeteoData>& ivec, const Date day_start, double &HS_daily_median, double &TSS_daily_median, double &RSWR_daily_10pc) const
{
	const Date day_end = day_start + 1;
	
	//extract values for the day
	std::vector<double> HS_dat, TSS_dat, RSWR_dat;
	for (size_t jj=0; jj<ivec.size(); jj++) {
		if (ivec[jj].date>=day_end) break;
		if (ivec[jj].date>=day_start) {
			HS_dat.push_back( ivec[jj](MeteoData::HS) );
			TSS_dat.push_back( ivec[jj](MeteoData::TSS) );
			RSWR_dat.push_back( ivec[jj](MeteoData::RSWR) );
		}
	}
	
	std::vector<double> quantiles;
	quantiles.push_back( 0.9 );
	std::vector<double> rswr_quantiles = Interpol1D::quantiles(RSWR_dat, quantiles);
	
	RSWR_daily_10pc = rswr_quantiles[0];
	HS_daily_median = Interpol1D::getMedian(HS_dat);
	TSS_daily_median = Interpol1D::getMedian(TSS_dat);
}

//daily values for TSS correction
void FilterSnowNosnow::getTSSDailyPpt(const std::vector<MeteoData>& ivec, const Date day_start, double &TSS_daily_min, double &TSS_daily_max, double &TSS_daily_mean) const
{
	const Date day_end = day_start + 1;
	
	//extract values for the day
	std::vector<double> TSS_dat;
	for (size_t jj=0; jj<ivec.size(); jj++) {
		if (ivec[jj].date>=day_end) break;
		if (ivec[jj].date>=day_start) {
			TSS_dat.push_back( ivec[jj](MeteoData::TSS) );
		}
	}
	
	TSS_daily_min = Interpol1D::min_element(TSS_dat);
	TSS_daily_max = Interpol1D::max_element(TSS_dat);
	TSS_daily_mean = Interpol1D::arithmeticMean(TSS_dat);
}

double FilterSnowNosnow::getDailyTSGVariance(const std::vector<MeteoData>& ivec, const Date day_start) const
{
	const Date day_end = day_start + 1;
	
	//extract values for the day
	std::vector<double> TSG_dat;
	for (size_t jj=0; jj<ivec.size(); jj++) {
		if (ivec[jj].date>=day_end) break;
		if (ivec[jj].date>=day_start) {
			TSG_dat.push_back( ivec[jj](MeteoData::TSG) );
		}
	}
	
	return Interpol1D::variance(TSG_dat);
}

/**
 * @brief For a given date, find the start of the day, considering that for midnight we return the day before!
 * (as is necessary for daily averages, sums, etc that can be provided at midnight for the day before)
 * @param resampling_date current date
 * @return start of the day or start of the day before in case of midnight
 */
Date FilterSnowNosnow::getDailyStart(const Date& resampling_date) const
{
	Date Start( resampling_date );
	Start.rnd(24*3600, Date::DOWN);
	if (Start==resampling_date) //if resampling_date=midnight GMT, the rounding lands on the exact same date
		Start -= 1.;
	
	return Start;
}


void FilterSnowNosnow::parse_args(std::vector<std::string> vec_args)
{
	if ( !vec_args.empty() ) //ie if there are arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);
}

} //end namespace
