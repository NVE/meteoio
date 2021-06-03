// SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <meteoio/meteoFilters/ProcShift.h>
#include <meteoio/meteoStats/libinterpol1D.h>

using namespace std;

namespace mio {

ProcShift::ProcShift(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg)
          : ProcessingBlock(vecArgs, name, cfg) //this has to match the class you are inheriting from! ie ProcessingBlock or WindowedFilter
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}

void ProcShift::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ivec.empty()) return;
	
	const size_t param_ta1 = ivec[0].getParameterIndex("TA_1");
	const size_t param_ta2 = ivec[0].getParameterIndex("TA_2");
	const double sampling_rate = getMedianSampling(param_ta1, ivec);
	const std::vector< std::pair<Date, double> > vecX( resampleVector(ivec, param_ta1, sampling_rate) );
	const std::vector< std::pair<Date, double> > vecY( resampleVector(ivec, param_ta2, sampling_rate) );
	
// 	//const Date dt_offset(2019, 4, 1, 0, 0, +1.);
// 	const Date dt_offset(2019, 1, 29, 4, 0, +1.);
// 	const double width_d = 2.;
// 	static const size_t width_idx = static_cast<size_t>( round( width_d / sampling_rate ) );
// 	size_t ii=0;
// 	for (; ii<ivec.size(); ii++) if (vecX[ii].first>=dt_offset) break;
// 	getOffset(vecX, vecY, ii, width_idx);
// 	return;
	
	/*std::cout << "ta1\n";
	for (size_t ii=0; ii<vecX.size(); ii++) {
		std::cout << "@ " << vecX[ii].first.toString(Date::ISO) << " " << vecX[ii].second << "\n";
	}
	std::cout << "\n";
	exit;*/
	
	static const double width_d = 3.;
	static const size_t width_idx = static_cast<size_t>( round( width_d / sampling_rate ) );
	for (size_t ii=0; ii<ivec.size(); ii++) {
		//if (vecX[ii].second == IOUtils::nodata) continue; //preserve nodata values

		getOffset(vecX, vecY, ii, width_idx);
	}
	
	
	
	/*for (size_t ii=0; ii<ovec.size(); ii++) {
		//here, implement what has to be done on each data point
		//for example:
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		if (tmp < 0.){ //delete all values less than zero
			tmp = IOUtils::nodata;
		}
	}*/
}

double ProcShift::getMedianSampling(const size_t& param, const std::vector<MeteoData>& ivec) const
{
	std::vector<double> vecSampling;
	
	for (size_t ii=0; ii<(ivec.size()-1); ii++){
		if (ivec[ii](param)!=IOUtils::nodata && ivec[ii+1](param)!=IOUtils::nodata)
			vecSampling.push_back( ivec[ii+1].date.getJulian(true) - ivec[ii].date.getJulian(true) );
	}
	
	return Interpol1D::getMedian(vecSampling, false);
}

std::vector< std::pair<Date, double> > ProcShift::resampleVector(const std::vector<MeteoData>& ivec, const size_t& param, const double& sampling_rate) const
{
	const size_t n_ivec = ivec.size();
	const Date dt_start( ivec.front().date );
	const size_t nrSteps = static_cast<size_t>(round( (ivec.back().date.getJulian(true) - dt_start.getJulian(true)) / sampling_rate ));
	std::vector< std::pair<Date, double> > vecResults;
	
	size_t jj=0; //position within ivec
	for (size_t ii=0; ii<nrSteps; ii++) {
		const Date dt( dt_start+(static_cast<double>(ii)*sampling_rate) );
		
		if (ivec[jj].date==dt) {
			vecResults.push_back( make_pair(dt, ivec[jj](param)) );
			jj++;
			continue;
		}
		
		//find the first element >= dt
		while (ivec[jj].date<dt) {
			if (jj==(n_ivec-1)) {
				return vecResults;
			}
			jj++;
		}
		
		if (jj>0) {
			if (jj==n_ivec) return vecResults;
			const double x1 = ivec[jj-1].date.getJulian(true);
			const double x2 = ivec[jj].date.getJulian(true);
			if ((x2-x1) > 2.*sampling_rate) { //only interpolate between nearby points
				vecResults.push_back( make_pair(dt, IOUtils::nodata) );
				continue;
			}
			
			const double y1 = ivec[jj-1](param);
			const double y2 = ivec[jj](param);
			if (x1==x2) throw IOException("Attempted division by zero", AT);
			
			if (y1!=IOUtils::nodata && y2!=IOUtils::nodata) {
				const double a = (y2 - y1) / (x2 - x1);
				const double b = y2 - a*x2;
				const double x = dt.getJulian(true);
				const double y = (a*x + b);
				
				vecResults.push_back( make_pair(dt, y) );
			} else {
				vecResults.push_back( make_pair(dt, IOUtils::nodata) );
			}
		} else {
			vecResults.push_back( make_pair(dt, IOUtils::nodata) );
		}
	}
}

//Pearson's coefficient between two datasets
//the sums for X are recomputed to make sure that if Y[ii] is nodata, no value is taken for the matching X
//indices given for vecX
double ProcShift::getPearson(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t& curr_idx, const size_t& width_idx, const int& offset) const
{
	//compute the required data window, accounting for offsets in vecY that could 
	//bring us outside vecY
	size_t startIdx = curr_idx-width_idx/2;
	if (-offset>(signed)startIdx) startIdx = static_cast<size_t>(-offset);
	size_t endIdx = curr_idx+width_idx/2;
	if (endIdx+static_cast<size_t>(offset)>vecX.size()) endIdx = vecX.size() - static_cast<size_t>(offset);
	
	//std::cout << "curr_idx=" << curr_idx << " width_idx=" << width_idx << " offset=" << offset << " startIdx=" << startIdx << " endIdx=" << endIdx << "\n";
	size_t count=0;
	double sumX=0., sumX2=0., sumY=0., sumY2=0., sumXY=0.;
	for (size_t ii=startIdx; ii<endIdx; ii++) {
		const double valueX = vecX[ii].second;
		const double valueY = vecY[static_cast<size_t>((signed)ii+offset)].second;
		if (valueX!=IOUtils::nodata && valueY!=IOUtils::nodata) {
			sumX += valueX;
			sumX2 += valueX * valueX;
			
			sumY += valueY;
			sumY2 += valueY * valueY;
			
			sumXY += valueX * valueY;
			count++;
		}
	}
	
	const double Xcontribution = static_cast<double>(count)*sumX2 - sumX*sumX;
	const double Ycontribution = static_cast<double>(count)*sumY2 - sumY*sumY;
	if (Xcontribution<=0. || Ycontribution<=0.) return IOUtils::nodata; //either count==0, count==1 or all values are identical
	
	const double pearson = (static_cast<double>(count)*sumXY - sumX*sumY) / (sqrt(Xcontribution) * sqrt(Ycontribution));
	return pearson;
}

int ProcShift::getOffsetFullScan(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t& curr_idx, const size_t& width_idx, const int& range_min, const int& range_max) const
{
	int offset_max = IOUtils::inodata;
	double pearson_max = IOUtils::nodata;
	unsigned int count = 0;
	for (int offset=range_min; offset<=range_max; offset++) {
		const double pearson = getPearson(vecX, vecY, curr_idx, width_idx, offset);
		if (pearson==IOUtils::nodata) continue;
		
		//std::cout << offset << " " << pearson << "\n";
		count++;
		if (pearson>pearson_max) {
			offset_max = offset;
			pearson_max = pearson;
		}
	}
	
	if (pearson_max!=IOUtils::inodata)
		std::cout << vecX[curr_idx].first.toString(Date::ISO) << " " << offset_max << " " << pearson_max << "\n";
	
	//we want to avoid looking for a maximum among too few points
	if (count<10 || pearson_max==IOUtils::nodata)
		return IOUtils::inodata;
	
	return offset_max;
}

double ProcShift::getOffset(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t& curr_idx, const size_t& width_idx) const
{
	//https://en.wikipedia.org/wiki/Golden-section_search
	static const double r = Cst::phi - 1.;
	static const double c = 1. - r;
	
	/*int range_min = -static_cast<int>(width_idx)/4;
	int range_max = static_cast<int>(width_idx)/4;*/
	int range_min = -12*6; //12 hours, times 10 minutes sampling rate
	int range_max = 12*6;
	int offset_c = IOUtils::inodata, offset_r = IOUtils::inodata;
	double pearson_c = IOUtils::nodata, pearson_r = IOUtils::nodata;
	size_t count = 0;
	
	//offset_c = getOffsetFullScan(vecX, vecY, curr_idx, width_idx, range_min, range_max);
	
	while (true) {
		if (offset_c==IOUtils::inodata) {
			offset_c = static_cast<int>( round(((double)range_max - (double)range_min)*c) ) + range_min;
			pearson_c = getPearson(vecX, vecY, curr_idx, width_idx, offset_c);
		}
		
		if (offset_r==IOUtils::inodata) {
			offset_r = static_cast<int>( round(((double)range_max - (double)range_min)*r) ) + range_min;
			pearson_r = getPearson(vecX, vecY, curr_idx, width_idx, offset_r);
		}
		
		if (pearson_c==IOUtils::nodata || pearson_r==IOUtils::nodata) { //we must do a full scan
			std::cout << "must do a full scan, curr_idx=" << curr_idx << "\n";
			offset_c = getOffsetFullScan(vecX, vecY, curr_idx, width_idx, range_min, range_max);
			break;
		}
		
		if (pearson_c > pearson_r) {
			range_max = offset_r;
			offset_r = offset_c;
			pearson_r = pearson_c;
			offset_c = IOUtils::inodata;
		} else {
			range_min = offset_c;
			offset_c = offset_r;
			pearson_c = pearson_r;
			offset_r = IOUtils::inodata;
		}
		
		//check convergence
		if ((range_max-range_min) < 2) 
			break; //convergence criteria: within 1 cell of optimum
		if (count>100) {
			std::cout << "No convergence at " << vecX[curr_idx].first.toString(Date::ISO) << " offset_c=" << offset_c*24.*60. << " offset_r=" << offset_r*24.*60. << "\n";
			break;
		}
		
		count++;
	}
	
	if (offset_c!=IOUtils::inodata)
		std::cout << vecX[curr_idx].first.toString(Date::ISO) << " " << offset_c << " " << pearson_c << "\n";
	
	return offset_c;
}

void ProcShift::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	//for a filter that does not take any arguments
	if ( !vecArgs.empty() ) //ie if there are arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for "+where, AT);
}

} //end namespace
