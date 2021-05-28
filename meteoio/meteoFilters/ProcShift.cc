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
	
	//getShiftOffset(ivec);
	getPearson(ivec, Date(2019, 4, 1, 0, 0, +1.), 2.);
	return;
	
	/*for (size_t ii=0; ii<ovec.size(); ii++){
		//here, implement what has to be done on each data point
		//for example:
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		if (tmp < 0.){ //delete all values less than zero
			tmp = IOUtils::nodata;
		}
	}*/
}

//resample vecY to the exact same timestamps as vecX with linear interpolation within a given gap width in days
std::vector< std::pair<Date, double> > ProcShift::resampleVector(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY) const
{
	//if (vecX.size() != vecY.size()) throw IOException("Both vectors must have the same size", AT);
	std::vector< std::pair<Date, double> > vecResults;
	size_t jj = 0; //index in vecY
	
	for (size_t ii=0; ii<vecX.size(); ii++) {
		if (jj==vecY.size()) { //fill with nodata since there is nothing else to do to keep all vector the same size
			vecResults.push_back( make_pair(vecX[ii].first, IOUtils::nodata) );
			continue;
		}
		
		if (vecY[jj].first==vecX[ii].first) {
			vecResults.push_back( vecY[jj] );
			jj++;
		} else { //either the point is before / after or we need to interpolate
			bool not_found=false;
			while (vecY[jj].first<vecX[ii].first) {
				if (jj==(vecY.size()-1)) {
					vecResults.push_back( make_pair(vecX[ii].first, IOUtils::nodata) );
					not_found = true;
					break;
				}
				jj++;
			}
			while (vecY[jj].first>vecX[ii].first) {
				if (jj==0) {//vecY starts after vecX, skipping all points until they overlap
					vecResults.push_back( make_pair(vecX[ii].first, IOUtils::nodata) );
					not_found = true;
					break;
				}
				jj--;
			}
			
			if (not_found) continue;
			
			//now vecY[jj].first<=vecX[ii].first
			if (vecY[jj].first==vecX[ii].first) {
				vecResults.push_back( vecY[jj] );
				jj++;
			} else { //we need to interpolate between jj and jj+1
				if (jj==vecY.size()-1) return vecResults;
				const double x1 = vecY[jj].first.getJulian();
				const double x2 = vecY[jj+1].first.getJulian();
				const double y1 = vecY[jj].second;
				const double y2 = vecY[jj+1].second;
				if (x1==x2) throw IOException("Attempted division by zero", AT);
				
				if (y1!=IOUtils::nodata && y2!=IOUtils::nodata) {
					const double a = (y2 - y1) / (x2 - x1);
					const double b = y2 - a*x2;
					const double x = vecX[ii].first.getJulian();
					const double y = (a*x + b);
					
					vecResults.push_back( make_pair(vecX[ii].first, y) );
				} else {
					vecResults.push_back( make_pair(vecX[ii].first, IOUtils::nodata) );
				}
			}
		}
	}
	
	return vecResults;
}

//Pearson's coefficient between two datasets, X and Y which for us are extracted from the same ivec but for different parameters
//the sums for X are recomputed to make sure that if Y[ii] is nodata, no value is taken for the matching X
double ProcShift::getPearson(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t &ii_start, const size_t &ii_end) const
{
	if (vecX.size() != vecY.size()) throw IOException("Both vectors must have the same size, "+IOUtils::toString(vecX.size())+" != "+IOUtils::toString(vecY.size()), AT);
	size_t count=0;
	double sumX=0., sumX2=0., sumY=0., sumY2=0., sumXY=0.;
	for (size_t ii=ii_start; ii<ii_end; ii++) {
		const double valueX = vecX[ii].second;
		const double valueY = vecY[ii].second;
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

double ProcShift::getPearson(const std::vector<MeteoData>& ivec, const Date& dt_start, const double& width_d) const
{
	if (ivec.empty()) return IOUtils::nodata;
	
	const int nrSteps = 24*6;
	const double step = width_d / static_cast<double>(nrSteps);
	
	//compute the reference vector vecX
	//find the indices to cover the right extend for TA_1
	size_t ii_Xstart=0;
	while (ii_Xstart<ivec.size() && ivec[ii_Xstart].date<dt_start) ii_Xstart++;
	if (ii_Xstart>0) ii_Xstart--; //to include the last timestamp before the period
	size_t ii_Xend=ii_Xstart;
	while (ii_Xend<ivec.size() && ivec[ii_Xend].date<(dt_start+width_d)) ii_Xend++;
	if (ii_Xend<ivec.size()-1) ii_Xend++; //to include the first timestamp after the period
	std::vector< std::pair<Date, double> > vecX(ii_Xend-ii_Xstart+1);
	for (size_t ii=0; ii<vecX.size(); ii++) vecX[ii] = make_pair(ivec[ii+ii_Xstart].date, ivec[ii+ii_Xstart]("TA_1"));
	
	//try nrSteps offsets, from -width_d/2 to +width_d/2
	//for each, compute the average Pearson's coefficient
	for (int ss=0; ss<nrSteps; ss++) {
		const double offset = static_cast<double>(ss)*step - 0.5*width_d;
		
		//compute the comparison vector vecY
		//find the indices to cover the right extend for TA_2
		size_t ii_Ystart=0;
		while (ii_Ystart<ivec.size() && ivec[ii_Ystart].date<(dt_start+offset)) ii_Ystart++;
		if (ii_Ystart>0) ii_Ystart--; //to include the last timestamp before the period
		size_t ii_Yend=ii_Ystart;
		while (ii_Yend<ivec.size() && ivec[ii_Yend].date<(dt_start+offset+width_d)) ii_Yend++;
		if (ii_Yend<ivec.size()-1) ii_Yend++; //to include the first timestamp after the period
		std::vector< std::pair<Date, double> > vecY(ii_Yend-ii_Ystart+1);
		for (size_t ii=0; ii<vecY.size(); ii++) vecY[ii] = make_pair(ivec[ii+ii_Ystart].date + offset, ivec[ii+ii_Ystart]("TA_2"));
		
		std::vector< std::pair<Date, double> > vecResults( resampleVector(vecX, vecY) );
		
		double sum=0.;
		size_t count=0;
		for (size_t ii=0; ii<vecResults.size(); ii++) {
			const Date dt_middle( 0.5*(vecX.front().first.getJulian() + vecX.back().first.getJulian()), vecX.front().first.getTimeZone() );
			const double Pearson = getPearson(vecX, vecResults, 0, vecX.size());
			if (Pearson!=IOUtils::nodata) {
				sum += Pearson;
				count++;
			}
		}
		
		if (count>0)
			std::cout << offset*24.*60. << " " << sum/static_cast<double>(count) << "\n";
			//std::cout << "Offset=" << offset*24.*60. << " Pearson_avg=" << sum/static_cast<double>(count) << "\n";
	}
}

double ProcShift::getShiftOffset(const std::vector<MeteoData>& ivec) const
{
	if (ivec.empty()) return IOUtils::nodata;
	static const double width_d = 1.; //in days
	const size_t n = ivec.size();
	
	//build vecX and vecY
	std::vector< std::pair<Date, double> > vecX(n), vecY(n);
	for (size_t ii=0; ii<n; ii++) vecX[ii] = make_pair(ivec[ii].date, ivec[ii]("TA_1"));
	for (size_t ii=0; ii<n; ii++) vecY[ii] = make_pair(ivec[ii].date, ivec[ii]("TA_2"));
	
	std::cout << "Pearson's coefficients\n"; 
	std::vector< std::pair<Date, double> > vecResults( resampleVector(vecX, vecY) );
	for (size_t ii=0; ii<n; ii++) {
		const Date dt_end( vecX[ii].first+width_d );
		size_t ii_end = ii;
		while (ii_end<n && vecX[ii_end].first<dt_end) ii_end++;
		if (ii_end==n) break;
		const Date dt_middle( 0.5*(vecX[ii].first.getJulian() + vecX[ii_end].first.getJulian()), vecX[ii].first.getTimeZone() );
		const double Pearson = getPearson(vecX, vecResults, ii, ii_end);
		if (Pearson!=IOUtils::nodata) 
			std::cout << dt_middle.toString(Date::ISO) << " " << fabs(Pearson) << "\n";
		else
			std::cout << dt_middle.toString(Date::ISO) << " " << IOUtils::nodata << "\n";
	}
	
	return IOUtils::nodata;
}

void ProcShift::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "Filters::"+block_name );
	//for a filter that does not take any arguments
	if ( !vecArgs.empty() ) //ie if there are arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for "+where, AT);
}

} //end namespace
