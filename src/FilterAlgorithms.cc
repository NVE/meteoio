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
#include "FilterAlgorithms.h"

 /**
 * @page filters Filters overview
 * The filtering infrastructure is described in FilterAlgorithms (for its API). The goal of this page is to give an overview of the available filters and their usage.
 * 
 * @section filters_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * @section filters_available Available filters
 * The filters that are currently available are the following:
 * - rate: calculates the derivative at the current point and reject the point if the derivative is above a given threshold.
 * - min_max: filter all values outside of the interval [min, max], soft operation supported
 * - min: filter out all values less than the given parameter, soft operation supported
 * - max: filter out all values greater than the given parameter, soft operation supported
 *
 * A few data transformations are also supported besides filtering:
 * - resample: resamples (if necessary) the data so that if a required time stamp is between two data points, an interpolated value will be returned (instead of nodata).
 * - accumulation: accumulates the data on a given period. A practical use is to return hourly precipitations from a sensor measuring precipitation on a 10 minutes interval.
 */

using namespace std;

#define PI 3.141592653589

std::map<std::string, FilterProperties> FilterAlgorithms::filterMap;
const bool FilterAlgorithms::__init = FilterAlgorithms::initStaticData();

bool FilterAlgorithms::initStaticData()
{
	filterMap["rate"]     = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::RateFilter);
	filterMap["resample"] = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::ResamplingFilter);
	filterMap["median_avg"] = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MedianAvgFilter);
	filterMap["mean_avg"] = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MeanAvgFilter);
	filterMap["wind_avg"] = FilterProperties(false, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::WindAvgFilter);
	filterMap["min_max"]  = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MinMaxFilter);
	filterMap["min"]      = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MinValueFilter);
	filterMap["max"]      = FilterProperties(true, (unsigned int)1, Date_IO(0.0), &FilterAlgorithms::MaxValueFilter);

	return true;
}

const FilterProperties& FilterAlgorithms::filterProperties(const std::string& filtername)
{
	std::map<std::string, FilterProperties>::const_iterator it;
	it = filterMap.find(filtername);

	if (it==filterMap.end())
		throw UnknownValueException("Unknown filter called: " + filtername, AT);

	return it->second;
}

void FilterAlgorithms::parseFilterArguments(const std::string& filtername, const std::vector<std::string>& vecArgs_in,
							    const unsigned int& minArgs, const unsigned int& maxArgs, 
							    bool& isSoft, std::vector<double>& vecArgs_out)
{
	isSoft = false;
	unsigned int argindex = 0;
	vecArgs_out.clear(); 
	double tmp;

	//Test for softness
	if (vecArgs_in.size() > 0)
		if (vecArgs_in[0] == "soft"){ isSoft = true; argindex=1;}

	try {
		for (unsigned int ii=argindex; ii<vecArgs_in.size(); ii++){
			IOUtils::convertString(tmp, vecArgs_in[ii]);
			vecArgs_out.push_back(tmp);
		}

		if ((vecArgs_out.size() < minArgs) || (vecArgs_out.size() > maxArgs))
			throw InvalidArgumentException("Wrong number of arguments for filter " + filtername, AT);
	} catch(std::exception& e){
		std::cerr << "[E] While processing arguments for filter " << filtername << std::endl;
		throw;
	}
}

void FilterAlgorithms::parseWindowFilterArguments(const std::string& filtername, 
							const std::vector<std::string>& vecArgs_in,
							const unsigned int& minArgs, const unsigned int& maxArgs,
							bool& isSoft, std::string& windowposition, std::vector<double>& vecArgs_out)
{
	std::vector<std::string> vecArgs_new;
	bool alreadyDefined = false;

	for (unsigned int ii=0; ii<vecArgs_in.size(); ii++){
		//Check for "left", "right" or "center" keywords
		if ((vecArgs_in[ii] == "center") || (vecArgs_in[ii] == "left") || (vecArgs_in[ii] == "right")){
			if (alreadyDefined) 
				throw InvalidArgumentException("Invalid or redefined arguments for filter " + filtername, AT);
			windowposition = vecArgs_in[ii];
			alreadyDefined = true;
		} else {
			vecArgs_new.push_back(vecArgs_in[ii]);
		}
	}

	parseFilterArguments(filtername, vecArgs_new, minArgs, maxArgs, isSoft, vecArgs_out);
}


/******************************************************************************
 * The following functions are implementations of different filter algorithms *
 ******************************************************************************/

bool FilterAlgorithms::RateFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						    const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						    const unsigned int& paramindex,
						    std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)date; (void)vecS; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
	parseFilterArguments("rate", _vecArgs, 1, 1, isSoft, vecArgs);
	
	Date_IO deltatime(vecArgs[0]/1440.0); //making a julian date out of the argument given in minutes

	int startposition = pos;
	for (int ii=vecFilteredM.size()-1; ii>=0; ii--){
		unsigned int mypos = startposition;

		double sum = 0.0;
		while((vecM[mypos].date + deltatime) > vecM[startposition].date){
			const double& val = vecM[mypos].param(paramindex);
			if (val != IOUtils::nodata)
				sum += vecM[mypos].param(paramindex);

			//cout << vecM[mypos].date << " HNW:" << vecM[mypos].param(paramindex) << "  SUM:" << sum << endl;		

			if (mypos>0) mypos--;
			else break;
		}
		
		//cout << "Accumulation for element (" << ii << "): " << sum << endl;
		vecFilteredM[ii].param(paramindex) = sum;

		if (startposition>0) startposition--; 
		else break;
	}

	return true;
}

bool FilterAlgorithms::ResamplingFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
							const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
							const unsigned int& paramindex,
							std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)_vecArgs;
	if ((vecFilteredM.size()==1) &&(date==vecFilteredM[0].date)){//Nothing to do
		return false; //Interpretation: filter not applied
	} else if ((vecFilteredM.size() < 2) || (vecFilteredM.size() > 3)){
		throw IOException("Not enough data to do resampling or index out of bounds", AT);
	} else if (vecFilteredM.size()==2){
		//add another element
		MeteoData newmd;
		newmd.date = date;
		newmd.setResampled(true);
		vecFilteredM.insert(vecFilteredM.begin() + 1, newmd);
		vecFilteredS.insert(vecFilteredS.begin() + 1, vecFilteredS[0]);
	}

	const MeteoData& tmpmd1 = vecFilteredM[0];
	MeteoData& resampledmd  = vecFilteredM[1];
	const MeteoData& tmpmd2 = vecFilteredM[2];

	double weight = (date.getJulian() - tmpmd1.date.getJulian()) / (tmpmd2.date.getJulian() - tmpmd1.date.getJulian());

	const double& val1 = tmpmd1.param(paramindex);
	const double& val2 = tmpmd2.param(paramindex);
	
	if ((val1 == IOUtils::nodata) || (val2 == IOUtils::nodata)){
		resampledmd.param(paramindex) = IOUtils::nodata;
	} else {
		resampledmd.param(paramindex) = Interpol1D::linearInterpolation(val1, val2, weight);
	}

	return true;
}

bool FilterAlgorithms::MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						 const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						 const unsigned int& paramindex,
						 std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
	parseFilterArguments("min_max", _vecArgs, 2, 2, isSoft, vecArgs);

	sort(vecArgs.begin(), vecArgs.end());

	//Run actual MinMax filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);

		if (value == IOUtils::nodata) continue;

		if (value<vecArgs[0]){
			if (isSoft) value=vecArgs[0];
			else value=IOUtils::nodata;				
			//cout << "Changed: " << value << endl;
		}
		if (value>vecArgs[1]){
			if (isSoft) value=vecArgs[1];
			else value=IOUtils::nodata;			
			//cout << "Changed: " << value << endl;	
		}
	}

	return true;
}

bool FilterAlgorithms::MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						   const unsigned int& paramindex,
						   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
	parseFilterArguments("min", _vecArgs, 1, 1, isSoft, vecArgs);

	//Run actual MinValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);
		if (value == IOUtils::nodata) continue;
		if (value<vecArgs[0]){
			if (isSoft) value=vecArgs[0];
			else value=IOUtils::nodata;				
		}
	}

	return true;
}

bool FilterAlgorithms::MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
						   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
						   const unsigned int& paramindex,
						   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs; 
	parseFilterArguments("max", _vecArgs, 1, 1, isSoft, vecArgs);

	//Run actual MaxValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);
		if (value == IOUtils::nodata) continue;
		if (value>vecArgs[0]){
			if (isSoft) value=vecArgs[0];
			else value=IOUtils::nodata;				
		}
	}

	return true;
}

bool FilterAlgorithms::MedianAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex,
				   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecS; (void)vecFilteredS;
	//Remarks:
	//1) nodata values are not considered when calculating the median
	//2) if there is an even number of window elements the arithmetic mean is used to calculate the median
	//3) Two arguments expected (both have to be fullfilled for the filter to start operating): 
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//4) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//5) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
	//  

	std::vector<double> vecWindow;
	if (!getWindowData("median_avg", vecM, pos, date, _vecArgs, paramindex, vecWindow))
		return false; //Not enough data to meet user configuration

	//Calculate median
	double median=IOUtils::nodata;
	unsigned int vecSize = vecWindow.size();
	unsigned int middle = (unsigned int)(vecSize/2);
	sort(vecWindow.begin(), vecWindow.end());

	if (vecSize == 0){
		return false; //only nodata values detected or other problem
	} else if ((vecSize % 2) == 1){ //uneven
		median = vecWindow.at(middle);
	} else { //use arithmetic mean of element n/2 and n/2-1
		median = Interpol1D::linearInterpolation(vecWindow.at(middle-1), vecWindow.at(middle), 0.5);
	}

	//cout << "Median value: " << median << endl;
	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		vecFilteredM[ii].param(paramindex) = median;
	}

	return true;
}

bool FilterAlgorithms::MeanAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex,
				   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecS; (void)vecFilteredS;
	//Remarks:
	//1) nodata values are not considered when calculating the median
	//2) Two arguments expected (both have to be fullfilled for the filter to start operating): 
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//3) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//4) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
	//  
	std::vector<double> vecWindow;
	if (!getWindowData("median_avg", vecM, pos, date, _vecArgs, paramindex, vecWindow))
		return false; //Not enough data to meet user configuration

	//Calculate mean
	double mean = IOUtils::nodata;
	unsigned int vecSize = vecWindow.size();

	if (vecSize == 0){
		return false; //only nodata values detected or other problem
	} else {
		double sum = 0.0;
		for (unsigned int ii=0; ii<vecSize; ii++){
			sum += vecWindow[ii];
		}
		mean = sum/vecSize;
	}

	//cout << "Mean value: " << mean << endl;
	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		vecFilteredM[ii].param(paramindex) = mean;
	}

	return true;
}

bool FilterAlgorithms::WindAvgFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				   const unsigned int& pos, const Date_IO& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex,
				   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecS; (void)vecFilteredS; (void)paramindex;
	//Remarks:
	//1) nodata values are not considered when calculating the mean
	//2) Two arguments expected (both have to be fullfilled for the filter to start operating): 
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//3) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//4) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
	//  
	std::vector<double> vecWindowVW, vecWindowDW;
	std::vector<Date_IO> vecDateVW, vecDateDW;
	if (!getWindowData("wind_avg", vecM, pos, date, _vecArgs, MeteoData::VW, vecWindowVW, &vecDateVW))
		return false; //Not enough data to meet user configuration

	if (!getWindowData("wind_avg", vecM, pos, date, _vecArgs, MeteoData::DW, vecWindowDW, &vecDateDW))
		return false; //Not enough data to meet user configuration

	if (vecWindowVW.size() != vecWindowDW.size()) //same amount of data points necessary
		return false;

	for (unsigned int ii=0; ii<vecDateVW.size(); ii++){ //The VW and DW data points have to correlate
		if (vecDateVW[ii] != vecDateDW[ii]) return false;
	}

	//Calculate mean
	double meanspeed     = IOUtils::nodata;
	double meandirection = IOUtils::nodata;
	unsigned int vecSize = vecWindowVW.size();

	if (vecSize == 0){
		return false; //only nodata values detected or other problem
	} else {
		//calculate ve and vn
		double ve=0.0, vn=0.0;
		for (unsigned int ii=0; ii<vecSize; ii++){
			ve += vecWindowVW[ii] * sin(vecWindowDW[ii] * PI / 180); //turn into radians
			vn += vecWindowVW[ii] * cos(vecWindowDW[ii] * PI / 180); //turn into radians
		}
		ve = (-1.0)*ve/vecSize;
		vn = (-1.0)*vn/vecSize;

		meanspeed = sqrt(ve*ve + vn*vn);
		meandirection = atan2(ve,vn) * 180 / PI + 180; // turn into degrees [0;360)
	}

	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		vecFilteredM[ii].param(MeteoData::VW) = meanspeed;
		vecFilteredM[ii].param(MeteoData::DW) = meandirection;
	}

	return true;
}

bool FilterAlgorithms::getWindowData(const std::string& filtername, const std::vector<MeteoData>& vecM, 
				   const unsigned int& pos, 
				   const Date_IO& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex, std::vector<double>& vecWindow, 
				   std::vector<Date_IO> *vecDate)
{
	vecWindow.clear();
	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> vecArgs; 
	parseWindowFilterArguments(filtername, _vecArgs, 2, 2, isSoft, windowposition, vecArgs);

	if (vecArgs[0] < 1) //the window size has to be at least 1
		throw InvalidArgumentException("Number of data points in window of " +filtername+ " filter cannot be < 1", AT);
	if (vecArgs[1] < 0) //the time window has to be at least 0 minutes
		throw InvalidArgumentException("Time span of window for filter " +filtername+ " cannot be < 0 minutes", AT);

	//Deal with first parameter: minimal number of data points
	unsigned int windowSize = (unsigned int)vecArgs[0];
	unsigned int increment = 0;
	if (windowposition=="right"){ 
		increment = (unsigned int)windowSize - 1; //window reaches to the right
		if ((date != vecM[pos].date) && (increment>0)) increment--;
	} else if (windowposition=="left") {
		increment = 0;                        //window reaches to the left
	} else {
		increment = (unsigned int)(windowSize/2);                   //window is centered
		if (increment>0) increment--; //reason: the element at pos might be of greater date
	}

	unsigned int startposition = pos + increment;

	if (startposition > (vecM.size()-1)){
		if (!isSoft) return false; //if "soft" is not defined then there have to be enough data elements to the right
		startposition = (vecM.size()-1);
	}

	unsigned int endposition = 0;
	if (startposition < (windowSize-1)){
		return false; //not enough data points available
	} else {
		endposition = startposition - (windowSize - 1);
	}

	//Now deal with the second argument: the time window
	Date_IO deltatime(vecArgs[1]/1440.0); //making a julian date out of the argument given in minutes	
	while(deltatime > (vecM[startposition].date - vecM[endposition].date)){
		//The time window is too small, so let's try enlargening it
		if ((windowposition == "right") && (!isSoft)){ //only expand to the right
			if (startposition<(vecM.size()-1)) startposition++;
			else return false; //time window not big enough
		} else if ((windowposition == "right") && (isSoft)){
			if (startposition<(vecM.size()-1)) startposition++;
			else if (endposition > 0) endposition--;
			else return false; //time window not big enough
		} else if ((windowposition == "left") && (!isSoft)){ //only expand to the left
			if (endposition > 0) endposition--;
			else return false; //time window not big enough
		} else if ((windowposition == "left") && (isSoft)){
			if (endposition > 0) endposition--;
			else if (startposition<(vecM.size()-1)) startposition++;
			else return false; //time window not big enough
		} else if ((windowposition == "center") && (!isSoft)){ //expand to left and right (equally)
			if ((endposition+startposition)%2 == 0) { //try to alternate when broadening window
				if (endposition > 0) endposition--;
				else return false; //time window not big enough
			} else {
				if (startposition<(vecM.size()-1)) startposition++;
				else return false; //time window not big enough
			}			
		} else if ((windowposition == "center") && (isSoft)){ //expand to left and right (whereever possible)
			if ((endposition+startposition)%2 == 0) { //try to alternate when broadening window
				if (endposition > 0) endposition--;
				else if (startposition<(vecM.size()-1)) startposition++;
				else return false; //time window not big enough
			} else {
				if (startposition<(vecM.size()-1)) startposition++;
				else if (endposition > 0) endposition--;
				else return false; //time window not big enough
			}			
		}
	} 
	
	//Date_IO gap(vecM[startposition].date - vecM[endposition].date);
	//cout << "The final gap is " << gap.toString() << "  windowposition: " << windowposition << endl;

	//Push all relevant data elements into the vector vecWindow
	//cout << "Start: " << startposition << "  End: " << endposition << endl;
	for (unsigned int ii=startposition; ii>=endposition; ii--){
		const double& tmp = vecM[ii].param(paramindex);
		if (tmp != IOUtils::nodata) {
			vecWindow.push_back(tmp);
			if (vecDate != NULL) (*vecDate).push_back(vecM[ii].date);
		}
		
		//cout << ii << ": pushed at vecM[" <<  ii << "] " << vecM[ii].date << " : " << tmp << endl;
	}

	return true;
}

