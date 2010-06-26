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
#include <meteoio/FilterAlgorithms.h>

using namespace std;

namespace mio {
 /**
 * @page filters Filters overview
 * The filtering infrastructure is described in FilterAlgorithms (for its API). The goal of this page is to give an overview of the available filters and their usage.
 *
 * @section filters_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * @section filters_available Available filters
 * The filters that are currently available are the following:
 * - rate: rate of change filter, see FilterAlgorithms::RateFilter
 * - min_max: range check filter, see FilterAlgorithms::MinMaxFilter
 * - min: minimum check filter, see FilterAlgorithms::MinValueFilter
 * - max: maximum check filter, see FilterAlgorithms::MaxValueFilter
 * - mad: median absolute deviation, see FilterAlgorithms::MedianAbsoluteDeviationFilter
 *
 * A few data transformations are also supported besides filtering:
 * - accumulate: data accumulates over a given period, see FilterAlgorithms::AccumulateProcess
 * - median_avg: running median average over a given window, see FilterAlgorithms::MedianAvgProcess
 * - mean_avg: running mean average over a given window, see FilterAlgorithms::MeanAvgProcess
 * - wind_avg: vector average over a given window, see FilterAlgorithms::WindAvgProcess
 */

std::map<std::string, FilterProperties> FilterAlgorithms::filterMap;
const bool FilterAlgorithms::__init = FilterAlgorithms::initStaticData();

bool FilterAlgorithms::initStaticData()
{
	filterMap["rate"]          = FilterProperties(true,  &FilterAlgorithms::RateFilter);
	filterMap["min_max"]       = FilterProperties(true,  &FilterAlgorithms::MinMaxFilter);
	filterMap["min"]           = FilterProperties(true,  &FilterAlgorithms::MinValueFilter);
	filterMap["max"]           = FilterProperties(true,  &FilterAlgorithms::MaxValueFilter);
	filterMap["mad"]           = FilterProperties(true,  &FilterAlgorithms::MedianAbsoluteDeviationFilter);
	filterMap["accumulate"]    = FilterProperties(false, &FilterAlgorithms::AccumulateProcess);
	filterMap["exp_smoothing"] = FilterProperties(false, &FilterAlgorithms::ExpSmoothingFilter);
	filterMap["median_avg"]    = FilterProperties(false, &FilterAlgorithms::MedianAvgProcess);
	filterMap["mean_avg"]      = FilterProperties(false, &FilterAlgorithms::MeanAvgProcess);
	filterMap["wind_avg"]      = FilterProperties(false, &FilterAlgorithms::WindAvgProcess);

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

bool FilterAlgorithms::getWindowData(const std::string& filtername, const std::vector<MeteoData>& vecM,
                                     const unsigned int& pos,
                                     const Date& date, const std::vector<std::string>& _vecArgs,
                                     const unsigned int& paramindex, std::vector<double>& vecWindow,
                                     std::vector<Date> *vecDate)
{
	vecWindow.clear();
	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> vecArgs;
	parseWindowFilterArguments(filtername, _vecArgs, 2, 2, isSoft, windowposition, vecArgs);

	if (vecArgs[0] < 1) //the window size has to be at least 1
		throw InvalidArgumentException("Number of data points in window of " +filtername+ " filter cannot be < 1", AT);
	if (vecArgs[1] < 0) //the time window has to be at least 0 seconds
		throw InvalidArgumentException("Time span of window for filter " +filtername+ " cannot be < 0 seconds", AT);

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
	Date deltatime(vecArgs[1]/(24.*3600.)); //making a julian date out of the argument given in seconds
	while((deltatime > (vecM[startposition].date - vecM[endposition].date)) && (deltatime != 0.0)){
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

	//Date gap(vecM[startposition].date - vecM[endposition].date);
	//cout << "The final gap is " << gap << "  windowposition: " << windowposition << endl;

	//Push all relevant data elements into the vector vecWindow
	//cout << "Start: " << startposition << "  End: " << endposition << endl;
	for (unsigned int ii=startposition+1; (ii--) > endposition; ){
		const double& tmp = vecM[ii].param(paramindex);
		if (tmp != IOUtils::nodata) {
			vecWindow.push_back(tmp);
			if (vecDate != NULL) (*vecDate).push_back(vecM[ii].date);
		}
		//cout << ii << ": pushed at vecM[" <<  ii << "] " << " : "  << endl;
	}

	return true;
}



/******************************************************************************
 * The following functions are implementations of different filter algorithms *
 ******************************************************************************/

/**
 * @brief Exponential smooting filter, exponential moving average
 * s_0 = x_0
 * s_n = alpha*x_(t-1) + (1-alpha)*s_t-1
 * - http://en.wikipedia.org/wiki/Exponential_smoothing
 * - alpha needs to be provided as argument, as well as the window size and centering 
 * - nodata values are excluded from the moving average calculation
 * - Three arguments expected (all have to be present and valid for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 *   - alpha value 0 < alpha < 1
 * - the arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added at the beginning, if the window position is allowed to be adjusted to the data present
 * @code
 *          TA::filter1 = exp_smoothing
 *          TA::arg1    = right 1 1800 0.8 (1800 seconds time span for the strictly right leaning window), alpha=0.8
 *          RH::filter1 = mean_avg
 *          RH::arg1    = 10 600 0.6 (strictly left window spanning 600 seconds and at least 10 points), alpha=0.6
 * @endcode
 */
void FilterAlgorithms::ExpSmoothingFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				   const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                       std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS;
	if (vecArgs.size() < 3)
		throw InvalidArgumentException("Wrong number of arguments for ExpSmoothingFilter", AT);

	vector<string> myArgs = vecArgs;

	//Take away the last argument (alpha value) and check validity
	double alpha = 0.5;
	IOUtils::convertString(alpha, myArgs.at(myArgs.size()-1)); 
	if ((alpha <= 0) || (alpha >= 1))
		throw InvalidArgumentException("ExpSmoothingFilter: alpha must be in [0;1]", AT);

	myArgs.pop_back(); //delete alpha from myArgs


	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> doubleArgs;
	parseWindowFilterArguments("exp_smoothing", myArgs, 2, 2, isSoft, windowposition, doubleArgs);

	//for every element in vecWindowM, get the Window and perform exponential smoothing
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		vector<MeteoData> vecTmpWindow;
		unsigned int position = IOUtils::seek(vecWindowM[ii].date, vecM);

		if (position != IOUtils::npos){
			unsigned int posfind = getWindowData("exp_smoothing", vecM, position, myArgs, vecTmpWindow);
			if (posfind == IOUtils::npos)
				continue;

			if (windowposition == "left"){
				vecTmpWindow.erase(vecTmpWindow.begin()+posfind+1, vecTmpWindow.end()); //delete all after posfind
				vecWindowM[ii].param(paramindex) = ExpSmoothingAlgorithm(vecTmpWindow, paramindex, alpha);
				//cout << "ExpSmoothing: " << vecWindowM[ii].param(paramindex) << endl;
			} else if (windowposition == "right"){
				vecTmpWindow.erase(vecTmpWindow.begin(), vecTmpWindow.begin()+posfind); //delete all before posfind
				std::reverse(vecTmpWindow.begin(), vecTmpWindow.end()); //reverse the vector, posfind most significant
				vecWindowM[ii].param(paramindex) = ExpSmoothingAlgorithm(vecTmpWindow, paramindex, alpha);
				//cout << "ExpSmoothing: " << vecWindowM[ii].param(paramindex) << endl;
			} else { //centered window - regroup according to time difference with posfind
				for (unsigned int jj=0; jj<vecTmpWindow.size(); jj++)
					vecTmpWindow[jj].date=Date(abs(vecWindowM[ii].date.getJulianDate() - vecTmpWindow[jj].date.getJulianDate()));
				std::sort(vecTmpWindow.begin(), vecTmpWindow.end(), compareMeteoData);
				vecWindowM[ii].param(paramindex) = ExpSmoothingAlgorithm(vecTmpWindow, paramindex, alpha);
				//cout << "ExpSmoothing: " << vecWindowM[ii].param(paramindex) << endl;
			}
		}	
	}
}

bool FilterAlgorithms::compareMeteoData (const MeteoData& m1, const MeteoData& m2)
{ 
	return (m1.date>m2.date);
}

unsigned int FilterAlgorithms::getWindowData(const std::string& filtername, const std::vector<MeteoData>& vecM, 
						    const unsigned int& pos, 
						    const std::vector<std::string>& _vecArgs, std::vector<MeteoData>& vecResult)
{
	Date date(vecM.at(pos).date);
	vecResult.clear();
	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> vecArgs;
	parseWindowFilterArguments(filtername, _vecArgs, 2, 2, isSoft, windowposition, vecArgs);

	if (vecArgs[0] < 1) //the window size has to be at least 1
		throw InvalidArgumentException("Number of data points in window of " +filtername+ " filter cannot be < 1", AT);
	if (vecArgs[1] < 0) //the time window has to be at least 0 seconds
		throw InvalidArgumentException("Time span of window for filter " +filtername+ " cannot be < 0 seconds", AT);

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
		if (!isSoft) return IOUtils::npos; //if "!isSoft" then there have to be enough data elements to the right
		startposition = (vecM.size()-1);
	}

	unsigned int endposition = 0;
	if (startposition < (windowSize-1)){
		return IOUtils::npos; //not enough data points available
	} else {
		endposition = startposition - (windowSize - 1);
	}

	//Now deal with the second argument: the time window
	Date deltatime(vecArgs[1]/(24.*3600.)); //making a julian date out of the argument given in seconds
	while((deltatime > (vecM[startposition].date - vecM[endposition].date)) && (deltatime != 0.0)){
		//The time window is too small, so let's try enlargening it
		if ((windowposition == "right") && (!isSoft)){ //only expand to the right
			if (startposition<(vecM.size()-1)) startposition++;
			else return IOUtils::npos; //time window not big enough
		} else if ((windowposition == "right") && (isSoft)){
			if (startposition<(vecM.size()-1)) startposition++;
			else if (endposition > 0) endposition--;
			else return IOUtils::npos; //time window not big enough
		} else if ((windowposition == "left") && (!isSoft)){ //only expand to the left
			if (endposition > 0) endposition--;
			else return IOUtils::npos; //time window not big enough
		} else if ((windowposition == "left") && (isSoft)){
			if (endposition > 0) endposition--;
			else if (startposition<(vecM.size()-1)) startposition++;
			else return IOUtils::npos; //time window not big enough
		} else if ((windowposition == "center") && (!isSoft)){ //expand to left and right (equally)
			if ((endposition+startposition)%2 == 0) { //try to alternate when broadening window
				if (endposition > 0) endposition--;
				else return IOUtils::npos; //time window not big enough
			} else {
				if (startposition<(vecM.size()-1)) startposition++;
				else return IOUtils::npos; //time window not big enough
			}
		} else if ((windowposition == "center") && (isSoft)){ //expand to left and right (whereever possible)
			if ((endposition+startposition)%2 == 0) { //try to alternate when broadening window
				if (endposition > 0) endposition--;
				else if (startposition<(vecM.size()-1)) startposition++;
				else return IOUtils::npos; //time window not big enough
			} else {
				if (startposition<(vecM.size()-1)) startposition++;
				else if (endposition > 0) endposition--;
				else return IOUtils::npos; //time window not big enough
			}
		}
	}

	//Date gap(vecM[startposition].date - vecM[endposition].date);
	//cout << "The final gap is " << gap << "  windowposition: " << windowposition << endl;

	//Push all relevant data elements into the vector vecResult
	//cout << "POS" << pos << "  start: " << startposition << "  end: " << endposition << endl;
	unsigned int posofposition = IOUtils::npos;
	unsigned int counter = 0;
	for (unsigned int ii=endposition; ii<startposition+1; ii++){
		vecResult.push_back(vecM[ii]);
		if (date == vecM[ii].date) posofposition=counter;
		counter++;
		//cout << ii << ": pushed at vecM[" <<  ii << "] " << " : "  << endl;
	}

	return posofposition;
}

/**
 * @brief Actual implementation of the smoothing algorithm
 * @param vecMeteo A vector of MeteoData
 * @param paramindex The meteo data parameter to be smoothed
 * @param alpha The alpha for the exponential smoothing, the smoothing factor
 * @return The smoothing value for the specified meteo parameter of the last element in vecMeteo
 */
double FilterAlgorithms::ExpSmoothingAlgorithm(const std::vector<MeteoData>& vecMeteo, const unsigned int& paramindex,
                                               const double& alpha)
{
	if (vecMeteo.size() == 0)
		return IOUtils::nodata;

	bool initCompleted = false;
	double expavg = IOUtils::nodata;
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		const double& currentval = vecMeteo[ii].param(paramindex);
		if (currentval != IOUtils::nodata){
			if (!initCompleted){
				expavg = currentval;
				initCompleted = true;
			} else {
				expavg = alpha*currentval + (1-alpha)*expavg;
			}
		}
	}

	return expavg;
}

/**
 * @brief Rate of change filter.
 * Calculate the change rate (ie: slope) between two points, if it is above a user given value, reject the point. Remarks:
 * - the maximum permissible rate of change (per seconds) has to be provided as an argument
 * @code
 * TA::filter1	= rate
 * TA::arg1	= 0.01
 * @endcode
 */
void FilterAlgorithms::RateFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
                                  const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                  std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS;

	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("rate", vecArgs, 1, 1, isSoft, doubleArgs);
	/*
	const double& maxRateOfChange = doubleArgs[0];

	//Run actual Rate filter over all relevant meteo data
	for (unsigned int ii=2; ii<vecWindowM.size(); ii++){
		MeteoData& current     = vecWindowM[ii];
		MeteoData& previous    = vecWindowM[ii-1];
		MeteoData& preprevious = vecWindowM[ii-2];
	}


	
	if(vecWindowM.size()>1) {
		//the request time step is NOT part of the data, it will have to be resampled
		const double prev_value = vecM[pos-vecWindowM.size()].param(paramindex);
		const double prev_time = vecM[pos-vecWindowM.size()].date.getJulianDate()*24.*3600.; //in seconds

		for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
			double& value = vecWindowM[ii].param(paramindex);
			const double curr_value = vecWindowM[ii].param(paramindex);
			const double curr_time = vecWindowM[ii].date.getJulianDate()*24.*3600.; //in seconds
			const double local_rate = abs((curr_value-prev_value)/(curr_time-prev_time));

			if( local_rate > doubleArgs[0] ) {
				value = IOUtils::nodata;
			}
		}
	} else {
		//const double tmp = vecWindowM.size();
		//std::cout << tmp << std::endl;
		//the request time step is part of the data
		if(pos>1) { //if we are at the start of the data set, we can not apply the filter...
			double& value = vecWindowM[0].param(paramindex);
			const double curr_value = vecM[pos].param(paramindex);
			const double curr_time = vecM[pos].date.getJulianDate()*24.*3600.; //in seconds
			const double prev_value = vecM[pos-1].param(paramindex);
			const double prev_time = vecM[pos-1].date.getJulianDate()*24.*3600.; //in seconds
			const double local_rate = abs((curr_value-prev_value)/(curr_time-prev_time));

			if( local_rate > doubleArgs[0] ) {
				value = IOUtils::nodata;
			}
		} else {
			//the filter can not be applied
		}
	}
	*/
}

/**
 * @brief Min/Max range filter.
 * Reject all values greater than the max or smaller than the min. Remarks:
 * - two arguments have to be provided, min and max (in SI)
 * - the keyword "soft" maybe added, in such a case all data greater than the max would be assigned
 * the maximum permissible value and all data smaller than the min would be assigned the minimum permissible value
 * @code
 * TA::filter1	= min_max
 * TA::arg1	= 230 330
 * @endcode
 */
void FilterAlgorithms::MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				                const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                    std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecM; (void)vecS; (void)vecWindowS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("min_max", vecArgs, 2, 2, isSoft, doubleArgs);

	sort(doubleArgs.begin(), doubleArgs.end()); //the two parameters are sorted ascending

	//Run actual MinMax filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii].param(paramindex);

		if (value == IOUtils::nodata) continue;

		if (value<doubleArgs[0]){
			if (isSoft) value=doubleArgs[0];
			else value=IOUtils::nodata;
			//cout << "Changed: " << value << endl;
		}
		if (value>doubleArgs[1]){
			if (isSoft) value=doubleArgs[1];
			else value=IOUtils::nodata;
			//cout << "Changed: " << value << endl;
		}
	}
}

/**
 * @brief Min range filter.
 * Reject all values smaller than the min. Remarks:
 * - the minimum permissible value has to be provided has an argument (in SI)
 * - the keyword "soft" maybe added, in such a case all data smaller than the min would be assigned
 * the minimum permissible value
 * @code
 * TA::filter1	= min
 * TA::arg1	= 230
 * @endcode
 */
void FilterAlgorithms::MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				                  const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                      std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecM; (void)vecS; (void)vecWindowS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("min", vecArgs, 1, 1, isSoft, doubleArgs);

	//Run actual MinValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii].param(paramindex);

		if (value == IOUtils::nodata) continue;

		if (value<doubleArgs[0]){
			if (isSoft) value=doubleArgs[0];
			else value=IOUtils::nodata;
		}
	}
}

/**
 * @brief Max range filter.
 * Reject all values greater than the max. Remarks:
 * - the maximum permissible value has to be provided has an argument (in SI)
 * - the keyword "soft" maybe added, in such a case all data greater than the max would be assigned
 * the maximum permissible value
 * @code
 * TA::filter1	= max
 * TA::arg1	= 330
 * @endcode
 */
void FilterAlgorithms::MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				                  const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                      std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecM; (void)vecS; (void)vecWindowS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("max", vecArgs, 1, 1, isSoft, doubleArgs);

	//Run actual MaxValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii].param(paramindex);

		if (value == IOUtils::nodata) continue;

		if (value>doubleArgs[0]){
			if (isSoft) value=doubleArgs[0];
			else value=IOUtils::nodata;
		}
	}
}

/**
 * @brief Median Absolute Deviation.
 * Values outside of median ± 3 σ_MAD are rejected. The σ_MAD is calculated as follow:\n
 * <center>\f$ \sigma_{MAD} = K \cdot \mathop{median_i} \left( \left| X_i - \mathop{median_j} ( X_j ) \right| \right) \f$ with \f$ K = \Phi^{-1}; \Phi = 0.6745 \f$ </center>\n
 * See http://en.wikipedia.org/wiki/Median_absolute_deviation
 * for more information about the Mean Absolute Deviation.
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = mad
 *          TA::arg1    = soft left 1 1800  (1800 seconds time span for the left leaning window)
 *          RH::filter1 = mad
 *          RH::arg1    = 10 600            (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */
void FilterAlgorithms::MedianAbsoluteDeviationFilter(const std::vector<MeteoData>& vecM, 
                                      const std::vector<StationData>& vecS, 
				                  const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                      std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS;

	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		unsigned int pos = IOUtils::seek(vecWindowM[ii].date, vecM, false);

		std::vector<double> vecWindow;
		if (!getWindowData("mad", vecM, pos, vecWindowM[ii].date, vecArgs, paramindex, vecWindow))
			return; //Not enough data to meet user configuration
	
		//Calculate MAD
		const double K = 1. / 0.6745;
		double mad     = IOUtils::nodata;
		double median  = IOUtils::nodata;

		try {
			median = Interpol1D::getMedian(vecWindow);
			mad    = Interpol1D::getMedianAverageDeviation(vecWindow);
		} catch(exception& e){
			return;
		}
	
		double sigma = mad * K;

		double& value = vecWindowM[ii].param(paramindex);
		if( (value>(median + 3.*sigma)) || (value<(median - 3.*sigma)) ) {
			value = IOUtils::nodata;
		}
	}
}


/**
 * @brief Accumulation over a user given period.
 * The input data is accumulated over a given time interval (given as filter argument, in minutes).
 * This is for example needed for converting rain gauges measurements read every 10 minutes to
 * hourly precipitation measurements. Remarks:
 * - the accumulation period has to be provided as an argument (in seconds)
 * @code
 * HNW::filter1 = accumulate
 * HNW::arg1	 = 3600
 * @endcode
 */
void FilterAlgorithms::AccumulateProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				             const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                 std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS;
	
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("accumulate", vecArgs, 1, 1, isSoft, doubleArgs);

	Date deltatime(doubleArgs[0]/(24.*3600.)); //making a julian date out of the argument given in seconds

	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		unsigned int pos = IOUtils::seek(vecWindowM[ii].date, vecM);

		if (pos != IOUtils::npos){
			unsigned int startpos = pos;
			double sum = 0.0;
			while((vecM[pos].date + deltatime) > vecM[startpos].date){
				const double& val = vecM[pos].param(paramindex);

				if (val != IOUtils::nodata)
					sum += val;

				if (pos > 0) pos--;
				else break;
			}
			vecWindowM[ii].param(paramindex) = sum;
			//cout << "sum: " << vecWindowM[ii].param(paramindex) << endl;
		} else {
			throw IOException("Could not find an element to start accumulation", AT);
		}
	}
}

/**
 * @brief Median averaging.
 * The median average filter returns the median value of all values within a user given time window. Remarks:
 * - nodata values are excluded from the median
 * - if there is an even number of window elements the arithmetic mean of the two central elements is used to calculate the median
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
 *
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = median_avg
 *          TA::arg1    = soft left 1 1800  (1800 seconds time span for the left leaning window)
 *          RH::filter1 = median_avg
 *          RH::arg1    = 10 600            (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */
void FilterAlgorithms::MedianAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				             const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                 std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS;
	//Remarks:
	//1) nodata values are not considered when calculating the median
	//2) if there is an even number of window elements the arithmetic mean is used to calculate the median
	//3) Two arguments expected (both have to be fullfilled for the filter to start operating):
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//4) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//5) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present

	//for every element in vecFilteredM, get the Window and calculate median average
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		std::vector<double> vecWindow;
		unsigned int position = IOUtils::seek(vecWindowM[ii].date, vecM);

		if (!getWindowData("median_avg", vecM, position, vecWindowM[ii].date, vecArgs, paramindex, vecWindow))
			continue; //Not enough data to meet user configuration

		double median = IOUtils::nodata;
		try {
			median = Interpol1D::getMedian(vecWindow);
			vecWindowM[ii].param(paramindex) = median;
		} catch(exception& e){
			continue; //the median calculation did not work out, filter is not applied, value unchanged
		}
	}
}

/**
 * @brief Mean averaging.
 * The mean average filter returns the mean value of all values within a user given time window. Remarks:
 * - nodata values are excluded from the mean
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
 *
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = mean_avg
 *          TA::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 *          RH::filter1 = mean_avg
 *          RH::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */
void FilterAlgorithms::MeanAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				             const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                 std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS;
	//Remarks:
	//1) nodata values are not considered when calculating the mean
	//2) Two arguments expected (both have to be fullfilled for the filter to start operating):
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//3) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//4) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present

	//for every element in vecWindowM, get the Window and calculate mean average
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		std::vector<double> vecWindow;
		unsigned int position = IOUtils::seek(vecWindowM[ii].date, vecM);

		if (!getWindowData("mean_avg", vecM, position, vecWindowM[ii].date, vecArgs, paramindex, vecWindow))
			continue; //Not enough data to meet user configuration

		double mean = IOUtils::nodata;
		try {
			mean = Interpol1D::arithmeticMean(vecWindow);
			vecWindowM[ii].param(paramindex) = mean;
		} catch(exception& e){
			continue; //the mean calculation did not work out, filter is not applied, value unchanged
		}
	}
}

/**
 * @brief Wind vector averaging.
 * This calculates the vector average over a user given time period. Each wind vector within this period
 * is added and the final sum is normalized by the number of vectors that have been added. A few more important information:
 * - nodata values are excluded from the mean
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
 * @code
 * Valid examples for the io.ini file:
 *          VW::filter1 = wind_avg
 *          VW::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 *          VW::filter1 = wind_avg
 *          VW::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */
void FilterAlgorithms::WindAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
				             const std::vector<std::string>& vecArgs, const MeteoData::Parameters& paramindex,
                                 std::vector<MeteoData>& vecWindowM, std::vector<StationData>& vecWindowS)
{
	(void)vecS; (void)vecWindowS; (void)paramindex;
	//Remarks:
	//1) nodata values are not considered when calculating the mean
	//2) Two arguments expected (both have to be fullfilled for the filter to start operating):
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//3) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//4) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
	//
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		unsigned int pos = IOUtils::seek(vecWindowM[ii].date, vecM);
		
		if (pos == IOUtils::npos)
			continue;

		std::vector<double> vecWindowVW, vecWindowDW;
		std::vector<Date> vecDateVW, vecDateDW;
		if (!getWindowData("wind_avg", vecM, pos, vecWindowM[ii].date, vecArgs, MeteoData::VW, vecWindowVW, &vecDateVW))
			continue; //Not enough data to meet user configuration

		if (!getWindowData("wind_avg", vecM, pos, vecWindowM[ii].date, vecArgs, MeteoData::DW, vecWindowDW, &vecDateDW))
			continue; //Not enough data to meet user configuration

		if (vecWindowVW.size() != vecWindowDW.size()) //same amount of data points necessary
			continue;

		for (unsigned int jj=0; jj<vecDateVW.size(); jj++){ //The VW and DW data points have to correlate
			if (vecDateVW[jj] != vecDateDW[jj]) continue;
		}

		//Calculate mean
		double meanspeed     = IOUtils::nodata;
		double meandirection = IOUtils::nodata;
		unsigned int vecSize = vecWindowVW.size();

		if (vecSize == 0){
			continue; //only nodata values detected or other problem
		} else {
			//calculate ve and vn
			double ve=0.0, vn=0.0;
			for (unsigned int jj=0; jj<vecSize; jj++){
				ve += vecWindowVW[jj] * sin(vecWindowDW[jj] * M_PI / 180.); //turn into radians
				vn += vecWindowVW[jj] * cos(vecWindowDW[jj] * M_PI / 180.); //turn into radians
			}
			ve /= vecSize;
			vn /= vecSize;

			meanspeed = sqrt(ve*ve + vn*vn);
			meandirection = fmod( atan2(ve,vn) * 180. / M_PI + 360. , 360.); // turn into degrees [0;360)
		}

		vecWindowM[ii].param(MeteoData::VW) = meanspeed;
		vecWindowM[ii].param(MeteoData::DW) = meandirection;
	}
}

} //namespace

