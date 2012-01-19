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
#include <cmath>
#include <meteoio/FilterAlgorithms.h>
#include <meteoio/meteolaws/Meteoconst.h> //for math constants

using namespace std;

namespace mio {

std::map<std::string, FilterProperties> FilterAlgorithms::filterMap;
const bool FilterAlgorithms::__init = FilterAlgorithms::initStaticData();

bool FilterAlgorithms::initStaticData()
{
	filterMap["rate"]          = FilterProperties(true,  &FilterAlgorithms::RateFilter);
	filterMap["min_max"]       = FilterProperties(true,  &FilterAlgorithms::MinMaxFilter);
	filterMap["min"]           = FilterProperties(true,  &FilterAlgorithms::MinValueFilter);
	filterMap["max"]           = FilterProperties(true,  &FilterAlgorithms::MaxValueFilter);
	filterMap["mad"]           = FilterProperties(true,  &FilterAlgorithms::MedianAbsoluteDeviationFilter);
	filterMap["stddev"]        = FilterProperties(true,  &FilterAlgorithms::StandardDeviationFilter);
	filterMap["Tukey53H"]      = FilterProperties(true,  &FilterAlgorithms::Tukey53HFilter);
	filterMap["accumulate"]    = FilterProperties(false, &FilterAlgorithms::AccumulateProcess);
	filterMap["exp_smoothing"] = FilterProperties(false, &FilterAlgorithms::ExpSmoothingProcess);
	filterMap["wma_smoothing"] = FilterProperties(false, &FilterAlgorithms::WMASmoothingProcess);
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
	} catch(const std::exception&){
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
                                     const Date& date, const std::vector<std::string>& i_vecArgs,
                                     const unsigned int& paramindex, std::vector<double>& vecWindow,
                                     std::vector<Date> *vecDate)
{
	vecWindow.clear();
	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> vecArgs;
	parseWindowFilterArguments(filtername, i_vecArgs, 2, 2, isSoft, windowposition, vecArgs);

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
		const double& tmp = vecM[ii](paramindex);
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
 * @brief Exponential smooting processing, exponential moving average
 * s_0 = x_0
 * s_n = alpha*x_(t-1) + (1-alpha)*s_t-1
 * - http://en.wikipedia.org/wiki/Exponential_smoothing
 * - alpha needs to be provided as argument, as well as the window size and centering
 * - nodata values are excluded from the moving average calculation
 * - Three arguments expected (all have to be present and valid for the algorithm to start operating):
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
void FilterAlgorithms::ExpSmoothingProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                          const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
//HACK: this should be a process, not a filter!!
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
				vecWindowM[ii](paramindex) = ExpSmoothingAlgorithm(vecTmpWindow, paramindex, alpha);
				//cout << "ExpSmoothing: " << vecWindowM[ii](paramindex) << endl;
			} else if (windowposition == "right"){
				vecTmpWindow.erase(vecTmpWindow.begin(), vecTmpWindow.begin()+posfind); //delete all before posfind
				std::reverse(vecTmpWindow.begin(), vecTmpWindow.end()); //reverse the vector, posfind most significant
				vecWindowM[ii](paramindex) = ExpSmoothingAlgorithm(vecTmpWindow, paramindex, alpha);
				//cout << "ExpSmoothing: " << vecWindowM[ii](paramindex) << endl;
			} else { //centered window - regroup according to time difference with posfind
				for (unsigned int jj=0; jj<vecTmpWindow.size(); jj++)
					vecTmpWindow[jj].date=Date(abs(vecWindowM[ii].date.getJulianDate() - vecTmpWindow[jj].date.getJulianDate()));
				std::sort(vecTmpWindow.begin(), vecTmpWindow.end(), compareMeteoData);
				vecWindowM[ii](paramindex) = ExpSmoothingAlgorithm(vecTmpWindow, paramindex, alpha);
				//cout << "ExpSmoothing: " << vecWindowM[ii](paramindex) << endl;
			}
		}
	}
}

/**
 * @brief WMA smooting processing, weighted moving average
 * - http://en.wikipedia.org/wiki/Simple_moving_average#Weighted_moving_average
 * - nodata values are excluded from the moving average calculation
 * - Two arguments expected (all have to be present and valid for the algorithm to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added at the beginning, if the window position is allowed to be adjusted to the data present
 * @code
 *          TA::filter1 = wma_smoothing
 *          TA::arg1    = right 1 1800 (1800 seconds time span for the strictly right leaning window)
 * @endcode
 */
void FilterAlgorithms::WMASmoothingProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                          const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	if (vecArgs.size() < 2)
		throw InvalidArgumentException("Wrong number of arguments for WMASmoothingFilter", AT);

	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> doubleArgs;
	parseWindowFilterArguments("wma_smoothing", vecArgs, 2, 2, isSoft, windowposition, doubleArgs);

	//for every element in vecWindowM, get the Window and perform weighted moving average smoothing
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		vector<MeteoData> vecTmpWindow;
		unsigned int position = IOUtils::seek(vecWindowM[ii].date, vecM);

		if (position != IOUtils::npos){
			unsigned int posfind = getWindowData("wma_smoothing", vecM, position, vecArgs, vecTmpWindow);
			if (posfind == IOUtils::npos)
				continue;

			if (windowposition == "left"){
				vecTmpWindow.erase(vecTmpWindow.begin()+posfind+1, vecTmpWindow.end()); //delete all after posfind
				vecWindowM[ii](paramindex) = WMASmoothingAlgorithm(vecTmpWindow, paramindex);
				//cout << "WMASmoothing: " << vecWindowM[ii](paramindex) << endl;
			} else if (windowposition == "right"){
				vecTmpWindow.erase(vecTmpWindow.begin(), vecTmpWindow.begin()+posfind); //delete all before posfind
				std::reverse(vecTmpWindow.begin(), vecTmpWindow.end()); //reverse the vector, posfind most significant
				vecWindowM[ii](paramindex) = WMASmoothingAlgorithm(vecTmpWindow, paramindex);
				//cout << "WMASmoothing: " << vecWindowM[ii](paramindex) << endl;
			} else { //centered window - regroup according to time difference with posfind
				for (unsigned int jj=0; jj<vecTmpWindow.size(); jj++)
					vecTmpWindow[jj].date=Date(abs(vecWindowM[ii].date.getJulianDate() - vecTmpWindow[jj].date.getJulianDate()));
				std::sort(vecTmpWindow.begin(), vecTmpWindow.end(), compareMeteoData);
				vecWindowM[ii](paramindex) = WMASmoothingAlgorithm(vecTmpWindow, paramindex);
				//cout << "WMASmoothing: " << vecWindowM[ii](paramindex) << endl;
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
                                             const std::vector<std::string>& i_vecArgs, std::vector<MeteoData>& vecResult)
{
	Date date(vecM.at(pos).date);
	vecResult.clear();
	bool isSoft = false;
	std::string windowposition = "center"; //the default is a centered window
	std::vector<double> vecArgs;
	parseWindowFilterArguments(filtername, i_vecArgs, 2, 2, isSoft, windowposition, vecArgs);

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
		const double& currentval = vecMeteo[ii](paramindex);
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
 * @brief Actual implementation of the smoothing algorithm
 * @param vecMyMeteo A vector of MeteoData
 * @param paramindex The meteo data parameter to be smoothed
 * @return The wma smoothing value for the specified meteo parameter of the last element in vecMeteo
 */
double FilterAlgorithms::WMASmoothingAlgorithm(const std::vector<MeteoData>& vecMyMeteo, const unsigned int& paramindex)
{
	vector<MeteoData> vecMeteo(vecMyMeteo);

	if (vecMeteo.size() == 0)
		return IOUtils::nodata;

	std::reverse(vecMeteo.begin(), vecMeteo.end()); //reverse the vector, posfind most significant

	double wma = IOUtils::nodata;
	unsigned int counter = 0;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		const double& currentval = vecMeteo[ii](paramindex);
		if (currentval != IOUtils::nodata){
			if (wma == IOUtils::nodata) wma = 0.0; //initialize wma here, otherwise it could be nodata

			wma += (vecMeteo.size()-counter) * currentval;
			counter++;
		}
	}

	if ((wma != IOUtils::nodata) && (counter > 0))
		wma /= (counter*(counter+1)/2);

	return wma;
}

/**
 * @brief Rate of change filter.
 * Calculate the change rate (ie: slope) between two points, if it is above a user given value, reject the point. Remarks:
 * - the maximum permissible rate of change (per seconds) has to be provided as an argument
 *
 * @code
 * TA::filter1	= rate
 * TA::arg1	= 0.01
 * @endcode
 */
void FilterAlgorithms::RateFilter(const std::vector<MeteoData>& /*vecM*/, const std::vector<std::string>& vecArgs,
                                  const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("rate", vecArgs, 1, 1, isSoft, doubleArgs);
	const double& maxRateOfChange = doubleArgs[0];

	int last_good=-1;
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++) {
		if(vecWindowM[ii](paramindex)!=IOUtils::nodata) {
			last_good=ii;
			break;
		}
	}
	if(last_good==-1)
		return; //can not find a good point to start
	const unsigned int start=last_good+1;

	for(unsigned int ii=start; ii<vecWindowM.size(); ii++) {
		double& value = vecWindowM[ii](paramindex);
		const double curr_time = vecWindowM[ii].date.getJulianDate();
		const double prev_value = vecWindowM[last_good](paramindex);
		const double prev_time = vecWindowM[last_good].date.getJulianDate();

		if(value==IOUtils::nodata) {
			continue;
		}
		const double local_rate = (value-prev_value) / ((curr_time-prev_time)*24.*3600.); //per seconds

		if( abs(local_rate) > maxRateOfChange ) {
			value = IOUtils::nodata;
		} else {
			last_good=ii;
		}
	}
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
void FilterAlgorithms::MinMaxFilter(const std::vector<MeteoData>& /*vecM*/, const std::vector<std::string>& vecArgs,
                                    const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("min_max", vecArgs, 2, 2, isSoft, doubleArgs);

	sort(doubleArgs.begin(), doubleArgs.end()); //the two parameters are sorted ascending

	//Run actual MinMax filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii](paramindex);

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
void FilterAlgorithms::MinValueFilter(const std::vector<MeteoData>& /*vecM*/, const std::vector<std::string>& vecArgs,
                                      const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("min", vecArgs, 1, 1, isSoft, doubleArgs);

	//Run actual MinValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii](paramindex);

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
void FilterAlgorithms::MaxValueFilter(const std::vector<MeteoData>& /*vecM*/, const std::vector<std::string>& vecArgs,
                                      const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> doubleArgs;
	parseFilterArguments("max", vecArgs, 1, 1, isSoft, doubleArgs);

	//Run actual MaxValue filter over all relevant meteo data
	for(unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii](paramindex);

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
void FilterAlgorithms::MedianAbsoluteDeviationFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                                     const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
//NOTE Problem: if a bunch of identical values ~median are part of the data set, mad will be 0
//and therefore we will reject all values that are != median
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii](paramindex);
		if (value == IOUtils::nodata) //No need to run filter for nodata points
			continue;

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
		} catch(const exception&){
			return;
		}

		double sigma = mad * K;

		if( (value>(median + 3.*sigma)) || (value<(median - 3.*sigma)) ) {
			value = IOUtils::nodata;
		}
	}
}


/**
 * @brief Standard deviation.
 * Values outside of mean ± 2 std_dev are rejected.
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = stddev
 *          TA::arg1    = soft left 1 1800  (1800 seconds time span for the left leaning window)
 *          RH::filter1 = stddev
 *          RH::arg1    = 10 6000            (strictly centered window spanning 6000 seconds and at least 10 points)
 * @endcode
 */
void FilterAlgorithms::StandardDeviationFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                               const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	for (unsigned int ii=0; ii<vecWindowM.size(); ii++){
		double& value = vecWindowM[ii](paramindex);
		if (value == IOUtils::nodata) //No need to run filter for nodata points
			continue;

		unsigned int pos = IOUtils::seek(vecWindowM[ii].date, vecM, false);

		std::vector<double> vecWindow;
		if (!getWindowData("stddev", vecM, pos, vecWindowM[ii].date, vecArgs, paramindex, vecWindow))
			return; //Not enough data to meet user configuration

		//Calculate deviation
		double mean     = IOUtils::nodata;
		double std_dev  = IOUtils::nodata;

		try {
			mean = Interpol1D::arithmeticMean(vecWindow);
			std_dev = Interpol1D::std_dev(vecWindow);
		} catch(const exception&){
			return;
		}

		if( abs(value-mean)>2.*std_dev) {
			value = IOUtils::nodata;
		}
	}
}

/**
 * @brief Tukey 53H method
 * A smooth time sequence is generated from the median, substracted from the original signal and compared with the standard deviation.
 * see <i>"Despiking Acoustic Doppler Velocimeter Data"</i>, Derek G. Goring and Vladimir L. Nikora, Journal of Hydraulic Engineering, <b>128</b>, 1, 2002
 * THIS CODE IS NOT ACTIVE YET
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = Tukey53H
 *          TA::arg1    = soft left 1 1800  (1800 seconds time span for the left leaning window)
 *          RH::filter1 = Tukey53H
 *          RH::arg1    = 10 6000            (strictly centered window spanning 6000 seconds and at least 10 points)
 * @endcode
 */
void FilterAlgorithms::Tukey53HFilter(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                      const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
	const double k = 3.;

	for (unsigned int ii=4; ii<vecWindowM.size()-4; ii++) {
		double& value = vecWindowM[ii](paramindex);
		if (value == IOUtils::nodata) //No need to run filter for nodata points
			continue;

		std::vector<double> vecWindow;
		unsigned int position = IOUtils::seek(vecWindowM[ii].date, vecM);

		if (!getWindowData("Tukey53H", vecM, position, vecWindowM[ii].date, vecArgs, paramindex, vecWindow))
			continue; //Not enough data to meet user configuration

		const double stddev = Interpol1D::std_dev(vecWindow);

		try {
			std::vector<double> u2;
			for(int i=-1; i<=1; i++) {
				std::vector<double> u1;
				for(int j=-1; j<=1; j++) {
					std::vector<double> u;
					for(int k=-2; k<=2; k++) {
						const int index = ii + k + j + i;
						const double value = vecWindowM[index](paramindex);
						if(value!=IOUtils::nodata)
							u.push_back( value );
					}
					if(u.size()>0)
						u1.push_back( Interpol1D::getMedian(u) );
				}
				if(u1.size()>0)
					u2.push_back( Interpol1D::getMedian(u1) );
				else
					u2.push_back( IOUtils::nodata );
			}
			//u3 = 1/4*( u2[0] + 2.*u2[1] + u2[2] )
			double u3=0.;
			unsigned int count=0;
			if(u2[0]!=IOUtils::nodata) {
				u3 += u2[0];
				count++;
			}
			if(u2[1]!=IOUtils::nodata) { //current timestep
				u3 += u2[1]*2.;
				count += 2;
			}
			if(u2[2]!=IOUtils::nodata) {
				u3 += u2[2];
				count++;
			}
			if(count>0) {
				u3 *= (double)(1/count);
				if( abs(value-u3) > k*stddev ) {
					value = IOUtils::nodata;
				}
			}
		} catch(const exception&){
			return;
		}

	}
}


/**
 * @brief Accumulation over a user given period.
 * The input data is accumulated over a given time interval (given as filter argument, in minutes).
 * This is for example needed for converting rain gauges measurements read every 10 minutes to
 * hourly precipitation measurements. Remarks: the accumulation period has to be provided as an argument (in seconds)
 * @code
 * HNW::filter1 = accumulate
 * HNW::arg1	 = 3600
 * @endcode
 */
void FilterAlgorithms::AccumulateProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                         const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
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
			bool exist=false;
			while((vecM[pos].date + deltatime) > vecM[startpos].date){
				const double& val = vecM[pos](paramindex);

				if (val != IOUtils::nodata) {
					sum += val;
					exist=true;
				}

				if (pos > 0) pos--;
				else break;
			}
			if(exist==true)
				vecWindowM[ii](paramindex) = sum;
			else
				vecWindowM[ii](paramindex) = IOUtils::nodata;
			//cout << "sum: " << vecWindowM[ii](paramindex) << endl;
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
void FilterAlgorithms::MedianAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                        const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
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
			vecWindowM[ii](paramindex) = median;
		} catch(const exception&){
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
void FilterAlgorithms::MeanAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                      const MeteoData::Parameters& paramindex, std::vector<MeteoData>& vecWindowM)
{
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
			vecWindowM[ii](paramindex) = mean;
		} catch(const exception&){
			continue; //the mean calculation did not work out, filter is not applied, value unchanged
		}
	}
}

/**
 * @brief Wind vector averaging.
 * This calculates the vector average over a user given time period. Each wind vector within this period
 * is added and the final sum is normalized by the number of vectors that have been added. Important information:
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
void FilterAlgorithms::WindAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<std::string>& vecArgs,
                                      const MeteoData::Parameters& /*paramindex*/, std::vector<MeteoData>& vecWindowM)
{
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
				ve += vecWindowVW[jj] * sin(vecWindowDW[jj] * Cst::PI / 180.); //turn into radians
				vn += vecWindowVW[jj] * cos(vecWindowDW[jj] * Cst::PI / 180.); //turn into radians
			}
			ve /= vecSize;
			vn /= vecSize;

			meanspeed = sqrt(ve*ve + vn*vn);
			meandirection = fmod( atan2(ve,vn) * 180. / Cst::PI + 360. , 360.); // turn into degrees [0;360)
		}

		vecWindowM[ii](MeteoData::VW) = meanspeed;
		vecWindowM[ii](MeteoData::DW) = meandirection;
	}
}

} //namespace

