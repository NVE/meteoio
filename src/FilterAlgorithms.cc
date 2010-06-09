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
 *
 * Two interpolation mechanism used for the resampling are implemented:
 * - linear: linear data resampling, see FilterAlgorithms::LinResamplingProcess
 * - nearest_neighbour:  data resampling, see FilterAlgorithms::NearestNeighbourResamplingProcess
 */

std::map<std::string, FilterProperties> FilterAlgorithms::filterMap;
const bool FilterAlgorithms::__init = FilterAlgorithms::initStaticData();

bool FilterAlgorithms::initStaticData()
{
	filterMap["rate"]     = FilterProperties(true, (unsigned int)1, Date(0.0), &FilterAlgorithms::RateFilter);
	filterMap["min_max"]  = FilterProperties(true, (unsigned int)1, Date(0.0), &FilterAlgorithms::MinMaxFilter);
	filterMap["min"]      = FilterProperties(true, (unsigned int)1, Date(0.0), &FilterAlgorithms::MinValueFilter);
	filterMap["max"]      = FilterProperties(true, (unsigned int)1, Date(0.0), &FilterAlgorithms::MaxValueFilter);
	filterMap["mad"]      = FilterProperties(true, (unsigned int)1, Date(0.0), &FilterAlgorithms::MedianAbsoluteDeviationFilter);
	filterMap["accumulate"] = FilterProperties(false, (unsigned int)1, Date(0.0), &FilterAlgorithms::AccumulateProcess);
	filterMap["linear"] = FilterProperties(false, (unsigned int)1, Date(0.0), &FilterAlgorithms::LinResamplingProcess);
	filterMap["nearest_neighbour"] = FilterProperties(false, (unsigned int)1, Date(0.0), &FilterAlgorithms::NearestNeighbourResamplingProcess);
	filterMap["median_avg"] = FilterProperties(false, (unsigned int)1, Date(0.0), &FilterAlgorithms::MedianAvgProcess);
	filterMap["mean_avg"] = FilterProperties(false, (unsigned int)1, Date(0.0), &FilterAlgorithms::MeanAvgProcess);
	filterMap["wind_avg"] = FilterProperties(false, (unsigned int)1, Date(0.0), &FilterAlgorithms::WindAvgProcess);

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
 * @brief Rate of change filter.
 * Calculate the change rate (ie: slope) between two points, if it is above a user given value, reject the point. Remarks:
 * - the maximum permissible rate of change (per seconds) has to be provided as an argument
 * @code
 * TA::filter1	= rate
 * TA::arg1	= 0.01
 * @endcode
 */
bool FilterAlgorithms::RateFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
						   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
						   const unsigned int& paramindex,
						   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)date; (void)vecFilteredS;
	(void)_vecArgs; (void)paramindex; (void)vecFilteredM;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs;
	parseFilterArguments("rate", _vecArgs, 1, 1, isSoft, vecArgs);

	//Run actual Rate filter over all relevant meteo data
	if(vecFilteredM.size()>1) {
		//the request time step is NOT part of the data, it will have to be resampled
		const double prev_value = vecM[pos-vecFilteredM.size()].param(paramindex);
		const double prev_time = vecM[pos-vecFilteredM.size()].date.getJulianDate()*24.*3600.; //in seconds

		for(unsigned int ii=0; ii<vecFilteredM.size(); ii++){
			double& value = vecFilteredM[ii].param(paramindex);
			const double curr_value = vecFilteredM[ii].param(paramindex);
			const double curr_time = vecFilteredM[ii].date.getJulianDate()*24.*3600.; //in seconds
			const double local_rate = abs((curr_value-prev_value)/(curr_time-prev_time));

			if( local_rate > vecArgs[0] ) {
				value = IOUtils::nodata;
			}
		}
	} else {
		//const double tmp = vecFilteredM.size();
		//std::cout << tmp << std::endl;
		//the request time step is part of the data
		if(pos>1) { //if we are at the start of the data set, we can not apply the filter...
			double& value = vecFilteredM[0].param(paramindex);
			const double curr_value = vecM[pos].param(paramindex);
			const double curr_time = vecM[pos].date.getJulianDate()*24.*3600.; //in seconds
			const double prev_value = vecM[pos-1].param(paramindex);
			const double prev_time = vecM[pos-1].date.getJulianDate()*24.*3600.; //in seconds
			const double local_rate = abs((curr_value-prev_value)/(curr_time-prev_time));

			if( local_rate > vecArgs[0] ) {
				value = IOUtils::nodata;
			}
		} else {
			//the filter can not be applied
			return false;
		}
	}

	return true;
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
bool FilterAlgorithms::MinMaxFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
						 const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
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
bool FilterAlgorithms::MinValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
						   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
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
bool FilterAlgorithms::MaxValueFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
						   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
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
bool FilterAlgorithms::MedianAbsoluteDeviationFilter(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
				   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex,
				   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecS; (void)vecFilteredS;

	std::vector<double> vecWindow;
	if (!getWindowData("mad", vecM, pos, date, _vecArgs, paramindex, vecWindow))
		return false; //Not enough data to meet user configuration
	
	//Calculate MAD
	const double K = 1. / 0.6745;
	double mad = IOUtils::nodata;
	double median = IOUtils::nodata;

	try {
		median = Interpol1D::getMedian(vecWindow);
		mad = Interpol1D::getMedianAverageDeviation(vecWindow);
	} catch(exception& e){
		return false;
	}
	
	double sigma = mad * K;

	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		double& value = vecFilteredM[ii].param(paramindex);
		if( (value>(median + 3.*sigma)) || (value<(median - 3.*sigma)) ) {
			value = IOUtils::nodata;
		}
	}

	return true;
}


/**
 * @brief Accumulation over a user given period.
 * The input data is accumulated over a given time interval (given as filter argument, in minutes).
 * This is for example needed for converting rain gauges measurements read every 10 minutes to
 * hourly precipitation measurements. Remarks:
 * - the accumulation period has to be provided as an argument (in seconds)
 * @code
 * HNW::filter1	= accumulate
 * HNW::arg1	= 3600
 * @endcode
 */
bool FilterAlgorithms::AccumulateProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
						    const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
						    const unsigned int& paramindex,
						    std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)date; (void)vecS; (void)vecFilteredS;
	//parse arguments and check whether they are valid
	bool isSoft = false;
	std::vector<double> vecArgs;
	parseFilterArguments("accumulate", _vecArgs, 1, 1, isSoft, vecArgs);

	Date deltatime(vecArgs[0]/(24.*3600.)); //making a julian date out of the argument given in seconds

	unsigned int index = vecFilteredM.size();
	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		if (vecFilteredM[ii].date >= date){
			index = ii;
			break;
		}
	}

	if (index >= vecFilteredM.size())
	  return false;

	int startposition = pos;
	for (int ii=index; (ii>=((int)index-1) && (ii>=0)); ii--){
		unsigned int mypos = startposition;

		double sum = 0.0;
		while((vecM[mypos].date + deltatime) > vecM[startposition].date){
			const double& val = vecM[mypos].param(paramindex);
			if (val != IOUtils::nodata)
				sum += vecM[mypos].param(paramindex);
			//cout << vecM[mypos].date.toString(Date::ISO) << " HNW:" << vecM[mypos].param(paramindex) << "  SUM:" << sum << endl;		
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

/**
 * @brief Nearest Neighbour data resampling: Find the nearest neighbour of a desired data point 
 *        that is not IOUtils::nodata and copy that value into the desired data point
 *        - If the data point itself is not IOUtils::nodata, nothing needs to be done
 *        - If two points have the same distance from the data point to be resampled, calculate mean and return it
 *        - no arguments are considered
 * @code
 * [Interpolations1D]
 * TA::resample = nearest_neighbour
 * @endcode
 */
bool FilterAlgorithms::NearestNeighbourResamplingProcess(const std::vector<MeteoData>& vecM, 
                                                           const std::vector<StationData>& vecS,
                                                           const unsigned int& pos, const Date& date, 
                                                           const std::vector<std::string>& _vecArgs,
                                                           const unsigned int& paramindex,
                                                           std::vector<MeteoData>& vecFilteredM, 
                                                           std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos; (void)_vecArgs;
	int indexBefore=-1, indexExact=-1, indexAfter=-1;

	if ((vecFilteredM.size()==1) && (date==vecFilteredM[0].date)){//Nothing to do
		return false; //Interpretation: filter not applied
	} else if (vecFilteredM.size() == 0){
		throw IOException("Not enough data to do resampling ...", AT);
	} else if (vecFilteredM.size()>=2){
		//check whether there is already an element with the correct date in vecFilteredM
		for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
			if (vecFilteredM.at(ii).date < date){
				indexBefore = (int)ii;
			} else if (vecFilteredM.at(ii).date == date){
				indexExact = (int)ii;
				if (indexExact < (int)(vecFilteredM.size()-1))
					indexAfter = indexExact + 1;
				break;
			} else if (vecFilteredM.at(ii).date > date){
				indexAfter = (int)ii;
				break;
			}
		}

		if ((indexExact == -1) && (indexAfter>0)){ //no MeteoData object with the correct date present
			MeteoData newmd(date);
			newmd.setResampled(true);

			vecFilteredM.insert(vecFilteredM.begin() + indexAfter, newmd);
			vecFilteredS.insert(vecFilteredS.begin() + indexAfter, vecFilteredS[0]);
			indexExact = indexAfter;
			indexAfter++;
		} else if ((indexExact == -1) && (indexAfter == -1) && (indexBefore>=0)){
			MeteoData newmd(date);
			newmd.setResampled(true);

			vecFilteredM.push_back(newmd);
			vecFilteredS.push_back(vecFilteredS.at(0));
			indexExact = vecFilteredM.size() - 1;
		} else if ((indexExact == -1) && (indexBefore == -1) && (indexAfter>=0)){
			MeteoData newmd(date);
			newmd.setResampled(true);

			vecFilteredM.insert(vecFilteredM.begin(), newmd);
			vecFilteredS.insert(vecFilteredS.begin(), vecFilteredS[indexAfter]);
			indexExact = 0;
			indexAfter++;
		} else if (indexExact != -1){
		  if (vecFilteredM[indexExact].param(paramindex) != IOUtils::nodata)
		    return false;
		}
	}

	//Try to find the nearest neighbour, if there are two equally distant, then return the arithmetic mean
	MeteoData m1, m2;
	bool found1=false, found2=false;
	for (unsigned int ii=indexExact+1; ii<vecFilteredM.size(); ii++){
		if (vecFilteredM[ii].param(paramindex) != IOUtils::nodata){
			m1 = vecFilteredM[ii];
			found1 = true;
			break;
		}
	}

	for (unsigned int ii=0; ii<(unsigned int)indexExact; ii++){
		if (vecFilteredM[ii].param(paramindex) != IOUtils::nodata){
			m2 = vecFilteredM[ii];
			found2 = true;
		}
	}

	if (found1 && !found2){
		vecFilteredM[indexExact].param(paramindex) = m1.param(paramindex);
	} else if (!found1 && found2){
		vecFilteredM[indexExact].param(paramindex) = m2.param(paramindex);
	} else if (!found1 && !found2){
		vecFilteredM[indexExact].param(paramindex) = IOUtils::nodata;
	} else {
		Date diff1 = m1.date - vecFilteredM[indexExact].date;
		Date diff2 = vecFilteredM[indexExact].date - m2.date;

		if (IOUtils::checkEpsilonEquality(diff1.getJulianDate(), diff2.getJulianDate(), 0.1/1440)){ //within 6 seconds
			vecFilteredM[indexExact].param(paramindex) = Interpol1D::linearInterpolation(m1.param(paramindex), 
																		  m2.param(paramindex), 0.5);
		} else if (diff1 < diff2){
			vecFilteredM[indexExact].param(paramindex) = m1.param(paramindex);
		} else if (diff1 > diff2){
			vecFilteredM[indexExact].param(paramindex) = m2.param(paramindex);
		}
	}

	return true;
}

/**
 * @brief Linear data resampling: If a point is requested that is in between two input data points, 
 *        the requested value is automatically calculated using a linear interpolation. Furthermore 
 *        if the argument extrapolate is provided there will be an attempt made to extrapolate the
 *        point if the interpolation fails
 * @code
 * [Interpolations1D]
 * TA::resample = linear
 * TA::args     = extrapolate
 * @endcode
 */
bool FilterAlgorithms::LinResamplingProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
							const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
							const unsigned int& paramindex,
							std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecM; (void)vecS; (void)pos;
	int indexBefore=-1, indexExact=-1, indexAfter=-1;
	bool extrapolate = false;

	//Check whether extrapolation is desired
	if (_vecArgs.size() > 0)
		if (_vecArgs.at(0) == "extrapolate")
			extrapolate = true;
	

	if ((vecFilteredM.size()==1) &&(date==vecFilteredM[0].date)){//Nothing to do
		return false; //Interpretation: filter not applied
	} else if (vecFilteredM.size() == 0){
		throw IOException("Not enough data to do resampling ...", AT);
	} else if (vecFilteredM.size()>=2){
		//check whether there is already an element with the correct date in vecFilteredM
		for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
			if (vecFilteredM.at(ii).date < date){
				indexBefore = (int)ii;
			} else if (vecFilteredM.at(ii).date == date){
				indexExact = (int)ii;
				if (indexExact < (int)(vecFilteredM.size()-1))
					indexAfter = indexExact + 1;
				break;
			} else if (vecFilteredM.at(ii).date > date){
				indexAfter = (int)ii;
				break;
			}
		}
	}

	//Now check whether a new element needs to be inserted and whether resampling is possible
	//cout <<"Size of vector: " << vecFilteredM.size()<<endl;
	//cout << "before="<<indexBefore<< "  exact=" << indexExact << "  after="<<indexAfter<< endl;

	if (indexExact == -1){ //no MeteoData object with the correct date present -> insertion
		MeteoData newmd(date);
		newmd.setResampled(true);

		if (indexAfter >= 0){
			vecFilteredM.insert(vecFilteredM.begin() + indexAfter, newmd);
			vecFilteredS.insert(vecFilteredS.begin() + indexAfter, vecFilteredS[0]);
			indexExact = indexAfter;
			indexAfter++;
		} else if (indexBefore >=0) {
			vecFilteredM.insert(vecFilteredM.begin() + indexBefore + 1, newmd);
			vecFilteredS.insert(vecFilteredS.begin() + indexBefore + 1, vecFilteredS[0]);
			indexExact = indexBefore + 1;
		}
	} else if (((indexExact != -1) && (indexAfter == -1)) || ((indexExact != -1) && (indexBefore == -1))){
		if (!extrapolate)
			return false; //No resampling possible, for lack of a right/left element
	} else if (indexExact != -1){
		if (vecFilteredM[indexExact].param(paramindex) != IOUtils::nodata)
			return false;
	}

	//Now find two points within the vecFilteredM (before and aft, that are not IOUtils::nodata)
	//If that condition cannot be met, simply add nodata for the resampled value

	bool found1=false, found2=false;
	for (unsigned int ii=indexBefore+1; (ii--) > 0; ){
		if (vecFilteredM[ii].param(paramindex) != IOUtils::nodata){
			indexBefore=ii;
			found1 = true;
			break;
		}
	}		

	for (unsigned int ii=indexAfter; ii<vecFilteredM.size(); ii++){
		if (vecFilteredM[ii].param(paramindex) != IOUtils::nodata){
			indexAfter = ii;
			found2 = true;
			break;
		}
	}

	if (extrapolate){
		if (!found1 && found2){ //only nodata values found before indexExact, try looking after indexAfter
			for (unsigned int ii=indexAfter; ii<vecFilteredM.size(); ii++){
				if (vecFilteredM[ii].param(paramindex) != IOUtils::nodata){
					indexBefore = ii;
					found1 = true;
					break;
				}
			}		
		} else if (found1 && !found2){ //only nodata found after indexExact, try looking before indexBefore
			for (unsigned int ii=indexBefore+1; (ii--) > 0; ){
				if (vecFilteredM[ii].param(paramindex) != IOUtils::nodata){
					indexAfter=ii;
					found2 = true;
					break;
				}
			}		
		}
	}

	MeteoData& resampledmd  = vecFilteredM.at((unsigned int) indexExact);
	const MeteoData& tmpmd1 = vecFilteredM.at((unsigned int) indexBefore);
	const MeteoData& tmpmd2 = vecFilteredM.at((unsigned int) indexAfter);

	const double& val1 = tmpmd1.param(paramindex);
	const double& val2 = tmpmd2.param(paramindex);

	if ((val1 == IOUtils::nodata) || (val2 == IOUtils::nodata)){
		resampledmd.param(paramindex) = IOUtils::nodata;
	} else {
		resampledmd.param(paramindex) = Interpol1D::linearInterpolation(tmpmd1.date.getJulianDate(), val1,
														    tmpmd2.date.getJulianDate(), val2, 
														    date.getJulianDate());
	}

	return true;
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
bool FilterAlgorithms::MedianAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
				   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
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

	double median = IOUtils::nodata;

	try {
		median = Interpol1D::getMedian(vecWindow);
	} catch(exception& e){
		return false; //the median calculation did not work out, filter is not applied
	}

	//cout << "Median value: " << median << endl;
	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		vecFilteredM[ii].param(paramindex) = median;
	}

	return true;
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
bool FilterAlgorithms::MeanAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
				   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
				   const unsigned int& paramindex,
				   std::vector<MeteoData>& vecFilteredM, std::vector<StationData>& vecFilteredS)
{
	(void)vecS; (void)vecFilteredS;
	//Remarks:
	//1) nodata values are not considered when calculating the mean
	//2) Two arguments expected (both have to be fullfilled for the filter to start operating):
	//   1. minimal number of points in window
	//   2. minimal time interval spanning the window
	//3) the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window
	//   position
	//4) the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
	//
	std::vector<double> vecWindow;
	if (!getWindowData("mean_avg", vecM, pos, date, _vecArgs, paramindex, vecWindow))
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
bool FilterAlgorithms::WindAvgProcess(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS,
				   const unsigned int& pos, const Date& date, const std::vector<std::string>& _vecArgs,
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
	std::vector<Date> vecDateVW, vecDateDW;
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
			ve += vecWindowVW[ii] * sin(vecWindowDW[ii] * M_PI / 180.); //turn into radians
			vn += vecWindowVW[ii] * cos(vecWindowDW[ii] * M_PI / 180.); //turn into radians
		}
		ve /= vecSize;
		vn /= vecSize;

		meanspeed = sqrt(ve*ve + vn*vn);
		meandirection = fmod( atan2(ve,vn) * 180. / M_PI + 360. , 360.); // turn into degrees [0;360)
	}

	for (unsigned int ii=0; ii<vecFilteredM.size(); ii++){
		vecFilteredM[ii].param(MeteoData::VW) = meanspeed;
		vecFilteredM[ii].param(MeteoData::DW) = meandirection;
	}

	return true;
}

} //namespace

