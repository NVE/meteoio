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
#include <meteoio/ResamplingAlgorithms.h>
#include <cmath>
#include <algorithm>

using namespace std;

namespace mio {
 /**
 * @page resampling Resampling overview
 * The resampling infrastructure is described in ResamplingAlgorithms (for its API).
 * The goal of this page is to give an overview of the available resampling algorithms and their usage.
 *
 * @section resampling_section Resampling section
 * The resampling is specified for each parameter in the [Interpol1D] section. This section contains
 * a list of the various meteo parameters with their associated choice of resampling algorithm and
 * optional parameters. If a meteo parameter is not listed in this section, a linear resampling would be
 * assumed. An example of such section is given below:
 * @code
 * [Interpolations1D]
 * TA::resample    = linear
 *
 * RH::resample    = linear
 *
 * VW::resample    = nearest_neighbour
 * VW::args        = extrapolate
 *
 * HNW::resample   = linear
 * @endcode
 *
 * The size (in seconds) of the biggest gap that can be interpolated is given with the key WINDOW_SIZE. Therefore if two valid points are less than
 * WINDOW_SIZE seconds apart, points in between will be interpolated. If they are further apart, all points in between will remain IOUtils::nodata.
 * If using the "extrapolate" optional argument, points at WINDOW_SIZE distance of only one valid point will be extrapolated, otherwise they will remain
 * IOUtils::nodata. Please keep in mind that allowing extrapolated values can lead to grossly out of range data: using the slope
 * between two hourly measurements to extrapolate a point 10 days ahead is obviously risky!
 *
 * By default, WINDOW_SIZE is set to 10 days. This key has a <b>dramatic impact on run time/performance</b>:
 * reducing WINDOW_SIZE from 10 days down to one day makes geting meteo data 8 times faster while increasing it to 20 days makes it twice slower.
 *
 * @section algorithms_available Available Resampling Algorithms
 * Two algorithms for the resampling are implemented:
 * - none: do not perform resampling, see ResamplingAlgorithms::NoResampling
 * - linear: linear data resampling, see ResamplingAlgorithms::LinearResampling
 * - nearest_neighbour:  data resampling, see ResamplingAlgorithms::NearestNeighbour
 * - accumulate: data re-accumulation as suitable for precipitations, see ResamplingAlgorithms::Accumulate
 */

std::map<std::string, ResamplingAlgorithms::resamplingptr> ResamplingAlgorithms::algorithmMap;
const bool ResamplingAlgorithms::__init = ResamplingAlgorithms::initStaticData();

bool ResamplingAlgorithms::initStaticData()
{
	algorithmMap["none"]              = &ResamplingAlgorithms::NoResampling;
	algorithmMap["linear"]            = &ResamplingAlgorithms::LinearResampling;
	algorithmMap["nearest_neighbour"] = &ResamplingAlgorithms::NearestNeighbour;
	algorithmMap["accumulate"]        = &ResamplingAlgorithms::Accumulate;

	return true;
}

const ResamplingAlgorithms::resamplingptr& ResamplingAlgorithms::getAlgorithm(const std::string& algoname)
{
	std::map<std::string, resamplingptr>::const_iterator it;
	it = algorithmMap.find(algoname);

	if (it==algorithmMap.end())
		throw UnknownValueException("Unknown resampling algorithm called: " + algoname, AT);

	return it->second;
}

/**********************************************************************************
 * The following functions are implementations of different resampling algorithms *
 **********************************************************************************/

/**
 * @brief No resampling: do not resample parameter but keep original sampling rate
 * @code
 * [Interpolations1D]
 * TA::resample = none
 * @endcode
 */

void ResamplingAlgorithms::NoResampling(const size_t& /*pos*/, const ResamplingPosition& /*position*/, const size_t& /*paramindex*/,
                                        const std::vector<std::string>& /*taskargs*/, const double& /*window_size*/, const std::vector<MeteoData>& /*vecM*/, MeteoData& md)
{
	return;
}


/**
 * @brief Nearest Neighbour data resampling
 * Find the nearest neighbour of a desired data point that is not IOUtils::nodata and copy that value into the desired data point
 *        - If the data point itself is not IOUtils::nodata, nothing needs to be done
 *        - If two points have the same distance from the data point to be resampled, calculate mean and return it
 *        - if the argument extrapolate is provided, points within WINDOW_SIZE seconds of only one valid point will receive the value of this point
 * @code
 * [Interpolations1D]
 * TA::resample = nearest_neighbour
 * @endcode
 */

void ResamplingAlgorithms::NearestNeighbour(const size_t& index, const ResamplingPosition& position, const size_t& paramindex, 
                                            const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	const Date& resampling_date = md.date;

	if (position == ResamplingAlgorithms::exact_match) {
		const double& value = vecM[index](paramindex);
		if (value != IOUtils::nodata) { 
			md(paramindex) = value; //propagate value
			return;
		}
	}

	//Check whether extrapolation is activated
	bool extrapolate = false;
	if ((taskargs.size()==1) && (taskargs[0]=="extrapolate"))
		extrapolate = true;

	size_t indexP1=IOUtils::npos, indexP2=IOUtils::npos;
	getNearestValidPts(index, paramindex, vecM, resampling_date, window_size, indexP1, indexP2);
	const bool foundP1=(indexP1!=IOUtils::npos), foundP2=(indexP2!=IOUtils::npos);

	//Try to find the nearest neighbour, if there are two equally distant, then return the arithmetic mean
	if(foundP1 && foundP2) { //standard behavior
		const Duration diff1 = resampling_date - vecM[indexP1].date; //calculate time interval to element at index
		const Duration diff2 = vecM[indexP2].date - resampling_date; //calculate time interval to element at index
		const double& val1 = vecM[indexP1](paramindex);
		const double& val2 = vecM[indexP2](paramindex);

		if (IOUtils::checkEpsilonEquality(diff1.getJulianDate(true), diff2.getJulianDate(true), 0.1/1440)){ //within 6 seconds
			md(paramindex) = Interpol1D::weightedMean(val1, val2, 0.5);
		} else if (diff1 < diff2){
			md(paramindex) = val1;
		} else if (diff1 > diff2){
			md(paramindex) = val2;
		}
	} else if(!extrapolate) {
		return;
	} else if(extrapolate) {
		if(foundP1 && !foundP2){ //nearest neighbour on found after index 'index'
			md(paramindex) = vecM[indexP1](paramindex);
		} else if (!foundP1 && foundP2){ //nearest neighbour on found before index 'index'
			md(paramindex) = vecM[indexP2](paramindex);
		} else { // no nearest neighbour with a value different from IOUtils::nodata
			return;
		}
	}
}

/**
 * @brief Linear data resampling: If a point is requested that is in between two input data points,
 *        the requested value is automatically calculated using a linear interpolation. Furthermore
 *        if the argument extrapolate is provided there will be an attempt made to extrapolate the
 *        point if the interpolation fails, by solving the line equation y = kx + d
 * @code
 * [Interpolations1D]
 * TA::resample = linear
 * TA::args     = extrapolate
 * @endcode
 */
void ResamplingAlgorithms::LinearResampling(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                                            const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	const Date& resampling_date = md.date;

	if (position == ResamplingAlgorithms::exact_match) {
		const double& value = vecM[index](paramindex);
		if (value != IOUtils::nodata) { 
			md(paramindex) = value; //propagate value
			return;
		}
	}

	//Check whether extrapolation is activated
	bool extrapolate = false;
	if ((taskargs.size()==1) && (taskargs[0]=="extrapolate"))
		extrapolate = true;

	if ((!extrapolate) && position == ResamplingAlgorithms::end) //no extrapolation wished for
		return;

	size_t indexP1=IOUtils::npos, indexP2=IOUtils::npos;

	getNearestValidPts(index, paramindex, vecM, resampling_date, window_size, indexP1, indexP2);
	bool foundP1=(indexP1!=IOUtils::npos), foundP2=(indexP2!=IOUtils::npos);

	//do nothing if we can't interpolate, and extrapolation is not explicitly activated
	if ((!extrapolate) && ((!foundP1) || (!foundP2)))
		return;
	//do nothing if not at least one value different from IOUtils::nodata has been found
	if (!foundP1 && !foundP2)
		return;

	//At this point we either have a valid indexP1 or indexP2 and we can at least try to extrapolate
	if (!foundP1 && foundP2){ //only nodata values found before index, try looking after indexP2
		for (size_t ii=indexP2+1; ii<vecM.size(); ii++){
			if (vecM[ii](paramindex) != IOUtils::nodata){
				indexP1 = ii;
				foundP1 = true;
				break;
			}
		}
	} else if (foundP1 && !foundP2){ //only nodata found after index, try looking before indexP1
		for (size_t ii=indexP1; (ii--) > 0; ){
			if (vecM[ii](paramindex) != IOUtils::nodata){
				indexP2=ii;
				foundP2 = true;
				break;
			}
		}
	}

	if (!foundP1 || !foundP2) //now at least two points need to be present
		return;

	//At this point indexP1 and indexP2 point to values that are different from IOUtils::nodata
	const double& val1 = vecM[indexP1](paramindex);
	const double jul1 = vecM[indexP1].date.getJulianDate(true);
	const double& val2 = vecM[indexP2](paramindex);
	const double jul2 = vecM[indexP2].date.getJulianDate(true);

	/*
	cout << "Resampling " << resampling_date.toString(Date::ISO) << " param " << md.getParameterName(paramindex) << endl
		<< "P1: " << vecM[indexP1].date.toString(Date::ISO) << "  val: " << val1 << endl
		<< "P1: " << vecM[indexP2].date.toString(Date::ISO) << "  val: " << val2 << endl;
	*/
	md(paramindex) = linearInterpolation(jul1, val1, jul2, val2, resampling_date.getJulianDate(true));
}

/**
 * @brief Accumulation over a user given period.
 * The input data is accumulated over a given time interval (given as filter argument, in seconds).
 * This is for example needed for converting rain gauges measurements read every 10 minutes to
 * hourly precipitation measurements. Remarks:
 * - the accumulation period has to be provided as an argument (in seconds)
 * @code
 * HNW::filter1 = accumulate
 * HNW::arg1	 = 3600
 * @endcode
 */
void ResamplingAlgorithms::Accumulate(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                                            const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size())
		throw IOException("The index of the element to be resampled is out of bounds", AT);

	const Date& resampling_date = md.date;

	md(paramindex) = IOUtils::nodata;

	//Read accumulation period
	double accumulate_period;
	if (taskargs.size()==1) {
		IOUtils::convertString(accumulate_period, taskargs[0]);
		if(accumulate_period<=0.) {
			std::stringstream tmp;
			tmp << "Invalid accumulation period (" << accumulate_period << ") ";
			tmp << "for parameter " << vecM.at(0).getNameForParameter(paramindex);
			throw InvalidArgumentException(tmp.str(), AT);
		}
	} else {
		std::stringstream tmp;
		tmp << "Please provide accumulation period (in seconds) for param" << vecM.at(0).getNameForParameter(paramindex);
		throw InvalidArgumentException(tmp.str(), AT);
	}

	//find start of accumulation period
	const Date dateStart(resampling_date.getJulianDate() - accumulate_period/(24.*3600.), resampling_date.getTimeZone());
	bool found_start=false;
	size_t start_idx; //this is the index of the first point of the window that will contain dateStart
	for (start_idx=index+1; (start_idx--) > 0; ) {
		if(vecM[start_idx].date <= dateStart) {
			found_start=true;
			break;
		}
	}
	if (!found_start) {
		cerr << "[W] Could not accumulate " << vecM.at(0).getNameForParameter(paramindex) << ": ";
		cerr << "not enough data for accumulation period at date " << resampling_date.toString(Date::ISO) << "\n";
		md(paramindex) = IOUtils::nodata;
		return;
	}

	//resample the starting point
	double sum = funcval(start_idx, paramindex, vecM, dateStart); //resampling the starting point
	if(sum==IOUtils::nodata) return;

	 //sum all data points until current position
	for(size_t idx=(start_idx+1); idx<index; idx++) { //HACK: <= on filtered_vec
		if(vecM[idx](paramindex)==IOUtils::nodata) return;
		sum += vecM[idx](paramindex);
	}

	//HACK resample end point (from filtered_vec)
	md(paramindex) = sum;
}

double ResamplingAlgorithms::funcval(const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM,
                                     const Date& date)
{
	size_t start = pos;
	if (vecM[start].isResampled()) { //HACK how could a point in vecM be resampled?
		if (start > 0){
			start--;
		} else {
			return IOUtils::nodata;
		}
	}
	const double& valstart = vecM[start](paramindex);

	if (!vecM[pos].isResampled() && (vecM[pos].date == date))
		return valstart;

	size_t end = pos+1;
	if ((vecM[end].isResampled())) //skip resampled value
		end++;

	//index either points to the element with date directly or to the last element with vecM[index] <= date
	if (end < vecM.size()){
		const double& valend = vecM[end](paramindex);
		if (valend == IOUtils::nodata)
			return IOUtils::nodata;

		const double jul1 = vecM[start].date.getJulianDate(true);
		const double jul2 = vecM[end].date.getJulianDate(true);

		return linearInterpolation(jul1, 0., jul2, valend, date.getJulianDate(true));
	}

	return IOUtils::nodata;
}

/**
 * @brief This function returns the last and next valid points around a given position
 * @param pos current position (index)
 * @param paramindex meteo parameter to use
 * @param vecM vector of MeteoData
 * @param window_size size of the search window
 * @param indexP1 index of point before the current position (IOUtils::npos if none could be found)
 * @param indexP2 index of point after the current position (IOUtils::npos if none could be found)
 */
void ResamplingAlgorithms::getNearestValidPts(const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& window_size, size_t& indexP1, size_t& indexP2)
{
	indexP1=IOUtils::npos;
	indexP2=IOUtils::npos;

	const Date dateStart = resampling_date - window_size;
	for (size_t ii=pos; (ii--) > 0; ) {
		if (vecM[ii].date < dateStart) break;
		if (vecM[ii](paramindex) != IOUtils::nodata){
			indexP1 = ii;
			break;
		}
	}

	Date dateEnd;
	if (indexP1 != IOUtils::npos) dateEnd = vecM[indexP1].date + window_size; //so the search window remains window_size
	else dateEnd = resampling_date + window_size;

	for (size_t ii=pos; ii<vecM.size(); ii++) {
		if (vecM[ii].date > dateEnd) break;
		if (vecM[ii](paramindex) != IOUtils::nodata) {
			indexP2 = ii;
			break;
		}
	}
}

/**
 * @brief This function solves the equation y = ax + b for two given points and returns y for a given x
 * @param x1 x-coordinate of first point
 * @param y1 y-coordinate of first point
 * @param x2 x-coordinate of second point
 * @param y2 y-coordinate of second point
 * @param x3 x-coordinate of desired point
 * @return y-coordinate of desired point
 */
double ResamplingAlgorithms::linearInterpolation(const double& x1, const double& y1,
                                       const double& x2, const double& y2, const double& x3)
{
	if (x1 == x2)
		throw IOException("Attempted division by zero", AT);

	//Solving y = ax + b
	const double a = (y2 - y1) / (x2 - x1);
	const double b = y2 - a*x2;

	return (a*x3 + b);
}

} //namespace

