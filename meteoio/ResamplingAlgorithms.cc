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
 * The search window size can be given with key WINDOW_SIZE that expresses (in seconds) how far a valid point can be searched for when
 * re-interpolating a missing value (up to WINDOW_SIZE/2 before and after the requested point, ie the search window is centered).
 *
 * @section algorithms_available Available Resampling Algorithms
 * Two algorithms for the resampling are implemented:
 * - none: do not perform resampling, see ResamplingAlgorithms::NoResampling
 * - linear: linear data resampling, see ResamplingAlgorithms::LinearResampling
 * - nearest_neighbour:  data resampling, see ResamplingAlgorithms::NearestNeighbour
 * - accumulate: data re-accumulation as suitable for precipitations, see ResamplingAlgorithms::Accumulate
 */

std::map<std::string, resamplingptr> ResamplingAlgorithms::algorithmMap;
const bool ResamplingAlgorithms::__init = ResamplingAlgorithms::initStaticData();

bool ResamplingAlgorithms::initStaticData()
{
	algorithmMap["none"]              = &ResamplingAlgorithms::NoResampling;
	algorithmMap["linear"]            = &ResamplingAlgorithms::LinearResampling;
	algorithmMap["nearest_neighbour"] = &ResamplingAlgorithms::NearestNeighbour;
	algorithmMap["accumulate"]        = &ResamplingAlgorithms::Accumulate;

	return true;
}

const resamplingptr& ResamplingAlgorithms::getAlgorithm(const std::string& algoname)
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

void ResamplingAlgorithms::NoResampling(const size_t& /*pos*/, const size_t& /*paramindex*/,
                                        const std::vector<std::string>& /*taskargs*/, std::vector<MeteoData>& /*vecM*/)
{
	return;
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

void ResamplingAlgorithms::NearestNeighbour(const size_t& pos, const size_t& paramindex,
                                            const std::vector<std::string>& /*taskargs*/, std::vector<MeteoData>& vecM)
{
	if (pos >= vecM.size())
		throw IOException("The position of the resampled element is out of bounds", AT);


	if (vecM[pos](paramindex) != IOUtils::nodata) //if there is a value, don't resample
		return;

	//Try to find the nearest neighbour, if there are two equally distant, then return the arithmetic mean
	MeteoData m1, m2;
	bool found1=false, found2=false;
	for (size_t ii=pos+1; ii<vecM.size(); ii++){
		if (vecM[ii](paramindex) != IOUtils::nodata){
			m1 = vecM[ii];
			found1 = true;
			break;
		}
	}

	for (size_t ii=0; ii<pos; ii++){
		if (vecM[ii](paramindex) != IOUtils::nodata){
			m2 = vecM[ii];
			found2 = true;
		}
	}

	if (found1 && !found2){ //nearest neighbour on found after index 'pos'
		vecM[pos](paramindex) = m1(paramindex);
	} else if (!found1 && found2){ //nearest neighbour on found before index 'pos'
		vecM[pos](paramindex) = m2(paramindex);
	} else if (!found1 && !found2){ // no nearest neighbour with a value different from IOUtils::nodata
		vecM[pos](paramindex) = IOUtils::nodata;
	} else {
		Duration diff1 = m1.date - vecM[pos].date; //calculate time interval to element at pos
		Duration diff2 = vecM[pos].date - m2.date; //calculate time interval to element at pos

		if (IOUtils::checkEpsilonEquality(diff1.getJulianDate(true), diff2.getJulianDate(true), 0.1/1440)){ //within 6 seconds
			vecM[pos](paramindex) = Interpol1D::weightedMean(m1(paramindex), m2(paramindex), 0.5);
		} else if (diff1 < diff2){
			vecM[pos](paramindex) = m1(paramindex);
		} else if (diff1 > diff2){
			vecM[pos](paramindex) = m2(paramindex);
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
void ResamplingAlgorithms::LinearResampling(const size_t& pos, const size_t& paramindex,
                                            const std::vector<std::string>& taskargs, std::vector<MeteoData>& vecM)
{
	if (pos >= vecM.size())
		throw IOException("The position of the resampled element is out of bounds", AT);

	if (vecM[pos](paramindex) != IOUtils::nodata) //if there is a value, don't resample
		return;

	//Check whether extrapolation is activated
	bool extrapolate = false;
	if ((taskargs.size()==1) && (taskargs[0]=="extrapolate"))
		extrapolate = true;

	//Now find two points within the vecM (before and aft, that are not IOUtils::nodata)
	//If that condition cannot be met, simply add nodata for the resampled value (exception: extrapolate)
	size_t indexP1=IOUtils::npos, indexP2=IOUtils::npos;
	bool foundP1=false, foundP2=false;

	for (size_t ii=pos; (ii--) > 0; ){
		if (vecM[ii](paramindex) != IOUtils::nodata){
			indexP1 = (size_t)ii;
			foundP1 = true;
			break;
		}
	}

	for (size_t ii=pos+1; ii<vecM.size(); ii++){
		if (vecM[ii](paramindex) != IOUtils::nodata){
			indexP2 = (size_t)ii;
			foundP2 = true;
			break;
		}
	}

	//do nothing if we can't interpolate, and extrapolation is not explicitly activated
	if ((!extrapolate) && ((!foundP1) || (!foundP2)))
		return;

	//do nothing if not at least one value different from IOUtils::nodata has been found
	if (!foundP1 && !foundP2)
		return;

	//At this point we either have a valid indexP1 or indexP2 and we can at least try to extrapolate
	if (!foundP1 && foundP2){ //only nodata values found before pos, try looking after indexP2
		for (size_t ii=indexP2+1; ii<vecM.size(); ii++){
			if (vecM[ii](paramindex) != IOUtils::nodata){
				indexP1 = ii;
				foundP1 = true;
				break;
			}
		}
	} else if (foundP1 && !foundP2){ //only nodata found after pos, try looking before indexP1
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
	const double& val2 = vecM[indexP2](paramindex);
	vecM[pos](paramindex) = linearInterpolation(vecM[indexP1].date.getJulianDate(true), val1,
	                                                  vecM[indexP2].date.getJulianDate(true), val2,
	                                                  vecM[pos].date.getJulianDate(true));
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
void ResamplingAlgorithms::Accumulate(const size_t& pos, const size_t& paramindex,
                                      const std::vector<std::string>& taskargs, std::vector<MeteoData>& vecM)
{
	/*
	 * HACK TODO: Overall check IOUtils::nodata data path and test all scenarios with good test cases
	 */

	if (pos >= vecM.size())
		throw IOException("The position of the resampled element is out of bounds", AT);

	//Get accumulation period
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
	bool found_start=false;
	size_t start_idx = pos+1;
	Date dateStart(vecM[pos].date.getJulianDate() - accumulate_period/(24.*3600.), vecM[pos].date.getTimeZone());

	for (start_idx=pos+1; (start_idx--) > 0; ){
		if(vecM[start_idx].date <= dateStart) {
			found_start=true;
			break;
		}
	}

	if (!found_start){
		cerr << "[W] Could not accumulate " << vecM.at(0).getNameForParameter(paramindex)
			<< ", not enough data for accumulation period at date " << vecM[pos].date.toString(Date::ISO)
			<< endl;
		vecM[pos](paramindex) = IOUtils::nodata;
		return;
	}

	//resample the starting point
	//HACK: we consider nodata to be 0. In fact, we should try to interpolate from valid points
	//if they are not too far away

	size_t interval_end   = start_idx + 1;

	double valstart = funcval(vecM, start_idx, dateStart, paramindex);
	double valend   = funcval(vecM, interval_end, vecM[interval_end].date, paramindex);
	double sum = IOUtils::nodata;

	if ((valend == IOUtils::nodata) || (valstart == IOUtils::nodata)){
		sum = 0.0; //HACK maybe it should be set it to IOUtils::nodata
	} else {
		if ((start_idx == (pos-1)) && (dateStart == vecM[start_idx].date))
			valstart = 0.0;

		sum = valend - valstart;
	}

	if (interval_end == pos){
		vecM[pos](paramindex) = sum;
		return;
	}

	if ((interval_end+1) == pos){
		valend = funcval(vecM, interval_end, vecM[pos].date, paramindex);
		if (valend != IOUtils::nodata)
			sum += valend;

		vecM[pos](paramindex) = sum;
		return;
	} else {
		for (size_t ii=interval_end+1; ii<pos; ii++){
			const double& val = vecM[ii](paramindex);
			if (val != IOUtils::nodata)
				sum += val;
		}
	}

	valend = funcval(vecM, pos, vecM[pos].date, paramindex);

	if (valend != IOUtils::nodata)
		sum += valend;

	vecM[pos](paramindex) = sum;
	//TODO:check if at least one point has been summed. If not -> nodata
}

double ResamplingAlgorithms::funcval(const std::vector<MeteoData>& vecM, const size_t& index,
                                     const Date& date, const size_t& paramindex)
{
	size_t start = index;
	if (vecM[start].isResampled()){
		if (start > 0){
			start--;
		} else {
			return IOUtils::nodata;
		}
	}

	size_t end   = index+1;
	const double& valstart = vecM[start](paramindex);

	if (!vecM[index].isResampled() && (vecM[index].date == date))
		return valstart;

	if ((vecM[end].isResampled())) //skip resampled value
		end++;

	//index either points to the element with date directly or to the last element with vecM[index] <= date
	if (end < vecM.size()){
		const double& valend   = vecM[end](paramindex);

		if (valend == IOUtils::nodata)
			return IOUtils::nodata;

		return linearInterpolation(vecM[start].date.getJulianDate(true), 0.0,
		                           vecM[end].date.getJulianDate(true), valend,
		                           date.getJulianDate(true));
	}

	return IOUtils::nodata;
}

/**
 * @brief This function solves the equation y = kx + d for two given points and returns y for a given x
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
		throw IOException("Attempted division by null", AT);

	//Solving y = kx +d
	double k = (y1 - y2) / (x1 - x2);
	double d = y2 - k*x2;

	return (k*x3 + d);
}

} //namespace

