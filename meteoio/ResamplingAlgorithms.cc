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
 * @section algorithms_available Available Resampling Algorithms
 * Two algorithms for the resampling are implemented:
 * - linear: linear data resampling, see ResamplingAlgorithms::LinearResampling
 * - nearest_neighbour:  data resampling, see ResamplingAlgorithms::NearestNeighbour
 */

std::map<std::string, resamplingptr> ResamplingAlgorithms::algorithmMap;
const bool ResamplingAlgorithms::__init = ResamplingAlgorithms::initStaticData();

bool ResamplingAlgorithms::initStaticData()
{
	algorithmMap["linear"]            = &ResamplingAlgorithms::LinearResampling;
	algorithmMap["nearest_neighbour"] = &ResamplingAlgorithms::NearestNeighbour;
	algorithmMap["accumulate"] = &ResamplingAlgorithms::Accumulate;

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

void ResamplingAlgorithms::NearestNeighbour(const unsigned int& pos, const MeteoData::Parameters& paramindex,
                                            const std::vector<std::string>& /*taskargs*/,
                                            std::vector<MeteoData>& vecM, std::vector<StationData>& /*vecS*/)
{
	if (pos >= vecM.size())
		throw IOException("The position of the resampled element is out of bounds", AT);


	if (vecM[pos].param(paramindex) != IOUtils::nodata) //if there is a value, don't resample
		return;

	//Try to find the nearest neighbour, if there are two equally distant, then return the arithmetic mean
	MeteoData m1, m2;
	bool found1=false, found2=false;
	for (unsigned int ii=pos+1; ii<vecM.size(); ii++){
		if (vecM[ii].param(paramindex) != IOUtils::nodata){
			m1 = vecM[ii];
			found1 = true;
			break;
		}
	}

	for (unsigned int ii=0; ii<(unsigned int)pos; ii++){
		if (vecM[ii].param(paramindex) != IOUtils::nodata){
			m2 = vecM[ii];
			found2 = true;
		}
	}

	if (found1 && !found2){ //nearest neighbour on found after index 'pos'
		vecM[pos].param(paramindex) = m1.param(paramindex);
	} else if (!found1 && found2){ //nearest neighbour on found before index 'pos'
		vecM[pos].param(paramindex) = m2.param(paramindex);
	} else if (!found1 && !found2){ // no nearest neighbour with a value different from IOUtils::nodata
		vecM[pos].param(paramindex) = IOUtils::nodata;
	} else {
		Date diff1 = m1.date - vecM[pos].date; //calculate time interval to element at pos
		Date diff2 = vecM[pos].date - m2.date; //calculate time interval to element at pos

		if (IOUtils::checkEpsilonEquality(diff1.getJulianDate(), diff2.getJulianDate(), 0.1/1440)){ //within 6 seconds
			vecM[pos].param(paramindex) = Interpol1D::linearInterpolation(m1.param(paramindex), m2.param(paramindex), 0.5);
		} else if (diff1 < diff2){
			vecM[pos].param(paramindex) = m1.param(paramindex);
		} else if (diff1 > diff2){
			vecM[pos].param(paramindex) = m2.param(paramindex);
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
void ResamplingAlgorithms::LinearResampling(const unsigned int& pos, const MeteoData::Parameters& paramindex,
                                            const std::vector<std::string>& taskargs,
                                            std::vector<MeteoData>& vecM, std::vector<StationData>& /*vecS*/)
{
	if (pos >= vecM.size())
		throw IOException("The position of the resampled element is out of bounds", AT);

	if (vecM[pos].param(paramindex) != IOUtils::nodata) //if there is a value, don't resample
		return;

	//Check whether extrapolation is activated
	bool extrapolate = false;
	if ((taskargs.size()==1) && (taskargs[0]=="extrapolate"))
		extrapolate = true;

	//Now find two points within the vecM (before and aft, that are not IOUtils::nodata)
	//If that condition cannot be met, simply add nodata for the resampled value (exception: extrapolate)
	unsigned int indexP1=IOUtils::npos, indexP2=IOUtils::npos;
	bool foundP1=false, foundP2=false;

	for (unsigned int ii=pos; (ii--) > 0; ){
		if (vecM[ii].param(paramindex) != IOUtils::nodata){
			indexP1=ii;
			foundP1 = true;
			break;
		}
	}		

	for (unsigned int ii=pos+1; ii<vecM.size(); ii++){
		if (vecM[ii].param(paramindex) != IOUtils::nodata){
			indexP2 = ii;
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
		for (unsigned int ii=indexP2+1; ii<vecM.size(); ii++){
			if (vecM[ii].param(paramindex) != IOUtils::nodata){
				indexP1 = ii;
				foundP1 = true;
				break;
			}
		}		
	} else if (foundP1 && !foundP2){ //only nodata found after pos, try looking before indexP1
		for (unsigned int ii=indexP1; (ii--) > 0; ){
			if (vecM[ii].param(paramindex) != IOUtils::nodata){
				indexP2=ii;
				foundP2 = true;
				break;
			}
		}		
	}

	if (!foundP1 || !foundP2) //now at least two points need to be present
		return;

	//At this point indexP1 and indexP2 point to values that are different from IOUtils::nodata
	const double& val1 = vecM[indexP1].param(paramindex);
	const double& val2 = vecM[indexP2].param(paramindex);
	vecM[pos].param(paramindex) = Interpol1D::linearInterpolation(vecM[indexP1].date.getJulianDate(), val1,
	                              vecM[indexP2].date.getJulianDate(), val2,
	                              vecM[pos].date.getJulianDate());
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
void ResamplingAlgorithms::Accumulate(const unsigned int& pos, const MeteoData::Parameters& paramindex,
                                      const std::vector<std::string>& taskargs,
                                      std::vector<MeteoData>& vecM, std::vector<StationData>& vecS)
{
	//Not ready for usage yet...
	throw IOException("Not finished...", AT);

	unsigned int npos = pos;
	if (pos >= vecM.size())
		throw IOException("The position of the resampled element is out of bounds", AT);

	//Get accumulation period
	double accumulate_period;
	if (taskargs.size()==1) {
		IOUtils::convertString(accumulate_period, taskargs[0]);
		if(accumulate_period<=0.) {
			std::stringstream tmp;
			tmp << "Invalid accumulation period (" << accumulate_period << ") ";
			tmp << "for parameter " << MeteoData::getParameterName(paramindex);
			throw InvalidArgumentException(tmp.str(), AT);
		}
	} else {
		std::stringstream tmp;
		tmp << "Please provide accumulation period (in seconds) for parameter " << MeteoData::getParameterName(paramindex);
		throw InvalidArgumentException(tmp.str(), AT);
	}
	
	//find start of accumulation period
	bool found_start=false;
	int ii=pos;
	const Date Start(vecM[pos].date.getJulianDate() - accumulate_period/(24.*3600.));
	for (; ii>=0; ii--) {
		if(vecM[(unsigned int)ii].date<=Start) {
			found_start=true;
			break;
		}
	}

	if(found_start==false) {
		std::cerr << "[W] Could not accumulate " << MeteoData::getParameterName(paramindex);
		std::cerr << ", not enough data for accumulation period\n";
		return;
	}

	//resample the starting point
	std::vector<std::string> dummy_args;
	const unsigned int org_size = vecM.size();
	LinearResampling(ii, paramindex, dummy_args, vecM, vecS);
	if(vecM.size()>org_size)
		npos++;

	//start accumulating
	double sum = 0.0;
	bool exist=false;
	while(vecM[ii].date<=vecM[npos].date) {
		//we will at least enter once into this loop since accumulate_period>0
		const double& val = vecM[npos].param(paramindex);
		if (val != IOUtils::nodata) {
			sum += val;
			exist=true;
		}
		ii++;
	}

	//add last remaining point after interpolating it
	//HACK: TODO


	//check if at least one point has been summed. If not -> nodata
}


} //namespace

