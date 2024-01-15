// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include "ARIMAResampling.h" // change include to make it uniform

#include <sstream>

namespace mio {

ARIMAResampling::ARIMAResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
             : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), gap_data(), data(), before_window_idx(0), after_window_idx(0) 
{
	const std::string where( "Interpolations1D::"+i_parname+"::"+i_algoname );
	if (!vecArgs.empty()) //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for \""+where+"\"", AT);
	
	before_window = after_window = 0.;
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="BEFORE_WINDOW") {
			IOUtils::parseArg(vecArgs[ii], where, before_window);
			before_window /= 86400.; //user uses seconds, internally julian day is used
		} else if (vecArgs[ii].first=="AFTER_WINDOW") {
			IOUtils::parseArg(vecArgs[ii], where, after_window);
			after_window /= 86400.; //user uses seconds, internally julian day is used
		} else if (vecArgs[ii].first=="MAX_P") {
			IOUtils::parseArg(vecArgs[ii], where, max_p);
		} else if (vecArgs[ii].first=="MAX_D") {
			IOUtils::parseArg(vecArgs[ii], where, max_d);
		} else if (vecArgs[ii].first=="MAX_Q") {
			IOUtils::parseArg(vecArgs[ii], where, max_q);
		} else if (vecArgs[ii].first=="MAX_P_SEASONAL") {
			IOUtils::parseArg(vecArgs[ii], where, max_P);
		} else if (vecArgs[ii].first=="MAX_D_SEASONAL") {
			IOUtils::parseArg(vecArgs[ii], where, max_D);
		} else if (vecArgs[ii].first=="MAX_Q_SEASONAL") {
			IOUtils::parseArg(vecArgs[ii], where, max_Q);
		} else if (vecArgs[ii].first=="SEASONAL_PERIOD") {
			IOUtils::parseArg(vecArgs[ii], where, s);
		} else if (vecArgs[ii].first=="LIK_METHOD") {
			IOUtils::parseArg(vecArgs[ii], where, method);
		} else if (vecArgs[ii].first=="OPTIMIZATION_METHOD") {
			IOUtils::parseArg(vecArgs[ii], where, opt_method);
		} else if (vecArgs[ii].first=="STEPWISE") {
			IOUtils::parseArg(vecArgs[ii], where, stepwise);
		} else if (vecArgs[ii].first=="APPROXIMATION") {
			IOUtils::parseArg(vecArgs[ii], where, approximation);
		} else if (vecArgs[ii].first=="NUM_MODELS") {
			IOUtils::parseArg(vecArgs[ii], where, num_models);
		} else if (vecArgs[ii].first=="SEASONAL") {
			IOUtils::parseArg(vecArgs[ii], where, seasonal);
		} else if (vecArgs[ii].first=="STATIONARY") {
			IOUtils::parseArg(vecArgs[ii], where, stationary);
		}
		else {
			throw InvalidArgumentException("Unknown argument \""+vecArgs[ii].first+"\" for \""+where+"\"", AT);
		}
	}

	if (!(before_window!=0 || after_window!=0)) throw InvalidArgumentException("Please provide a PERIOD for "+where, AT);
	
}

std::string ARIMAResampling::toString() const
{
	//this should help when debugging, so output relevant parameters for your algorithm
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
	return ss.str();
}

void ARIMAResampling::resample(const std::string& /*stationHash*/, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                            const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size()) throw IOException("The index of the element to be resampled is out of bounds", AT);

	// the first time the resample function is called, the data is processed
	if (gap_data.empty()) {
		
	}

	// either read the data, as it is within the same sampling rate
	// or interpolate the data, as it is not within the same sampling rate
	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		//implement here how value==nodata should be handled
	}

	return;
}

} //namespace
