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
             : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), gap_data(), filled_data(), sampling_rate(0)
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

	if (!(before_window!=0 || after_window!=0)) throw InvalidArgumentException("Please provide a ARIMA window for "+where, AT);
	
}

std::string ARIMAResampling::toString() const
{
	//this should help when debugging, so output relevant parameters for your algorithm
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
	return ss.str();
}


double ARIMAResampling::interpolVecAt(const std::vector<MeteoData> &vecMet,const size_t &idx, const Date &date, const size_t &paramindex) {
	MeteoData p1 = vecMet[idx];
	MeteoData p2 = vecMet[idx+1];
    return linearInterpolation(p1.date.getJulian(true), p1(paramindex), p2.date.getJulian(true), p2(paramindex), date.getJulian(true));
}

double ARIMAResampling::interpolVecAt(const std::vector<double> &data, const std::vector<Date> &datesVec, const size_t &pos, const Date &date) {
	return linearInterpolation(datesVec[pos].getJulian(), data[pos], datesVec[pos+1].getJulian(), data[pos+1], date.getJulian(true));
}

void fillGapWithPrediction(std::vector<double>& data, const std::string& direction, const size_t &startIdx, const size_t &length, const int &period) {
	InterpolARIMA arima(data, length, direction, period);
    std::vector<double> predictions = arima.predict();
    std::copy(predictions.begin(), predictions.end(), data.begin() + startIdx);
}

void ARIMAResampling::resample(const std::string& /*stationHash*/, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                            const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (index >= vecM.size()) throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		if (value != IOUtils::nodata) {
			md(paramindex) = value; //propagate value
			return;
		}
	}

	Date resampling_date = md.date;
	// check wether given position is in a known gap, if it is either return the 
	// exact value or linearly interpolate, to get the correct value
	for (int ii = 0; ii < gap_data.size(); ii++) {
		ARIMA_GAP gap = gap_data[ii];
		std::vector<Date> gap_dates = all_dates[ii];
		std::vector<double> gap_data = filled_data[ii];
		if (resampling_date >= gap.startDate && resampling_date <= gap.endDate) {
			auto exactTime = [&resampling_date](Date curr_date) {return requal(curr_date, resampling_date);};
			auto it = std::find_if(gap_dates.begin(), gap_dates.end(), exactTime);

			// if there is an exact match, return the data
			if (it != gap_dates.end()) {
				size_t idx = std::distance(gap_dates.begin(), it);
				md(paramindex) = gap_data[idx];
				return;
			} else {
				// otherwise linearly interpolate the data
				size_t idx = std::distance(gap_dates.begin(), std::lower_bound(gap_dates.begin(), gap_dates.end(), resampling_date));
				md(paramindex) = interpolVecAt(gap_data, gap_dates, idx, resampling_date);
				return;
			}
		}
	}

	// if it is not in a known gap, cache the gap, and interpolate it for subsequent calls
	size_t gap_start = IOUtils::npos;
	size_t gap_end = IOUtils::npos;
	ARIMA_GAP new_gap;
	computeGap(new_gap,index, paramindex, vecM, resampling_date, window_size, gap_start, gap_end);

	if (new_gap.isGap()) {
		Date data_start_date = gap_start - before_window;
		Date data_end_date = gap_end + after_window;
		new_gap.sampling_rate = computeSamplingRate(data_start_date, data_end_date, vecM);
		

		// data vector is of length (data_end_date - data_start_date) * sampling_rate
		int length = (data_end_date - data_start_date).getJulian(true) * sampling_rate;
		std::vector<double> data(length);
		std::vector<Date> dates(length);
		// get a vector of meteodata that contains the data before and after the gap
		std::vector<MeteoData> data_vec_before;
		std::vector<MeteoData> data_vec_after;
		// get the data before the gap
		for (int ii=0; ii<vecM.size(); ii++) {
			if (vecM[ii].date >= data_start_date && ii < gap_start) {
				data_vec_before.push_back(vecM[ii]);
			}
			if (ii >= gap_end && vecM[ii].date <= data_end_date) {
				data_vec_after.push_back(vecM[ii]);
			}
		}

		if (data_vec_before.size() < 10 && data_vec_after.size() < 10) {
			std::cerr << "Not enough data to interpolate the gap" << std::endl;
			std::cerr << "Datapoints before the gap: " << data_vec_before.size() << std::endl;
			std::cerr << "Datapoints after the gap: " << data_vec_after.size() << std::endl;
			return;
		}

		// resample to the desired sampling rate
		size_t length_gap_interpol = IOUtils::npos;
		size_t startIdx_interpol = IOUtils::npos;
		for (int i =0; i < length; i++) {
			Date date = data_start_date + i * sampling_rate;
			dates[i] = date;
			if (date >=data_start_date && date < new_gap.startDate) {
				if (requal(date, data_vec_before[i].date)) {
					data[i] = data_vec_before[i](paramindex);
				} else {
					// linearly interpolate the data
					data[i] = interpolVecAt(data_vec_before, i, date, paramindex);
				}
			} else if (date >= new_gap.startDate && date <= new_gap.endDate) {
				if (startIdx_interpol == IOUtils::npos) {
					startIdx_interpol = i;
				}
				// the data is missing, so set it to IOUtils::nodata
				data[i] = IOUtils::nodata;
			} else if (date > new_gap.endDate && date <= data_end_date) {
				if (length_gap_interpol == IOUtils::npos) {
					length_gap_interpol = i - 	startIdx_interpol-1;
				}
				if (requal(date, data_vec_after[i].date)) {
					data[i] = data_vec_after[i](paramindex);
				} else {
					// linearly interpolate the data
					data[i] = interpolVecAt(data_vec_after, i, date, paramindex);
				}
			}
		}

		// now fill the data with the arima model
		// either by interpolating or predicting forward or backward
		if (data_vec_before.size()<10) {
			fillGapWithPrediction(data, "backward", startIdx_interpol, length_gap_interpol, s);
		} else if (data_vec_after.size()<10) {
			fillGapWithPrediction(data, "forward", startIdx_interpol, length_gap_interpol, s);
		} else {
			InterpolARIMA arima(data, startIdx_interpol, length_gap_interpol, s);
			arima.fillGap();
		}
		
		gap_data.push_back(new_gap);
		std::vector<double> interpolated_data(data.begin() + startIdx_interpol, data.begin() + startIdx_interpol + length_gap_interpol);
		std::vector<MeteoData> interpolated_dates(dates.begin() + startIdx_interpol, dates.begin() + startIdx_interpol + length_gap_interpol);
		filled_data.push_back(interpolated_data);
		all_dates.push_back(dates);
	} else {
		// linearly interpolate the point
		double start_value = vecM[gap_start-1](paramindex);
		double end_value = vecM[gap_end+1](paramindex);
		Date start_date = vecM[gap_start-1].date;
		Date end_date = vecM[gap_end+1].date;
		md(paramindex) = linearInterpolation(start_date.getJulian(true), start_value, end_date.getJulian(true), end_value, resampling_date.getJulian(true));
	}
	return;
}
} //namespace
