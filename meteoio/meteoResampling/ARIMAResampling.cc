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
#include <unistd.h>

#include <sstream>

namespace mio {

    ARIMAResampling::ARIMAResampling(const std::string &i_algoname, const std::string &i_parname, const double &dflt_window_size,
                                     const std::vector<std::pair<std::string, std::string>> &vecArgs)
        : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs), gap_data(), filled_data(), all_dates(), before_window(),
          after_window() {
        const std::string where("Interpolations1D::" + i_parname + "::" + i_algoname);
        if (vecArgs.empty()) // incorrect arguments, throw an exception
            throw InvalidArgumentException("Wrong number of arguments for \"" + where + "\"", AT);

        before_window = after_window = 0.;
        for (size_t ii = 0; ii < vecArgs.size(); ii++) {
            if (vecArgs[ii].first == "BEFORE_WINDOW") {
                IOUtils::parseArg(vecArgs[ii], where, before_window);
                before_window /= 86400.; // user uses seconds, internally julian day is used
            } else if (vecArgs[ii].first == "AFTER_WINDOW") {
                IOUtils::parseArg(vecArgs[ii], where, after_window);
                after_window /= 86400.; // user uses seconds, internally julian day is used
            } else if (vecArgs[ii].first == "MAX_P") {
                IOUtils::parseArg(vecArgs[ii], where, max_p);
            } else if (vecArgs[ii].first == "MAX_D") {
                IOUtils::parseArg(vecArgs[ii], where, max_d);
            } else if (vecArgs[ii].first == "MAX_Q") {
                IOUtils::parseArg(vecArgs[ii], where, max_q);
            } else if (vecArgs[ii].first == "MAX_P_SEASONAL") {
                IOUtils::parseArg(vecArgs[ii], where, max_P);
            } else if (vecArgs[ii].first == "MAX_D_SEASONAL") {
                IOUtils::parseArg(vecArgs[ii], where, max_D);
            } else if (vecArgs[ii].first == "MAX_Q_SEASONAL") {
                IOUtils::parseArg(vecArgs[ii], where, max_Q);
            } else if (vecArgs[ii].first == "SEASONAL_PERIOD") {
                IOUtils::parseArg(vecArgs[ii], where, period);
                period /= 86400.;
            } else if (vecArgs[ii].first == "LIK_METHOD") {
                IOUtils::parseArg(vecArgs[ii], where, method);
            } else if (vecArgs[ii].first == "OPTIMIZATION_METHOD") {
                IOUtils::parseArg(vecArgs[ii], where, opt_method);
            } else if (vecArgs[ii].first == "STEPWISE") {
                IOUtils::parseArg(vecArgs[ii], where, stepwise);
            } else if (vecArgs[ii].first == "APPROXIMATION") {
                IOUtils::parseArg(vecArgs[ii], where, approximation);
            } else if (vecArgs[ii].first == "NUM_MODELS") {
                IOUtils::parseArg(vecArgs[ii], where, num_models);
            } else if (vecArgs[ii].first == "SEASONAL") {
                IOUtils::parseArg(vecArgs[ii], where, seasonal);
            } else if (vecArgs[ii].first == "STATIONARY") {
                IOUtils::parseArg(vecArgs[ii], where, stationary);
            } else {
                throw InvalidArgumentException("Unknown argument \"" + vecArgs[ii].first + "\" for \"" + where + "\"", AT);
            }
        }

        if (!(before_window != 0 || after_window != 0))
            throw InvalidArgumentException("Please provide a ARIMA window for " + where, AT);
        if (before_window + after_window > window_size)
            throw InvalidArgumentException("The ARIMA window is larger than the resampling window for " + where, AT);
    }

    std::string ARIMAResampling::toString() const {
        // this should help when debugging, so output relevant parameters for your algorithm
        std::ostringstream ss;
        ss << std::right << std::setw(10) << parname << "::" << std::left << std::setw(15) << algo << "[ ]" << std::endl;
        ss << "  Amount of found gaps " << gap_data.size() << std::endl;
        ss << "  Amount of filled data " << filled_data.size() << std::endl;
        ss << "  Amount of dates " << all_dates.size() << std::endl;
        ss << "  interpolated data: " << std::endl;
        ss << convertVectorsToString(filled_data);
        return ss.str();
    }

    double ARIMAResampling::interpolVecAt(const std::vector<MeteoData> &vecMet, const size_t &idx, const Date &date,
                                          const size_t &paramindex) {
        if (idx >= vecMet.size())
            throw IOException("The index of the element to be resampled is out of bounds", AT);
        else if (idx == vecMet.size() - 1) {
#ifdef DEBUG
            std::cout << "Extrapolation needed to gather ARIMA data, will use constant value";
#endif
            return vecMet[idx](paramindex);
        }
        MeteoData p1 = vecMet[idx];
        MeteoData p2 = vecMet[idx + 1];
        return linearInterpolation(p1.date.getJulian(true), p1(paramindex), p2.date.getJulian(true), p2(paramindex), date.getJulian(true));
    }

    double ARIMAResampling::interpolVecAt(const std::vector<double> &data, const std::vector<Date> &datesVec, const size_t &pos,
                                          const Date &date) {
        if (pos >= data.size())
            throw IOException("The index of the element to be resampled is out of bounds, for date: " + date.toString(Date::ISO), AT);
        else if (pos == data.size() - 1) {
#ifdef DEBUG
            std::cout << "Extrapolation needed to gather ARIMA data, will use constant value";
#endif
            return data[pos];
        }
        return linearInterpolation(datesVec[pos].getJulian(), data[pos], datesVec[pos + 1].getJulian(), data[pos + 1],
                                   date.getJulian(true));
    }

    std::vector<double> ARIMAResampling::fillGapWithPrediction(std::vector<double> &data, const std::string &direction,
                                                               const size_t &startIdx, const int &length, const int &period,
                                                               ResamplingAlgorithms::ResamplingPosition re_position) {
        InterpolARIMA arima(data, startIdx, length, direction, period);
        setMetaData(arima);
        std::vector<double> predictions = arima.predict();
        std::copy(predictions.begin(), predictions.end(), data.begin() + startIdx);
        return predictions;
    }

    std::vector<Date>::iterator findDate(std::vector<Date> &gap_dates, Date &resampling_date) {
        auto exactTime = [&resampling_date](Date curr_date) { return requal(curr_date, resampling_date); };
        return std::find_if(gap_dates.begin(), gap_dates.end(), exactTime);
    }

    void ARIMAResampling::setMetaData(InterpolARIMA &arima) {
        arima.setAutoArimaMetaData(max_p, max_d, max_q, start_p, start_q, max_P, max_D, max_Q, start_P, start_Q, seasonal, stationary);
        arima.setOptMetaData(method, opt_method, stepwise, approximation, num_models);
    }

    void ARIMAResampling::resample(const std::string & /*stationHash*/, const size_t &index, const ResamplingPosition &position,
                                   const size_t &paramindex, const std::vector<MeteoData> &vecM, MeteoData &md) {
        if (index >= vecM.size())
            throw IOException("The index of the element to be resampled is out of bounds", AT);

        if (position == ResamplingAlgorithms::exact_match) {
            const double value = vecM[index](paramindex);
            if (value != IOUtils::nodata) {
                md(paramindex) = value; // propagate value
                return;
            }
        }

		if (!checked_vecM) {
			// check if there are any zero values in the data, to see if a prediction of zero makes sense
			double sum = 0.0, sumOfSquares = 0.0, mean, standardDeviation = 0.0;

			for (size_t ii = 0; ii < vecM.size(); ++ii) {
				double value = vecM[ii](paramindex);
				sum += value;
				sumOfSquares += value * value;
			}

			mean = sum / vecM.size();
			standardDeviation = sqrt(sumOfSquares / vecM.size() - mean * mean);

			for (size_t ii = 0; ii < vecM.size(); ii++) {
				if (vecM[ii](paramindex) <= standardDeviation) {
					is_zero_possible = true;
					break;
				}
			}

			checked_vecM = true;
		}

        Date resampling_date = md.date;

        // check wether given position is in a known gap, if it is either return the
        // exact value or linearly interpolate, to get the correct value
        for (size_t ii = 0; ii < gap_data.size(); ii++) {
            ARIMA_GAP gap = gap_data[ii];
            std::vector<Date> gap_dates = all_dates[ii];
            std::vector<double> data_in_gap = filled_data[ii];
            bool is_valid_data = is_valid_gap_data[ii];
            if (resampling_date >= gap.startDate && resampling_date <= gap.endDate) {
                if (!is_valid_data) {
                    if (!warned_about_gap[ii]) {
                        std::cerr << "Could not find a useful ARIMA model, try other parameters, or another interpolation algorithm for "
                                     "data in between "
                                  << gap.startDate.toString(Date::ISO) << "---" << gap.endDate.toString(Date::ISO) << "||" << std::endl;
                        warned_about_gap[ii] = true;
                    }
                    return;
                }
                auto it = findDate(gap_dates, resampling_date);
                // if there is an exact match, return the data
                if (it != gap_dates.end()) {
                    size_t idx = std::distance(gap_dates.begin(), it);
                    double value = data_in_gap[idx];
                    if (value != IOUtils::nodata) {
                        md(paramindex) = value; // propagate value
                        return;
                    }
                    break;
                } else {
                    // otherwise linearly interpolate the data
                    size_t idx = std::distance(gap_dates.begin(), std::lower_bound(gap_dates.begin(), gap_dates.end(), resampling_date));
                    if (idx == gap_dates.size() - 1)
                        break;
                    else if ((data_in_gap[idx] == IOUtils::nodata || data_in_gap[idx + 1] == IOUtils::nodata))
                        break;
                    md(paramindex) = interpolVecAt(data_in_gap, gap_dates, idx, resampling_date);
                    return;
                }
            } else if (position == ResamplingAlgorithms::end && resampling_date > gap.endDate && resampling_date >= gap.startDate &&
                       gap.startDate == vecM[vecM.size() - 1].date) {
                if (!is_valid_data) {
                    if (!warned_about_gap[ii]) {
                        std::cerr << "Could not find a useful ARIMA model, try other parameters, or another interpolation algorithm for "
                                     "data in between "
                                  << gap.startDate.toString(Date::ISO) << "---" << gap.endDate.toString(Date::ISO) << "||" << std::endl;
                        warned_about_gap[ii] = true;
                    }
                    return;
                }
                if (!gave_warning_end) {
                    std::cerr << "Extrapolating more than 25 steps into the future is pointless, last known data point: "
                              << gap.startDate.toString(Date::ISO) << "||" << std::endl;
                    gave_warning_end = true;
                }
                return;
            }
        }

        // if it is not in a known gap, cache the gap, and interpolate it for subsequent calls
        size_t gap_start = IOUtils::npos;
        size_t gap_end = IOUtils::npos;
        ARIMA_GAP new_gap;
        Date data_start_date;
        Date data_end_date;
        if (position == ResamplingAlgorithms::end) {
            new_gap.startDate = vecM[vecM.size() - 1].date;
            new_gap.start = vecM.size() - 1;
            data_start_date = new_gap.startDate - window_size;
            new_gap.sampling_rate = computeSamplingRate(data_start_date, new_gap.startDate, vecM);
            new_gap.endDate = resampling_date + 25 / new_gap.sampling_rate;
            new_gap.end = vecM.size() - 1;
            data_end_date = new_gap.endDate;
            data_start_date = adjustStartDate(vecM, new_gap, data_start_date, data_end_date);
        } else {
            computeARIMAGap(new_gap, index, paramindex, vecM, resampling_date, gap_start, gap_end, before_window, after_window, window_size,
                            data_start_date, data_end_date);
        }

        // check if gap ended up being bigger than the window size
        if (new_gap.endDate - new_gap.startDate > window_size) {
            double difference = (new_gap.endDate - new_gap.startDate).getJulian(true) - window_size;
            difference *= 86400;
            throw IOException(
                "The window size is smaller than the data gap to be interpolated, please increase window size, by at least: " +
                    std::to_string(difference) + "s",
                AT);
        }

        if (new_gap.isGap()) {
            // data vector is of length (data_end_date - data_start_date) * sampling_rate
            int length = static_cast<int>((data_end_date - data_start_date).getJulian(true) * new_gap.sampling_rate) +
                         1; // otherwise end date is not included
            std::vector<double> data(length);
            std::vector<Date> dates(length);
            // get a vector of meteodata that contains the data before and after the gap
            std::vector<MeteoData> data_vec_before;
            std::vector<MeteoData> data_vec_after;
            // get the data before and after the gap
            for (size_t ii = 0; ii < vecM.size(); ii++) {
                Date date_v = vecM[ii].date;
                if (date_v >= data_start_date && date_v <= new_gap.startDate) {
                    data_vec_before.push_back(vecM[ii]);
                }
                if (date_v >= new_gap.endDate && date_v <= data_end_date) {
                    data_vec_after.push_back(vecM[ii]);
                }
            }

            bool has_data_before = data_vec_before.size() > 1;
            bool has_data_after = data_vec_after.size() > 1;

            if (data_vec_before.size() < 10 && data_vec_after.size() < 10 && !gave_warning_interpol) {
                std::cerr << "Not enough data to interpolate the gap" << std::endl;
                std::cerr << new_gap.toString() << std::endl;
                std::cerr << "Datapoints before the gap: " << data_vec_before.size() << std::endl;
                std::cerr << "Datapoints after the gap: " << data_vec_after.size() << "||" << std::endl;
                gave_warning_interpol = true;
                return;
            }

            // resample to the desired sampling rate
            int length_gap_interpol = 0;
            int endIdx_interpol = IOUtils::npos;
            size_t startIdx_interpol = (data_start_date == new_gap.startDate) ? 0 : IOUtils::npos;
            for (int i = 0; i < length; i++) {
                Date date = data_start_date + i / new_gap.sampling_rate;
                dates[i] = date;
                if (date >= data_start_date && date <= new_gap.startDate) {
                    if (has_data_before && i < data_vec_before.size()) {
                        if (requal(date, data_vec_before[i].date)) {
                            data[i] = data_vec_before[i](paramindex);
                        } else {
                            // linearly interpolate the data
                            data[i] = interpolVecAt(data_vec_before, i, date, paramindex);
                        }
                    } else {
                        // the data is missing, so set it to IOUtils::nodata
                        data[i] = IOUtils::nodata;
                    }
                } else if (date > new_gap.startDate && date < new_gap.endDate) {
                    if (startIdx_interpol == IOUtils::npos) { // TODO: this is not reached when the sampling rate is unfitting
                        startIdx_interpol = i;
                    }
                    // the data is missing, so set it to IOUtils::nodata
                    data[i] = IOUtils::nodata;
                } else if (date >= new_gap.endDate && date <= data_end_date) {
                    if (length_gap_interpol == 0) {
                        length_gap_interpol = i - static_cast<int>(startIdx_interpol);
                        endIdx_interpol = i;
                    }
                    int after_id = i - endIdx_interpol;
                    if (has_data_after && after_id < data_vec_after.size()) {
                        if (requal(date, data_vec_after[after_id].date)) {
                            data[i] = data_vec_after[after_id](paramindex);
                        } else {
                            // linearly interpolate the data
                            data[i] = interpolVecAt(data_vec_after, after_id, date, paramindex);
                        }
                    } else {
                        data[i] = IOUtils::nodata;
                    }
                }
            }

            if (data_end_date == new_gap.endDate)
                length_gap_interpol = data.size() - startIdx_interpol;

            // now fill the data with the arima model
            // either by interpolating or predicting forward or backward
            std::vector<double> interpolated_data;
            if (data_vec_before.size() < 8 && data_vec_after.size() > 8) {
#ifdef DEBUG
                std::cout << "predicting backward" << std::endl;
#endif
                interpolated_data = fillGapWithPrediction(data, "backward", startIdx_interpol, length_gap_interpol, period, position);
            } else if (data_vec_after.size() < 8 && data_vec_before.size() > 8) {
#ifdef DEBUG
                std::cout << "predicting forward" << std::endl;
#endif
                interpolated_data = fillGapWithPrediction(data, "forward", startIdx_interpol, length_gap_interpol, period, position);
            } else if (data_vec_before.size() < 8 && data_vec_after.size() < 8) {
                throw IOException("Could not accumulate enough data for parameter estimation; Increasing window sizes might help");
            } else {
                InterpolARIMA arima(data, startIdx_interpol, length_gap_interpol, period);
                setMetaData(arima);
                arima.interpolate();
                interpolated_data = arima.getInterpolatedData();
            }

            std::vector<Date> interpolated_dates(dates.begin() + startIdx_interpol,
                                                 dates.begin() + startIdx_interpol + length_gap_interpol);
            // need to push the first value after the gap to the interpolated data and dates, otherwise the interpolation will fail in case
            // of requested value between data points
            if (data_end_date != new_gap.endDate) {
                interpolated_data.push_back(data[endIdx_interpol]);
                interpolated_dates.push_back(dates[endIdx_interpol]);
            }

            // check whether the interpolated data has only zeros, or if it has zeros when the original data does not
            bool only_zeros = true;
            bool any_zeros = false;
            for (size_t ii = 0; ii < interpolated_data.size(); ii++) {
                if (interpolated_data[ii] != 0.0)
                    only_zeros = false;
                if (interpolated_data[ii] == 0.0)
                    any_zeros = true;
            }

            gap_data.push_back(new_gap);
            if (only_zeros) {
                is_valid_gap_data.push_back(false);
            } else if (any_zeros && !is_zero_possible) {
                is_valid_gap_data.push_back(false);
            } else {
                is_valid_gap_data.push_back(true);
            }
            warned_about_gap.push_back(false);
            filled_data.push_back(interpolated_data);
            all_dates.push_back(interpolated_dates);

            if (!is_valid_gap_data[gap_data.size() - 1]) {
                if (!warned_about_gap[warned_about_gap.size() - 1]) {
                    std::cerr << "Could not find a useful ARIMA model, try other parameters, or another interpolation algorithm for data "
                                 "in between "
                              << new_gap.startDate.toString(Date::ISO) << "---" << new_gap.endDate.toString(Date::ISO) << "||" << std::endl;
                    warned_about_gap[warned_about_gap.size() - 1] = true;
                }
                return;
            }

            auto it = findDate(interpolated_dates, resampling_date);
            // if there is an exact match, return the data
            if (it != interpolated_dates.end()) {
                size_t idx = std::distance(interpolated_dates.begin(), it);
                md(paramindex) = interpolated_data[idx];
                return;
            } else {
                // otherwise linearly interpolate the data
                size_t idx = std::distance(interpolated_dates.begin(),
                                           std::lower_bound(interpolated_dates.begin(), interpolated_dates.end(), resampling_date));
                md(paramindex) = interpolVecAt(interpolated_data, interpolated_dates, idx, resampling_date);
                return;
            }

        } else {
            // linearly interpolate the point
#ifdef DEBUG
            std::cout << "linearly interpolating the point" << std::endl;
#endif
            double start_value = vecM[gap_start - 1](paramindex);
            double end_value = vecM[gap_end + 1](paramindex);
            Date start_date = vecM[gap_start - 1].date;
            Date end_date = vecM[gap_end + 1].date;
            md(paramindex) = linearInterpolation(start_date.getJulian(true), start_value, end_date.getJulian(true), end_value,
                                                 resampling_date.getJulian(true));
            return;
        }
        std::cerr << "something went seriously wrong in the arima interpolations, please contact the developers!" << std::endl;
        return;
    }
} // namespace
