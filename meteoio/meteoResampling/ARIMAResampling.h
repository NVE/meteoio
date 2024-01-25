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
#ifndef ARIMARESAMPLING_H
#define ARIMARESAMPLING_H

#include "ARIMAutils.h"    // change the includes to make it uniform
#include "InterpolARIMA.h" // change the includes to make it uniform
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <vector>

namespace mio {

    /**
     * @brief This class is designed to handle interpolation (resampling) of data using the ARIMA (AutoRegressive Integrated Moving Average) model
     * 
     * @details It uses the InterpolARIMA class, see InterpolARIMA, to perform interpolation using ARIMA.
     *
     * Gaps in the data are detected, and when possible data before and after the gap is used to interpolate the missing values. Otherwise, either 
     * only the data before or after the gap (depending on what is available) is used to predict the missing values. 
     * 
     * The ARIMA model needs constant sampling rates, therefore the most likely rate is calculated in the data used to fit, and the data is resampled 
     * to that sampling rate. If a requested point falls in between available data, it will be linearly interpolated.
     * 
     * Missing values in the data used to fit the ARIMA Model, will be linearly interpolated as well.
     * 
     * A gap is defined as a period of missing data, that has at least 2 data points with the most likely sampling rate. Only 1 missing data point is
     * linearly interpolated (should maybe just return instead).
     *
     * Mandatory parameters:
     * - `BEFORE_WINDOW` : The time before a gap that will be used to accumulate data to fit the ARIMA model. 
     * - `AFTER_WINDOW` : The time after a gap that will be used to accumulate data to fit the ARIMA model.
     * (BEFORE_WINDOW + AFTER_WINDOW < window_size)
     * 
     * Optional parameters:
     * - `MAX_P` : The maximum number of AR coefficients to use in the ARIMA model. Default: 8
     * - `MAX_D` : The maximum number of differences to use in the ARIMA model. Default: 3
     * - `MAX_Q` : The maximum number of MA coefficients to use in the ARIMA model. Default: 8
     * - `START_P` : The starting number of AR coefficients to use in the ARIMA model. Default: 2
     * - `START_Q` : The starting number of MA coefficients to use in the ARIMA model. Default: 2
     * - `MAX_P_SEASONAL` : The maximum number of seasonal AR coefficients to use in the ARIMA model. Default: 2
     * - `MAX_D_SEASONAL` : The maximum number of seasonal differences to use in the ARIMA model. Default: 1
     * - `MAX_Q_SEASONAL` : The maximum number of seasonal MA coefficients to use in the ARIMA model. Default: 2
     * - `START_P_SEASONAL` : The starting number of seasonal AR coefficients to use in the ARIMA model. Default: 1
     * - `START_Q_SEASONAL` : The starting number of seasonal MA coefficients to use in the ARIMA model. Default: 1
     * - `SEASONAL_PERIOD` : The period of the seasonal component. Default: 0 (no seasonal component)
     * - `LIK_METHOD` : The method used to fit the ARIMA model. Default: CSS-MLE
     * - `OPT_METHOD` : The optimization method used to fit the ARIMA model. Default: BFGS
     * - `STEPWISE` : Whether to use stepwise search of the best ARIMA model. Default: true (faster)
     * - `APPROXIMATION` : Whether to use approximation to determin the Information Criteria and the Likelihood. Default: true
     * - `NUM_MODELS` : The number of models to try when using stepwise search. Default: 94
     * - `SEASONAL` : Whether to use a seasonal component in the ARIMA model. Default: true
     * - `STATIONARY` : Whether to use a stationary ARIMA model. Default: false
     * 
     * 
     * @code
     * [Interpolations1D]
     * TA::resample = ARIMA
     * TA::ARIMA::BEFORE_WINDOW = 86400
     * TA::ARIMA::AFTER_WINDOW = 86400
     * 
     * @endcode
     * 
     * @note In the case that only random/random walk arima models are found, the missing values will not be filled (It would be just the mean otherwise)
     * 
     * @author Patrick Leibersperger
     * @date 2024-01-25
     * 
     * TODO: - do i need to use getJulian(true) or does it not matter?
     *       - how do i avoid the whole size_t vs int problem?
     *       - should i do the linear interpolation, or just return?
     *       - vecM is filled with the new value right? (avoid having nodata in before data)
     *       - should i set a limit on the prediction steps?
     *
     */
    class ARIMAResampling : public ResamplingAlgorithms {
    public:
        ARIMAResampling(const std::string &i_algoname, const std::string &i_parname, const double &dflt_window_size,
                        const std::vector<std::pair<std::string, std::string>> &vecArgs);

        // Performs ARIMA interpolation the first time it is called, as this is also the first time the data is available
        void resample(const std::string &stationHash, const size_t &index, const ResamplingPosition &position, const size_t &paramindex,
                      const std::vector<MeteoData> &vecM, MeteoData &md);
        std::string toString() const;

    private:
        // ARIMA related data
        std::vector<ARIMA_GAP> gap_data;
        std::vector<std::vector<double>> filled_data;
        std::vector<std::vector<Date>> all_dates;

        // Window parameters
        double before_window, after_window;

        // User defined Metadata
        int max_p = 8, max_d = 3, max_q = 8;
        int start_p = 2, start_q = 2;
        int max_P = 2, max_D = 1, max_Q = 2;
        int start_P = 1, start_Q = 1;
        int period = 0;
        std::string method = "CSS-MLE", opt_method = "BFGS";
        bool stepwise = true, approximation = true;
        int num_models = 94;
        bool seasonal = true, stationary = false;

        // Flags
        bool is_zero_possible = false;
        bool checked_vecM = false;
        bool gave_warning_end = false;
        bool gave_warning_start = false;
        bool gave_warning_interpol = false;
        std::vector<bool> is_valid_gap_data;
        std::vector<bool> warned_about_gap;

        // Private methods
        void setMetaData(InterpolARIMA &arima);
        double interpolVecAt(const std::vector<MeteoData> &vecM, const size_t &idx, const Date &date, const size_t &paramindex);
        double interpolVecAt(const std::vector<double> &data, const std::vector<Date> &dates, const size_t &pos, const Date &date);
        std::vector<double> predictData(std::vector<double>& data, const std::string& direction, size_t startIdx_interpol, int length_gap_interpol, int period, ResamplingAlgorithms::ResamplingPosition re_position);
    };
} // end namespace mio

#endif
