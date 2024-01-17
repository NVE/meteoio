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

#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include "InterpolARIMA.h" // change the includes to make it uniform
#include "ARIMAutils.h" // change the includes to make it uniform
#include <vector>


namespace mio {

/**
 * @brief Brief description
 * @details
 * Longer description of the algorithm as well as example of use
 * 
 * Assuming a constant sampling rate a data vector will be created with the available data, 
 * all gaps will be located, and the gaps will be filled with the ARIMA model.
 * Points that fall in between interpolated data will be interpolated linearly.
 * a gap is defined as at least 3 missing values in a row
 * 
 * Also assuming that the meteovector fills all necessary fields with IOUtils::npos beforehand
 * so i dont need to work with dates
 * 
 * missing values in the data used to interpolate will be linearly interpolated
 * 
 * TODO: - do i need to use getJulian(true) or does it not matter?
 *       - predictions should only happen if the gap is at the end or the beginning
 *       - if interpolation is not possible because the window size is smaller than the gap, just predict one step ahead
 *       - be careful to not interpolate, if the gap still has missing data
 * 
 */
class ARIMAResampling : public ResamplingAlgorithms {
	public:
		ARIMAResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs);

        // it will do the ARIMA interpolation the first time it is called, as this is also the first time the data is available
		void resample(const std::string& stationHash, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
    
    private:

        std::vector<ARIMA_GAP> gap_data;
        std::vector<std::vector<double>> filled_data;
        std::vector<std::vector<Date>> all_dates;

        // depending of before and after windows are not available can predict in both directions
        double before_window, after_window;

        // User defined Metadata
        int max_p = 8, max_d = 3, max_q = 8;
        int start_p = 2, start_q = 2;
        int max_P = 2, max_D = 1, max_Q = 2;
        int start_P = 1, start_Q = 1;
        int period = 0;
        std::string method = "CSS-MLE", opt_method = "BFGS";
        bool stepwise = true, approximation = false;
        int num_models = 94;
        bool seasonal = true, stationary = false;  
        double interpolVecAt(const std::vector<MeteoData> &vecM, const size_t &idx, const Date &date, const size_t &paramindex);
        double interpolVecAt(const std::vector<double> &data, const std::vector<Date> &dates, const size_t &pos, const Date &date);
};

} //end namespace mio

#endif
