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

#ifndef INTERPOLARIMA_H
#define INTERPOLARIMA_H

#include <map>
#include <meteoio/thirdParty/ctsa.h>
#include <string>
#include <vector>
#include <iostream>


namespace mio {
/**
 * @class InterpolARIMA
 * 
 * @brief This class is used for interpolating or predicting missing data in a time series using the Auto ARIMA algorithm.
 * 
 * @details Depending on the constructor that is used, the data and auto ARIMA models are set up for either interpolation or prediction. 
 * The interpolation is interpolate(). The prediction methods is predict(). 
 * 
 * Interpolate will fill a gap in the data, whose start is specified by gap_loc 
 * and whose length is specified by N_gap. Data is assumed to be of equal time steps, and is split into two parts, data_before and data_after. 
 * So in the end data should be of size data_before + data_after + N_gap. The interpolation is done by fitting one ARIMA model to data_before and 
 * one to data_after. The ARIMA models are fitted using the auto.arima algorithm from the <a href="https://github.com/rafat/ctsa/tree/master">ctsa</a> (BSD-3 Clause, see below). The ARIMA models are then used to predict 
 * the missing data forward and backward in time. The final prediction is a weighted average of the two, where the weighting is done so more information 
 * comes from the closer data. 
 * 
 * Predict will predict the next n_steps values in the time series. It can either be forward in time (direction = "forward") or backward in time (direction = "backward").
 * For forward prediction data[0:gap_loc] is used to fit the ARIMA model, and for backward prediction data[gap_loc + N_gap:] is used.
 * 
 * For more Information concerning ARIMA see, https://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average, and https://otexts.com/fpp2/arima.html, and the final 
 * interpolation algorithm: https://www.tandfonline.com/doi/abs/10.1080/02664769624332?casa_token=fEVPFRYrr7sAAAAA:ozZFAcUWX4mKaUI8tvOn6R-3giOHefH0p8vaRDFCN1ORGy0d9evP7Hn9aLbMWsUQsIKrKEKxP-M
 * 
 * 
 * @note Interpolate is meant to only be used, when there is actually backward data available. If there is no backward data, then predict should be used instead. 
 *          Where predict is meant to be used in conjunction with the according constructor.
 *      Currently prediction forward or backward, when providing both is not implemented, but can be easily added on demand.
 * 
 * @author Patrick Leibersperger
 * @date 2024-01-25
 * 
 * 
 * Copyright (c) 2014, Rafat Hussain
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1.  Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
 	in the documentation and/or other materials provided with the distribution.

 3.  The name of the author may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 OF SUCH DAMAGE.
 */
    class InterpolARIMA {
    public:
        InterpolARIMA();
        InterpolARIMA(std::vector<double> data_in, size_t gap_loc, int N_gap, int s = 0);
        InterpolARIMA(std::vector<double> data_in, size_t gap_loc, int N_gap, std::vector<double> xreg_vec, int s = 0);
        InterpolARIMA(std::vector<double> data_in, size_t gap_loc, int n_predictions, std::string direction = "forward", int s = 0);

        // Setters
        void setAutoArimaMetaData(int max_p_param = 8, int max_d_param = 3, int max_q = 8, int start_p = 2, int start_q = 2, int max_P = 2,
                                  int max_D = 1, int max_Q = 2, int start_P = 1, int start_Q = 1, bool seasonal = true,
                                  bool stationary = false);
        void setOptMetaData(std::string method = "CSS-MLE", std::string opt_method = "BFGS", bool stepwise = true,
                            bool approximation = false, int num_models = 94);

        // Interpolation methods
        std::vector<double> simulate(int n_steps, int seed = 0);
        void fillGap();
        void interpolate();
        std::vector<double> predict(int n_steps = 0);
        std::vector<double> ARIMApredict(int n_steps);

        // Getters
        std::vector<double> getData() { return data; };
        std::vector<double> getForwardData() { return data_forward; };
        std::vector<double> getBackwardData() { return data_backward; };
        std::vector<double> getInterpolatedData();


        // Copy constructor
        InterpolARIMA(const InterpolARIMA &other)
            : gap_loc(other.gap_loc), N_gap(other.N_gap), time(other.time), pred_forward(other.pred_forward),
              pred_backward(other.pred_backward), data(other.data), xreg_vec_f(other.xreg_vec_f), xreg_vec_b(other.xreg_vec_b), 
                data_forward(other.data_forward), data_backward(other.data_backward), new_xreg_vec_f(other.new_xreg_vec_f), 
                new_xreg_vec_b(other.new_xreg_vec_b), N_data_forward(other.N_data_forward), N_data_backward(other.N_data_backward),
                max_p(other.max_p), max_d(other.max_d), max_q(other.max_q), start_p(other.start_p), start_q(other.start_q),
                max_P(other.max_P), max_D(other.max_D), max_Q(other.max_Q), start_P(other.start_P), start_Q(other.start_Q),
                r(other.r), s(other.s), method(other.method), opt_method(other.opt_method), stepwise(other.stepwise),
                approximation(other.approximation), num_models(other.num_models), seasonal(other.seasonal), stationary(other.stationary) 
        {    
            std::cout << "copy constructor called" << std::endl;
            auto_arima_forward = auto_arima_copy(other.auto_arima_forward);
            auto_arima_backward = auto_arima_copy(other.auto_arima_backward);
            xreg_f = (xreg_vec_f.empty()) ? NULL : &xreg_vec_f[0]; 
            xreg_b = (xreg_vec_b.empty()) ? NULL : &xreg_vec_b[0];
            new_xreg_f = (new_xreg_vec_f.empty()) ? NULL : &new_xreg_vec_f[0];
            new_xreg_b = (new_xreg_vec_b.empty()) ? NULL : &new_xreg_vec_b[0];
        }

        // Copy assignment operator
        InterpolARIMA& operator=(const InterpolARIMA &other)
        {
            std::cout << "copy assignment operator called" << std::endl;
            if (this != &other) // protect against invalid self-assignment
            {
                auto_arima_forward = auto_arima_copy(other.auto_arima_forward);
                auto_arima_backward = auto_arima_copy(other.auto_arima_backward);

                // 3: copy all the other fields from the other object
                gap_loc = other.gap_loc;
                N_gap = other.N_gap;
                time = other.time;
                pred_forward = other.pred_forward;
                pred_backward = other.pred_backward;
                data = other.data;
                xreg_vec_f = other.xreg_vec_f;
                xreg_vec_b = other.xreg_vec_b;
                data_forward = other.data_forward;
                data_backward = other.data_backward;
                new_xreg_vec_f = other.new_xreg_vec_f;
                new_xreg_vec_b = other.new_xreg_vec_b;
                N_data_forward = other.N_data_forward;
                N_data_backward = other.N_data_backward;
                max_p = other.max_p;
                max_d = other.max_d;
                max_q = other.max_q;
                start_p = other.start_p;
                start_q = other.start_q;
                max_P = other.max_P;
                max_D = other.max_D;
                max_Q = other.max_Q;
                start_P = other.start_P;
                start_Q = other.start_Q;
                r = other.r;
                s = other.s;
                method = other.method;
                opt_method = other.opt_method;
                stepwise = other.stepwise;
                approximation = other.approximation;
                num_models = other.num_models;
                seasonal = other.seasonal;
                stationary = other.stationary;

                // 4: handle the pointers to the vectors
                xreg_f = (xreg_vec_f.empty()) ? NULL : &xreg_vec_f[0]; 
                xreg_b = (xreg_vec_b.empty()) ? NULL : &xreg_vec_b[0];
                new_xreg_f = (new_xreg_vec_f.empty()) ? NULL : &new_xreg_vec_f[0];
                new_xreg_b = (new_xreg_vec_b.empty()) ? NULL : &new_xreg_vec_b[0];
            }
            // by convention, always return *this
            return *this;
        }

        ~InterpolARIMA() {
            auto_arima_free(auto_arima_forward);
            auto_arima_free(auto_arima_backward);
            delete xreg_f;
            delete xreg_b;
            delete new_xreg_f;
            delete new_xreg_b;
        }
        std::string toString();

    private:
        // Interpolation variables
        std::vector<double> data;
        size_t gap_loc;
        int N_gap;
        std::vector<double> time;
        std::vector<double> pred_forward, pred_backward;

        // Auto Arima variables
        // const doesnt work wiht c
        std::vector<double> xreg_vec_f, xreg_vec_b, data_forward, data_backward, new_xreg_vec_f, new_xreg_vec_b;
        double *xreg_f;
        double *xreg_b;
        double *new_xreg_f;
        double *new_xreg_b;
        std::vector<double> amse_forward, amse_backward;
        size_t N_data_forward, N_data_backward;
        int max_p = 8, max_d = 3, max_q = 8;
        int start_p = 2, start_q = 2;
        int max_P = 2, max_D = 1, max_Q = 2;
        int start_P = 1, start_Q = 1;
        int r = 0, s = 0;
        std::string method = "CSS-MLE", opt_method = "BFGS";
        bool stepwise = true, approximation = true;
        int num_models = 94;
        bool seasonal = true, stationary = false;

        std::map<std::string, int> method_map = {{"CSS-MLE", 0}, {"MLE", 1}, {"CSS", 2}};
        std::map<std::string, int> opt_method_map = {{"Nelder-Mead", 0},
                                                     {"Newton Line Search", 1},
                                                     {"Newton Trust Region - Hook Step", 2},
                                                     {"Newton Trust Region - Double Dog-Leg", 3},
                                                     {"Conjugate Gradient", 4},
                                                     {"BFGS", 5},
                                                     {"Limited Memory BFGS", 6},
                                                     {"BFGS Using More Thuente Method", 7}};

        bool consistencyCheck();
        auto_arima_object initAutoArima(int N_data);

    // last to be initialized
    public:
        auto_arima_object auto_arima_forward;
        auto_arima_object auto_arima_backward;

    };


} // namespace mio

#endif // INTERPOLARIMA_H