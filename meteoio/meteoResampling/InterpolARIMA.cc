#include "InterpolARIMA.h"
#include "ARIMAutils.h"
#include <cmath>
#include <cstdlib> // for std::rand and std::srand
#include <cstring>
#include <iostream>
#include <meteoio/MeteoIO.h>
#include <unistd.h>

namespace mio {
    // ------------------- Constructor ------------------- //
    // Default constructor
    InterpolARIMA::InterpolARIMA()
        : data(5, 0.0), gap_loc(0), N_gap(5), time(), pred_forward(), pred_backward(), xreg_vec_f(), xreg_vec_b(), data_forward(),
          data_backward(), new_xreg_vec_f(), new_xreg_vec_b(), xreg_f(NULL), xreg_b(NULL), new_xreg_f(NULL), new_xreg_b(NULL),
          amse_forward(), amse_backward(), N_data_forward(5), N_data_backward(5) {
        // dont initialize auto arima objects to not accidentally use "empty ones"
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_data_forward);
        auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_data_backward);
    }

    // only need s when it is known
    InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t gap_location, int gap_length, int period)
        : data(data_in), gap_loc(gap_location), N_gap(gap_length), time(arange(0, N_gap)), pred_forward(N_gap), pred_backward(N_gap),
          xreg_vec_f(0), xreg_vec_b(0), data_forward(slice(data, 0, static_cast<int>(gap_loc))),
          data_backward(slice(data, gap_loc + static_cast<size_t>(N_gap))), new_xreg_vec_f(0), new_xreg_vec_b(0), xreg_f(NULL),
          xreg_b(NULL), new_xreg_f(NULL), new_xreg_b(NULL), amse_forward(N_gap), amse_backward(N_gap), N_data_forward(data_forward.size()),
          N_data_backward(data_backward.size()), s(period) {
        // reverse the backward data
        reverseVector(data_backward); // TODO: can be done with std::reverse in C++17

        // initialize auto_arima objects
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        std::vector<int> pqdmax_b = {max_p, max_d, max_q};
        std::vector<int> PQDmax_b = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_data_forward);
        auto_arima_backward = auto_arima_init(pqdmax_b.data(), PQDmax_b.data(), s, r, N_data_backward);
    }

    InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t gap_location, int gap_length, std::vector<double> xreg_vec_in,
                                 int period)
        : data(data_in), gap_loc(gap_location), N_gap(gap_length), time(arange(0, N_gap)), pred_forward(N_gap), pred_backward(N_gap),
          xreg_vec_f(slice(xreg_vec_in, 0, gap_loc)), xreg_vec_b(reverseVectorReturn(slice(xreg_vec_in, gap_loc + N_gap))),
          data_forward(slice(data, 0, gap_loc)), data_backward(slice(data, gap_loc + N_gap)),
          new_xreg_vec_f(xreg_vec_f.size() == 0 ? 0 : N_gap), new_xreg_vec_b(xreg_vec_b.size() == 0 ? 0 : N_gap),
          xreg_f(xreg_vec_f.size() == 0 ? NULL : &xreg_vec_f[0]), xreg_b(xreg_vec_b.size() == 0 ? NULL : &xreg_vec_b[0]),
          new_xreg_f(xreg_vec_f.size() == 0 ? NULL : &new_xreg_vec_f[0]), new_xreg_b(xreg_vec_b.size() == 0 ? NULL : &new_xreg_vec_b[0]),
          amse_forward(N_gap), amse_backward(N_gap), N_data_forward(data_forward.size()), N_data_backward(data_backward.size()),
          r(xreg_vec_in.size() == 0 ? 0 : xreg_vec_f.size() / (N_data_forward)), s(period) {
        // reverse the backward data
        reverseVector(data_backward); // TODO: Can be done with std::reverse in C++17

        // initialize auto_arima objects
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        std::vector<int> pqdmax_b = {max_p, max_d, max_q};
        std::vector<int> PQDmax_b = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_data_forward);
        auto_arima_backward = auto_arima_init(pqdmax_b.data(), PQDmax_b.data(), s, r, N_data_backward);
    }

    InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t data_end, int n_predictions, std::string direction, int period)
        : data(data_in), gap_loc(data_end), N_gap(n_predictions), time(arange(0, static_cast<int>(data.size()))),
          pred_forward(n_predictions), pred_backward(n_predictions), xreg_vec_f(0), xreg_vec_b(0),
          data_forward(decideDirection(data_in, direction, true, gap_loc, n_predictions)),
          data_backward(decideDirection(data_in, direction, false, gap_loc, n_predictions)), new_xreg_vec_f(0), new_xreg_vec_b(0),
          xreg_f(NULL), xreg_b(NULL), new_xreg_f(NULL), new_xreg_b(NULL), amse_forward(N_gap), amse_backward(N_gap),
          N_data_forward(data_forward.size()), N_data_backward(data_backward.size()), s(period) {
        // initialize auto_arima objects

        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        std::vector<int> pqdmax_b = {max_p, max_d, max_q};
        std::vector<int> PQDmax_b = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_data_forward);
        auto_arima_backward = auto_arima_init(pqdmax_b.data(), PQDmax_b.data(), s, r, N_data_backward);
    }

    std::string InterpolARIMA::toString() {
        std::stringstream ss;

        ss << "\n---------------- Auto ARIMA Model Information ----------------\n";

        // Base Arima variables
        ss << "\nBase Arima Variables:\n";
        ss << std::left << std::setw(10) << "Max p:" << max_p << "\n";
        ss << std::left << std::setw(10) << "Max d:" << max_d << "\n";
        ss << std::left << std::setw(10) << "Max q:" << max_q << "\n";
        ss << std::left << std::setw(10) << "Start p:" << start_p << "\n";
        ss << std::left << std::setw(10) << "Start q:" << start_q << "\n";

        ss << "\n-------------------------------------------------------------\n";

        // Seasonal Arima variables
        ss << "\nSeasonal Arima Variables:\n";
        ss << std::left << std::setw(10) << "Max P:" << max_P << "\n";
        ss << std::left << std::setw(10) << "Max D:" << max_D << "\n";
        ss << std::left << std::setw(10) << "Max Q:" << max_Q << "\n";
        ss << std::left << std::setw(10) << "Start P:" << start_P << "\n";
        ss << std::left << std::setw(10) << "Start Q:" << start_Q << "\n";
        ss << std::left << std::setw(10) << "r:" << r << "\n";
        ss << std::left << std::setw(10) << "s:" << s << "\n";

        ss << "\n-------------------------------------------------------------\n";

        // Data information
        ss << "\nData Information:\n";
        ss << std::left << std::setw(20) << "Data size:" << data.size() << "\n";
        ss << std::left << std::setw(20) << "Gap location:" << gap_loc << "\n";
        ss << std::left << std::setw(20) << "N_gap:" << N_gap << "\n";
        ss << std::left << std::setw(20) << "Data forward size:" << data_forward.size() << "\n";
        ss << std::left << std::setw(20) << "Data backward size:" << data_backward.size() << "\n";

        ss << "\n-------------------------------------------------------------\n";
        ss << "Forward Model:\n";
        auto_arima_summary(auto_arima_forward);
        ss << "\n-------------------------------------------------------------\n";
        ss << "Backward Model:\n";
        auto_arima_summary(auto_arima_backward);

        ss << "\n-------------------------------------------------------------\n";

        return ss.str();
    }

    // Set the metadata for the auto arima objects
    void InterpolARIMA::setAutoArimaMetaData(int max_p_param, int max_d_param, int max_q_param, int start_p_param, int start_q_param,
                                             int max_P_param, int max_D_param, int max_Q_param, int start_P_param, int start_Q_param,
                                             bool seasonal_param, bool stationary_param) {
        this->max_p = max_p_param;
        this->max_d = max_d_param;
        this->max_q = max_q_param;
        this->start_p = start_p_param;
        this->start_q = start_q_param;
        this->max_P = max_P_param;
        this->max_D = max_D_param;
        this->max_Q = max_Q_param;
        this->start_P = start_P_param;
        this->start_Q = start_Q_param;
        this->seasonal = seasonal_param;
        this->stationary = stationary_param;
        auto_arima_backward->pmax = max_p;
        auto_arima_forward->pmax = max_p;
        auto_arima_backward->dmax = max_d;
        auto_arima_forward->dmax = max_d;
        auto_arima_backward->qmax = max_q;
        auto_arima_forward->qmax = max_q;
        auto_arima_backward->Pmax = max_P;
        auto_arima_forward->Pmax = max_P;
        auto_arima_backward->Dmax = max_D;
        auto_arima_forward->Dmax = max_D;
        auto_arima_backward->Qmax = max_Q;
        auto_arima_forward->Qmax = max_Q;
        auto_arima_backward->p_start = start_p;
        auto_arima_forward->p_start = start_p;
        auto_arima_backward->q_start = start_q;
        auto_arima_forward->q_start = start_q;
        auto_arima_backward->P_start = start_P;
        auto_arima_forward->P_start = start_P;
        auto_arima_backward->Q_start = start_Q;
        auto_arima_forward->Q_start = start_Q;
        auto_arima_backward->seasonal = seasonal;
        auto_arima_forward->seasonal = seasonal;
        auto_arima_backward->stationary = stationary;
        auto_arima_forward->stationary = stationary;
    }

    // Set the metadata for the auto arima objects optimization
    // options for method: "css-mle", "ml", "css"
    // options for opt_method: "Nelder-Mead", "Newton Line Search", "Newton Trust Region - Hook Step", "Newton Trust Region - Double
    // Dog-Leg", "Conjugate Gradient", "BFGS", "Limited Memory BFGS", "BFGS Using More Thuente Method"
    void InterpolARIMA::setOptMetaData(std::string method_param, std::string opt_method_param, bool stepwise_param,
                                       bool approximation_param, int num_models_param) {
        this->method = method_param;
        this->opt_method = opt_method_param;
        this->stepwise = stepwise_param;
        this->approximation = approximation_param;
        this->num_models = num_models_param;
        auto_arima_backward->method = method_map[method];
        auto_arima_forward->method = method_map[method];
        auto_arima_backward->optmethod = opt_method_map[opt_method];
        auto_arima_forward->optmethod = opt_method_map[opt_method];
        auto_arima_backward->stepwise = stepwise;
        auto_arima_forward->stepwise = stepwise;
    }

    // ------------------- Getters ------------------- //
    // Get the interpolated data
    std::vector<double> InterpolARIMA::getInterpolatedData() {
        std::vector<double> interpolated_data(data.begin() + gap_loc, data.begin() + gap_loc + N_gap);
        return interpolated_data;
    }

    // ------------------- Interpolation methods ------------------- //
    // Simulate n_steps into the future
    std::vector<double> InterpolARIMA::simulate(int n_steps, int seed) {
        std::vector<double> sim(n_steps);
        seed++;
        // use the equations to simulate with random errors
        std::cerr << "not implemented, and not needed for now\n";
        return sim;
    }

    std::vector<double> InterpolARIMA::predict(int n_steps) {
        if (n_steps == 0) {
            n_steps = N_gap;
        } else {
            pred_forward.resize(n_steps);
            amse_forward.resize(n_steps);
        }


        auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
        // check if the models are valid (p and q should not be zero at the same time)
        if ((auto_arima_forward->p == 0 && auto_arima_forward->q == 0 ) && (auto_arima_forward->P == 0 && auto_arima_forward->Q == 0)) {
            bool current_stepwise = auto_arima_forward->stepwise;
            auto_arima_setStepwise(auto_arima_forward, !current_stepwise);
            auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
        }
        auto_arima_predict(auto_arima_forward, data_forward.data(), xreg_f, n_steps, new_xreg_f, pred_forward.data(), amse_forward.data());
        return pred_forward;
    }

    // Fill the gap using the auto arima objects
    void InterpolARIMA::fillGap() {
        bool isRandom_f;
        bool isRandom_b;
        for (int meth_Id = 0; meth_Id <3 ; meth_Id++) {
            // fit the models
            auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
            auto_arima_exec(auto_arima_backward, data_backward.data(), xreg_b);

            isRandom_b = false;
            isRandom_f = false;
            // weighting should not be done if one of the models ends up being a random walk
            if ((auto_arima_forward->p == 0 && auto_arima_forward->q == 0) && (auto_arima_forward->P == 0 && auto_arima_forward->Q == 0)) {
                isRandom_f = true;
            } 
            if ((auto_arima_backward->p == 0 && auto_arima_backward->q == 0) && (auto_arima_backward->P == 0 && auto_arima_backward->Q == 0)){
                isRandom_b = true;
            }

            if (isRandom_b && isRandom_f && meth_Id == 0) {
                isRandom_b = false;
                isRandom_f = false;
                bool current_stepwise = auto_arima_forward->stepwise;
                auto_arima_setStepwise(auto_arima_forward, !current_stepwise);
                auto_arima_setStepwise(auto_arima_backward, !current_stepwise);
                auto_arima_exec(auto_arima_forward, data_forward.data(), xreg_f);
                auto_arima_exec(auto_arima_backward, data_backward.data(), xreg_b);
                if ((auto_arima_forward->p == 0 && auto_arima_forward->q == 0) && (auto_arima_forward->P == 0 && auto_arima_forward->Q == 0)) {
                    isRandom_f = true;
                } 
                if ((auto_arima_backward->p == 0 && auto_arima_backward->q == 0) && (auto_arima_backward->P == 0 && auto_arima_backward->Q == 0)){
                    isRandom_b = true;
                }
            }
            if (isRandom_b && isRandom_f)
                break;
            // predict the gap
            // forward
            auto_arima_predict(auto_arima_forward, data_forward.data(), xreg_f, N_gap, new_xreg_f, pred_forward.data(),
                                amse_forward.data());

            // backward
            auto_arima_predict(auto_arima_backward, data_backward.data(), xreg_b, N_gap, new_xreg_b, pred_backward.data(),
                                amse_backward.data());
            // interpolate with the weighting according to
            assert(pred_forward.size() == pred_backward.size());
            assert(pred_forward.size() == N_gap);

            reverseVector(pred_backward); // TODO: can be done with std::reverse in C++17

            bool consistency = consistencyCheck();
            std::cout << "consistency check " << consistency << std::endl;
            if (consistency)
                break;
        }
        // W1 = sqrt(T-t/T)
        // W2 = sqrt(t/T)
        for (int id = 0; id < N_gap; id++) {
            double weight_f, weight_b;

            if (isRandom_f) {
                weight_f = 0;
                weight_b = 1;
            } else if (isRandom_b) {
                weight_f = 1;
                weight_b = 0;
            } else {
                weight_f = std::sqrt((N_gap - time[id]) / N_gap);
                weight_b = std::sqrt(time[id] / N_gap);
            }
            data[gap_loc + id] = weight_f * pred_forward[id] + weight_b * pred_backward[id];
        }
    }

    bool InterpolARIMA::consistencyCheck() {
        double mean_before = calcVecMean(data_forward);
        double std_before = stdDev(data_forward);
        double mean_after = calcVecMean(data_backward);
        double std_after = stdDev(data_backward);
        double max_interpolated = findMinMax(slice(data, gap_loc, N_gap), false);

        std::cout << "mean before: " << mean_before << std::endl;
        std::cout << "std before: " << std_before << std::endl;
        std::cout << "mean after: " << mean_after << std::endl;
        std::cout << "std after: " << std_after << std::endl;
        std::cout << "max interpolated: " << max_interpolated << std::endl;

        // needs more checks
        if (max_interpolated > mean_before + 4 * std_before && max_interpolated > mean_after + 4 * std_after)
            return false;
        else if (max_interpolated < mean_before - 4 * std_before && max_interpolated < mean_after - 4 * std_after)
            return false;
        return true;
    }

    void InterpolARIMA::interpolate() {
        bool fit = true;
        if (N_data_backward == 0 || N_data_forward == 0) {
            throw NoDataException("No data to interpolate: forward datapoints " + std::to_string(N_data_forward) +
                                  ", backward datapoints " + std::to_string(N_data_backward) + "\n");
            return;
        }
        while (fit) {
            fillGap();
            int retval_f = auto_arima_forward->retval;
            int retval_b = auto_arima_backward->retval;

            if (retval_f == 0 || retval_b == 0) {
                std::string where = (retval_f == 0) ? "forward data" : "backward data";
                throw AccessException("Interpolation Input data is erroneous in " + where);
            } else if (retval_f == 15 || retval_b == 15) {
                std::string where = (retval_f == 15) ? "forward data" : "backward data";
                throw InvalidFormatException("Interpolation Input data has Inf/Nan values in " + where);
            } else if (retval_f == 4 || retval_b == 4) {
                std::string where = (retval_f == 4) ? "forward data" : "backward data";
                if (method != "CSS-MLE" && opt_method != "BFGS") {
                    throw IOException("Optimization of ARIMA did not converge in " + where +
                                      ".\n Please try another method and optimization method");
                } else {
                    std::string new_opt_method = "Nelder-Mead";
                    setOptMetaData(method, new_opt_method);
                }
            } else {
                fit = false;
            }
        }
    }

} // end namespace mio
