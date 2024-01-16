#include "InterpolARIMA.h"
#include "ARIMAutils.h"
#include <cmath>
#include <iostream>
#include <meteoio/MeteoIO.h>

namespace mio {
// ------------------- Constructor ------------------- //
// Default constructor
InterpolARIMA::InterpolARIMA()
    : gap_loc(0), N_gap(0), time(), pred_forward(), pred_backward(), data(), xreg_vec(), 
    data_forward(), data_backward(), new_xreg_vec(), xreg(NULL), new_xreg(NULL), 
    amse_forward(), amse_backward(), N_data_forward(0), N_data_backward(0), max_p(8), 
    max_d(3), max_q(8), start_p(2), start_q(2), max_P(2), max_D(1), max_Q(2), start_P(1), 
    start_Q(1), r(0), s(0), method("css-mle"), opt_method("bfgs"), stepwise(true), 
    approximation(false), num_models(94), seasonal(false), stationary(false) 
    {// initialize auto_arima objects
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
        auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
    }

// only need s when it is known
InterpolARIMA::InterpolARIMA(std::vector<double> data, int gap_loc, int N_gap, int s)
    : data(data), time(arange(0, N_gap)), gap_loc(gap_loc), N_gap(N_gap), data_forward(slice(data, 0, gap_loc)),
      N_data_forward(data_forward.size()), data_backward(slice(data, gap_loc + N_gap - 1)), N_data_backward(data_backward.size()),
      xreg_vec(0), new_xreg_vec(0), xreg(NULL), r(0), s(s), pred_forward(N_gap), pred_backward(N_gap), amse_forward(N_gap),
      amse_backward(N_gap) {
    // reverse the backward data
    reverseVector(data_backward); // TODO: can be done with std::reverse in C++17

    // initialize auto_arima objects
    std::vector<int> pqdmax = {max_p, max_d, max_q};
    std::vector<int> PQDmax = {max_P, max_D, max_Q};
    auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
    auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
}

InterpolARIMA::InterpolARIMA(std::vector<double> data, int gap_loc, int N_gap, std::vector<double> xreg_vec, int s)
    : data(data), time(arange(0, N_gap)), gap_loc(gap_loc), N_gap(N_gap), xreg_vec(xreg_vec), data_forward(slice(data, 0, gap_loc)),
      N_data_forward(data_forward.size()), data_backward(slice(data, gap_loc + N_gap - 1)), N_data_backward(data_backward.size()),
      xreg(xreg_vec.size() == 0 ? NULL : &xreg_vec[0]), new_xreg_vec(xreg_vec.size() == 0 ? 0 : N_gap),
      new_xreg(xreg_vec.size() == 0 ? NULL : &new_xreg_vec[0]),
      r(xreg_vec.size() == 0 ? 0 : xreg_vec.size() / (N_data_backward + N_data_forward)), s(s), pred_forward(N_gap), pred_backward(N_gap),
      amse_forward(N_gap), amse_backward(N_gap) {
    // reverse the backward data
    reverseVector(data_backward); // TODO: Can be done with std::reverse in C++17

    // initialize auto_arima objects
    std::vector<int> pqdmax = {max_p, max_d, max_q};
    std::vector<int> PQDmax = {max_P, max_D, max_Q};
    auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
    auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
}


InterpolARIMA::InterpolARIMA(std::vector<double> data, int n_predictions, std::string direction, int s)
    : data(data), time(arange(0, data.size())), gap_loc(data.size()-1), N_gap(n_predictions), data_forward(decideDirection(data, direction, true)),
      N_data_forward(data_forward.size()), data_backward(decideDirection(data,direction,false)), N_data_backward(data_backward.size()),
      xreg_vec(0), new_xreg_vec(0), xreg(NULL), r(0), s(s), pred_forward(N_gap), pred_backward(N_gap), amse_forward(N_gap),
      amse_backward(N_gap) {
        // initialize auto_arima objects
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
        auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, 0);
      }

// Set the metadata for the auto arima objects
void InterpolARIMA::setAutoArimaMetaData(int max_p, int max_d, int max_q, int start_p, int start_q, int max_P, int max_D, int max_Q,
                                         int start_P, int start_Q, bool seasonal, bool stationary) {
    this->max_p = max_p;
    this->max_d = max_d;
    this->max_q = max_q;
    this->start_p = start_p;
    this->start_q = start_q;
    this->max_P = max_P;
    this->max_D = max_D;
    this->max_Q = max_Q;
    this->start_P = start_P;
    this->start_Q = start_Q;
    this->seasonal = seasonal;
    this->stationary = stationary;
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
// options for opt_method: "Nelder-Mead", "Newton Line Search", "Newton Trust Region - Hook Step", "Newton Trust Region - Double Dog-Leg",
// "Conjugate Gradient", "BFGS", "Limited Memory BFGS", "BFGS Using More Thuente Method"
void InterpolARIMA::setOptMetaData(std::string method, std::string opt_method, bool stepwise, bool approximation, int num_models) {
    this->method = method;
    this->opt_method = opt_method;
    this->stepwise = stepwise;
    this->approximation = approximation;
    this->num_models = num_models;
    auto_arima_backward->method = method_map[method];
    auto_arima_forward->method = method_map[method];
    auto_arima_backward->optmethod = opt_method_map[opt_method];
    auto_arima_forward->optmethod = opt_method_map[opt_method];
    auto_arima_backward->stepwise = stepwise;
    auto_arima_forward->stepwise = stepwise;
}

// ------------------- Interpolation methods ------------------- //
// Simulate n_steps into the future
std::vector<double> InterpolARIMA::simulate(int n_steps, int seed) {
    std::vector<double> sim(n_steps);
    // use the equations to simulate with random errors
    std::cerr << "not implemented, and not needed for now\n";
    return sim;
}

std::vector<double> InterpolARIMA::predict() {
    auto_arima_exec(auto_arima_forward, data_forward.data(), xreg);
    auto_arima_predict(auto_arima_forward, data_forward.data(), xreg, N_gap, new_xreg, pred_forward.data(), amse_forward.data());
    return pred_forward;
}

// Fill the gap using the auto arima objects
void InterpolARIMA::fillGap() {
    int max_iter = 5;
    for (int i = 0; i < max_iter; i++) {
        // fit the models
        auto_arima_exec(auto_arima_forward, data_forward.data(), xreg);
        auto_arima_exec(auto_arima_backward, data_backward.data(), xreg);

        // predict the gap
        // forward
        auto_arima_predict(auto_arima_forward, data_forward.data(), xreg, N_gap, new_xreg, pred_forward.data(), amse_forward.data());

        // backward
        auto_arima_predict(auto_arima_backward, data_backward.data(), xreg, N_gap, new_xreg, pred_backward.data(), amse_backward.data());
        // interpolate with the weighting according to
        assert(pred_forward.size() == pred_backward.size());

        reverseVector(pred_backward); // TODO: can be done with std::reverse in C++17

        // W1 = sqrt(T-t/T)
        // W2 = sqrt(t/T)
        for (int i = 0; i < N_gap; i++) {
            double weight_f = std::sqrt((N_gap - time[i]) / N_gap);
            double weight_b = std::sqrt(time[i] / N_gap);
            data[gap_loc + i] = weight_f * pred_forward[i] + weight_b * pred_backward[i];
        }
        if (consistencyCheck())
            break;
    }
}

bool InterpolARIMA::consistencyCheck() {
    double mean_before = calcVecMean(data_forward);
    double std_before = stdDev(data_forward);
    double mean_after = calcVecMean(data_backward);
    double std_after = stdDev(data_backward);
    double max_interpolated = findMinMax(slice(data, gap_loc, N_gap), false);

    if (max_interpolated>mean_before+2*std_before && max_interpolated>mean_after+2*std_after)
        return false;
    return true;
}

void InterpolARIMA::interpolate() {
    bool fit = true;
    if (N_data_backward == 0 || N_data_forward == 0) {
        throw mio::NoDataException("No data to interpolate");
        return;
    }
    while (fit) {
        fillGap();
        int retval_f = auto_arima_forward->retval;
        int retval_b = auto_arima_backward->retval;

        if (retval_f == 0 || retval_b == 0) {
            std::string where = (retval_f == 0) ? "forward data" : "backward data";
            throw mio::AccessException("Interpolation Input data is erroneous in " + where);
        } else if (retval_f == 15 || retval_b == 15) {
            std::string where = (retval_f == 15) ? "forward data" : "backward data";
            throw mio::InvalidFormatException("Interpolation Input data has Inf/Nan values in " + where);
        } else if (retval_f == 4 || retval_b == 4) {
            std::string where = (retval_f == 4) ? "forward data" : "backward data";
            if (method != "css-mle" && opt_method != "BFGS") {
                throw mio::IOException("Optimization of ARIMA did not converge in " + where +
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
