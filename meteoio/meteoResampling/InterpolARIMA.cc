#include "InterpolARIMA.h"
#include "ARIMAutils.h"
#include <cmath>
#include <iostream>
#include <meteoio/MeteoIO.h>

namespace mio {
// ------------------- Constructor ------------------- //
// Default constructor
InterpolARIMA::InterpolARIMA()
    :  auto_arima_forward(), auto_arima_backward(), data(), gap_loc(0), N_gap(0), time(), pred_forward(), pred_backward(), xreg_vec(), 
    data_forward(), data_backward(), new_xreg_vec(), xreg(NULL), new_xreg(NULL), 
    amse_forward(), amse_backward(), N_data_forward(0), N_data_backward(0)
    {// initialize auto_arima objects
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
        auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
    }

// only need s when it is known
InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t gap_location, int gap_length, int period)
    : auto_arima_forward(), auto_arima_backward(), data(data_in), gap_loc(gap_location), N_gap(gap_length), time(arange(0, N_gap)), pred_forward(N_gap), pred_backward(N_gap),
    xreg_vec(0), data_forward(slice(data, 0, static_cast<int>(gap_loc))), data_backward(slice(data, gap_loc + static_cast<size_t>(N_gap) - 1)), 
    new_xreg_vec(0), xreg(NULL), new_xreg(NULL), amse_forward(N_gap), amse_backward(N_gap), 
    N_data_forward(data_forward.size()), N_data_backward(data_backward.size()), s(period){
    // reverse the backward data
    reverseVector(data_backward); // TODO: can be done with std::reverse in C++17

    // initialize auto_arima objects
    std::vector<int> pqdmax = {max_p, max_d, max_q};
    std::vector<int> PQDmax = {max_P, max_D, max_Q};
    auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
    auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
}

InterpolARIMA::InterpolARIMA(std::vector<double> data_in, size_t gap_location, int gap_length, std::vector<double> xreg_vec_in, int period)
    : auto_arima_forward(), auto_arima_backward(), data(data_in), gap_loc(gap_location), N_gap(gap_length), time(arange(0, N_gap)), pred_forward(N_gap), pred_backward(N_gap),
    xreg_vec(xreg_vec_in), data_forward(slice(data, 0, gap_loc)), data_backward(slice(data, gap_loc + N_gap - 1)), 
    new_xreg_vec(xreg_vec.size() == 0 ? 0 : N_gap), xreg(xreg_vec.size() == 0 ? NULL : &xreg_vec[0]), 
    new_xreg(xreg_vec.size() == 0 ? NULL : &new_xreg_vec[0]), amse_forward(N_gap), amse_backward(N_gap), 
    N_data_forward(data_forward.size()), N_data_backward(data_backward.size()), r(xreg_vec.size() == 0 ? 0 : xreg_vec.size() / (N_data_backward + N_data_forward)), s(period) {
    // reverse the backward data
    reverseVector(data_backward); // TODO: Can be done with std::reverse in C++17

    // initialize auto_arima objects
    std::vector<int> pqdmax = {max_p, max_d, max_q};
    std::vector<int> PQDmax = {max_P, max_D, max_Q};
    auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
    auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
}


InterpolARIMA::InterpolARIMA(std::vector<double> data_in, int n_predictions, std::string direction, int period)
    :auto_arima_forward(), auto_arima_backward(), data(data_in), gap_loc(data.size()-1), N_gap(n_predictions), time(arange(0, static_cast<int>(data.size()))), pred_forward(n_predictions), pred_backward(n_predictions), 
     xreg_vec(0), data_forward(decideDirection(data, direction, true)), 
    data_backward(decideDirection(data,direction,false)), new_xreg_vec(0), xreg(NULL), new_xreg(NULL),
    amse_forward(N_gap), amse_backward(N_gap), N_data_forward(data_forward.size()), 
    N_data_backward(data_backward.size()), s(period) {
        // initialize auto_arima objects
        std::vector<int> pqdmax = {max_p, max_d, max_q};
        std::vector<int> PQDmax = {max_P, max_D, max_Q};
        auto_arima_forward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, N_gap);
        auto_arima_backward = auto_arima_init(pqdmax.data(), PQDmax.data(), s, r, 0);
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

    ss << "\n-------------------------------------------------------------\n";
    ss << "Backward Model:\n";

    ss << "\n-------------------------------------------------------------\n";

    return ss.str();
}

// Set the metadata for the auto arima objects
void InterpolARIMA::setAutoArimaMetaData(int max_p_param, int max_d_param, int max_q_param, int start_p_param, int start_q_param, int max_P_param, int max_D_param, int max_Q_param,
                                         int start_P_param, int start_Q_param, bool seasonal_param, bool stationary_param) {
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
// options for opt_method: "Nelder-Mead", "Newton Line Search", "Newton Trust Region - Hook Step", "Newton Trust Region - Double Dog-Leg",
// "Conjugate Gradient", "BFGS", "Limited Memory BFGS", "BFGS Using More Thuente Method"
void InterpolARIMA::setOptMetaData(std::string method_param, std::string opt_method_param, bool stepwise_param, bool approximation_param, int num_models_param) {
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

// ------------------- Interpolation methods ------------------- //
// Simulate n_steps into the future
std::vector<double> InterpolARIMA::simulate(int n_steps, int seed) {
    std::vector<double> sim(n_steps);
    seed++;
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
    for (int _i = 0; _i < max_iter; _i++) {
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
        for (int id = 0; id < N_gap; id++) {
            double weight_f = std::sqrt((N_gap - time[id]) / N_gap);
            double weight_b = std::sqrt(time[id] / N_gap);
            data[gap_loc + id] = weight_f * pred_forward[id] + weight_b * pred_backward[id];
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
            if (method != "CSS-MLE" && opt_method != "BFGS") {
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
