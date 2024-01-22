#ifndef INTERPOLARIMA_H
#define INTERPOLARIMA_H

#include <meteoio/thirdParty/ctsa.h>
#include <map>
#include <string>
#include <vector>

namespace mio {


    // It is assumed that data is a vector of length N_data_forward + N_gap + N_data_backward, i.e. containing missing values, or similar in
    // the gap with linearly spaced entries in time
    class InterpolARIMA {
    public:
        InterpolARIMA();
        InterpolARIMA(std::vector<double> data_in, size_t gap_loc, int N_gap, int s = 0);
        InterpolARIMA(std::vector<double> data_in, size_t gap_loc, int N_gap, std::vector<double> xreg_vec, int s = 0);
        InterpolARIMA(std::vector<double> data_in, size_t gap_loc, int n_predictions, std::string direction = "forward", int s = 0);

        void setAutoArimaMetaData(int max_p_param = 8, int max_d_param = 3, int max_q = 8, int start_p = 2, int start_q = 2, int max_P = 2,
                                  int max_D = 1, int max_Q = 2, int start_P = 1, int start_Q = 1, bool seasonal = true,
                                  bool stationary = false);
        void setOptMetaData(std::string method = "css-mle", std::string opt_method = "BFGS", bool stepwise = true,
                            bool approximation = false, int num_models = 94);

        auto_arima_object auto_arima_forward;
        auto_arima_object auto_arima_backward;

        std::vector<double> simulate(int n_steps, int seed = 0);
        void fillGap();
        void interpolate();
        std::vector<double> predict();
        std::vector<double> getData() { return data; };
        std::vector<double> getInterpolatedData();


        // Swap function
        void swap(InterpolARIMA& first, InterpolARIMA& second) 
        {
            using std::swap;
            std::swap(first.auto_arima_forward, second.auto_arima_forward);
            std::swap(first.auto_arima_backward, second.auto_arima_backward);
            std::swap(first.xreg_f, second.xreg_f);
            std::swap(first.xreg_b, second.xreg_b);
            std::swap(first.new_xreg_f, second.new_xreg_f);
            std::swap(first.new_xreg_b, second.new_xreg_b);
            std::swap(first.gap_loc, second.gap_loc);
            std::swap(first.N_gap, second.N_gap);
            std::swap(first.time, second.time);
            std::swap(first.pred_forward, second.pred_forward);
            std::swap(first.pred_backward, second.pred_backward);
            std::swap(first.data, second.data);
            std::swap(first.xreg_vec_f, second.xreg_vec_f);
            std::swap(first.xreg_vec_b, second.xreg_vec_b);
            std::swap(first.data_forward, second.data_forward);
            std::swap(first.data_backward, second.data_backward);
            std::swap(first.new_xreg_vec_f, second.new_xreg_vec_f);
            std::swap(first.new_xreg_vec_b, second.new_xreg_vec_b);
            std::swap(first.N_data_forward, second.N_data_forward);
            std::swap(first.N_data_backward, second.N_data_backward);
            std::swap(first.max_p, second.max_p);
            std::swap(first.max_d, second.max_d);
            std::swap(first.max_q, second.max_q);
            std::swap(first.start_p, second.start_p);
            std::swap(first.start_q, second.start_q);
            std::swap(first.max_P, second.max_P);
            std::swap(first.max_D, second.max_D);
            std::swap(first.max_Q, second.max_Q);
            std::swap(first.start_P, second.start_P);
            std::swap(first.start_Q, second.start_Q);
            std::swap(first.r, second.r);
            std::swap(first.s, second.s);
            std::swap(first.method, second.method);
            std::swap(first.opt_method, second.opt_method);
            std::swap(first.stepwise, second.stepwise);
            std::swap(first.approximation, second.approximation);
            std::swap(first.num_models, second.num_models);
            std::swap(first.seasonal, second.seasonal);
            std::swap(first.stationary, second.stationary);

        }

        // Copy assignment operator
        InterpolARIMA& operator=(InterpolARIMA other) 
        {
            swap(*this, other);
            return *this;
        }

        // Copy constructor
        InterpolARIMA(const InterpolARIMA& other) 
        : InterpolARIMA() // Default-construct this object
        {
            *this = other; // Use copy assignment
        }

        // Destructor
        ~InterpolARIMA() {
            delete auto_arima_forward;
            delete auto_arima_backward;
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
        double* xreg_f;
        double* xreg_b;
        double* new_xreg_f;
        double* new_xreg_b;
        std::vector<double> amse_forward, amse_backward;
        size_t N_data_forward, N_data_backward;
        int max_p = 8, max_d = 3, max_q = 8;
        int start_p = 2, start_q = 2;
        int max_P = 2, max_D = 1, max_Q = 2;
        int start_P = 1, start_Q = 1;
        int r = 0, s = 0;
        std::string method = "CSS-MLE", opt_method = "BFGS";
        bool stepwise = true, approximation = false;
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
    };

} // namespace mio

#endif // INTERPOLARIMA_H