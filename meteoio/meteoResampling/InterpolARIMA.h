#ifndef INTERPOLARIMA_H
#define INTERPOLARIMA_H

#include <map>
#include <meteoio/thirdParty/ctsa.h>
#include <string>
#include <vector>
#include <iostream>
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
        std::vector<double> predict(int n_steps = 0);
        std::vector<double> getData() { return data; };
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