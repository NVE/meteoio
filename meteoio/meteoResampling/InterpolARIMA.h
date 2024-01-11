#ifndef INTERPOLARIMA_H
#define INTERPOLARIMA_H

#include <vector>
#include <string>
#include <map>
#include <ctsa.h>

// It is assumed that data is a vector of length N_data_forward + N_gap + N_data_backward, i.e. containing missing values, or similar in the gap
// with linearly spaced entries in time
class InterpolARIMA {
    public:
        InterpolARIMA();
        InterpolARIMA(std::vector<double> data, int gap_loc, int N_gap, int s = 0);
        InterpolARIMA(std::vector<double> data, int gap_loc, int N_gap, std::vector<double> xreg_vec, int s = 0);

        void setAutoArimaMetaData(int max_p = 8, int max_d =3, int max_q = 8, int start_p = 2, int start_q = 2, int max_P = 2, int max_D = 1, int max_Q = 2, int start_P = 1, int start_Q = 1, bool seasonal = false, bool stationary = false);
        void setOptMetaData(std::string method = "css-mle", std::string opt_method = "BFGS", bool stepwise = true, bool approximation = false, int num_models = 94);
        
        auto_arima_object auto_arima_forward;
        auto_arima_object auto_arima_backward;

        std::vector<double> simulate(int n_steps, int seed = 0);
        void fillGap();
        void interpolate();
        std::vector<double> getData() {return data;};


    private:
        // Interpolation variables
        const int gap_loc, N_gap;
        const std::vector<double> time;
        std::vector<double> pred_forward, pred_backward;

        // Auto Arima variables
        // const doesnt work wiht c
        std::vector<double> data, xreg_vec, data_forward, data_backward, new_xreg_vec;
        double* xreg, *new_xreg; 
        std::vector<double> amse_forward, amse_backward;
        const int N_data_forward, N_data_backward;
        int max_p = 8, max_d = 3, max_q = 8;
        int start_p = 2, start_q = 2;
        int max_P = 2, max_D = 1, max_Q = 2;
        int start_P = 1, start_Q = 1;
        int r = 0, s = 0;
        std::string method = "css-mle", opt_method = "bfgs";
        bool stepwise = true, approximation = false;
        int num_models = 94;
        bool seasonal = false, stationary = false;        

        std::map<std::string, int> method_map = {{"css-mle", 0}, {"mle", 1}, {"css", 2}};
        std::map<std::string, int> opt_method_map = {
            {"Nelder-Mead", 0}, 
            {"Newton Line Search", 1}, 
            {"Newton Trust Region - Hook Step", 2}, 
            {"Newton Trust Region - Double Dog-Leg", 3}, 
            {"Conjugate Gradient", 4}, 
            {"BFGS", 5}, 
            {"Limited Memory BFGS", 6}, 
            {"BFGS Using More Thuente Method", 7}
        };

        bool consistencyCheck();
};

#endif //INTERPOLARIMA_H