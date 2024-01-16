#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <meteoio/MeteoIO.h>
#include <cassert>


namespace mio {

// slice a vector from start to start+N
std::vector<double> slice(const std::vector<double> &vec, int start, int N);

// slice a vector from start to end
std::vector<double> slice(const std::vector<double> &vec, int start);

// np.arange for c++
std::vector<double> arange(int start, int N);

template <typename T>
T findMinMax(const std::vector<T>& vec, bool findMin) {
    assert(!vec.empty()); // Ensure the vector is not empty

    T extremeValue = vec[0];
    for(const auto& value : vec) {
        if(findMin ? value < extremeValue : value > extremeValue) {
            extremeValue = value;
        }
    }
    return extremeValue;
}

//calculate the of a vector
double calcVecMean(std::vector<double> vec);

//calculate the standard deviation of a vector
double stdDev(std::vector<double> vec);

// TODO: Test this
template <typename T>
void reverseVector(std::vector<T>& vec) {
    int start = 0;
    int end = vec.size() - 1;

    while (start < end) {
        std::swap(vec[start], vec[end]);
        start++;
        end--;
    }
}

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(std::vector<MeteoData> vecM, const std::string &paramname);

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(std::vector<MeteoData> vecM, const int &paramindex);

// helper to parse direction argument for interpolarima
std::vector<double> decideDirection(std::vector<double> data, std::string direction, bool forward);

// a class to cache information about a gap
struct ARIMA_GAP {
        ARIMA_GAP() : start(), end(), startDate(), endDate(), sampling_rate() {}
        void extend(const size_t& idx, const std::vector<MeteoData>& vecM) {if (idx<start) setStart(idx, vecM); if (idx>end) setEnd(idx, vecM);}
        void setStart(const size_t& idx, const std::vector<MeteoData>& vecM) {if (idx>=vecM.size()) return; start=idx; startDate=vecM[idx].date;}
        void setEnd(const size_t& idx, const std::vector<MeteoData>& vecM) {if (idx>=vecM.size()) return; end=idx; endDate=vecM[idx].date;}
        std::string toString() const {std::ostringstream os; os << startDate.toString(Date::ISO) << " (" << start << ") - " << endDate.toString(Date::ISO) << " (" << end << ")"; return os.str();}
        void reset() {start=IOUtils::npos; end=IOUtils::npos; startDate=Date(); endDate=Date();}
        size_t start, end;
        Date startDate, endDate;
        double sampling_rate;
        bool isGap(){return (endDate-startDate).getJulian(true)*sampling_rate > 3;};
};

//return true if a valid point could be found backward from pos
size_t searchBackward(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size);

//return true if a valid point could be found forward from pos
size_t searchForward(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size, const size_t& indexP1);

void computeGap(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size, size_t& indexP1, size_t& indexP2);

// roughly equal between two dates, given a tolerance level
static bool requal(Date &date1, Date &date2) {
    double tolerance = 1e-6; // Define your tolerance level
    return std::abs((date1 - date2).getJulian(true)) <= tolerance;
}


// returns the most often accuring value in a vector
double mostLikelyValue(const std::vector<double>& vec);

// compute the most often occuring sampling rate rounded to 1e-6
double computeSamplingRate(Date data_start_date, Date data_end_date, std::vector<MeteoData> vecM);

} // namespace mio
#endif //UTILS_H