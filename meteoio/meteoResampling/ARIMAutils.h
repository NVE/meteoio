#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <meteoio/MeteoIO.h>
#include <cassert>


namespace mio {

// slice a vector from start to start+N
std::vector<double> slice(const std::vector<double> &vec, size_t start, int N);

// slice a vector from start to end
std::vector<double> slice(const std::vector<double> &vec, size_t start);

// np.arange for c++
std::vector<double> arange(size_t start, int N);

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
    if (vec.empty()) {
        throw std::invalid_argument("Cannot reverse an empty vector");
    }
    int start = 0;
    int end = int(vec.size()) - 1;

    while (start < end) {
        std::swap(vec[start], vec[end]);
        start++;
        end--;
    }
}

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(std::vector<MeteoData> vecM, const std::string &paramname);

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(std::vector<MeteoData> vecM, const size_t &paramindex);

// helper to parse direction argument for interpolarima
std::vector<double> decideDirection(std::vector<double> data, std::string direction, bool forward, size_t gap_loc);

// a class to cache information about a gap
struct ARIMA_GAP {
        ARIMA_GAP() : start(), end(), startDate(), endDate(), sampling_rate() {}
        void extend(const size_t& idx, const std::vector<MeteoData>& vecM) {if (idx<start) setStart(idx, vecM); if (idx>end) setEnd(idx, vecM);}
        void setStart(const size_t& idx, const std::vector<MeteoData>& vecM) {if (idx>=vecM.size()) return; start=idx; startDate=vecM[idx].date;}
        void setEnd(const size_t& idx, const std::vector<MeteoData>& vecM) {if (idx>=vecM.size()) return; end=idx; endDate=vecM[idx].date;}
        std::string toString() const {std::ostringstream os; os << startDate.toString(Date::ISO) << " (" << start << ") - " << endDate.toString(Date::ISO) << " (" << end << ")\n With sampling rate: " << sampling_rate ; return os.str();}
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

void computeARIMAGap(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date, size_t& indexP1, size_t& indexP2, double& before_window, double& after_window, double& window_size, Date& data_start_date, Date& data_end_date);

// roughly equal between two dates, given a tolerance level
bool requal(Date &date1, Date &date2);


// returns the most often accuring value in a vector
double mostLikelyValue(const std::vector<double>& vec);

// compute the most often occuring sampling rate rounded to 1e-6
double computeSamplingRate(Date data_start_date, Date data_end_date, std::vector<MeteoData> vecM);

template <typename T>
void printVectors(const std::vector<std::vector<T>>& vecs) {
    size_t maxSize = 0;
    for (const auto& vec : vecs) {
        maxSize = std::max(maxSize, vec.size());
    }

    // Print headers
    for (size_t i = 0; i < vecs.size(); i++) {
        std::cout << std::left << std::setw(10) << "Vector" + std::to_string(i+1);
    }
    std::cout << std::endl;
    std::cout << std::string(vecs.size() * 10, '-') << std::endl;

    for (size_t i = 0; i < maxSize; i++) {
        for (const auto& vec : vecs) {
            // Print elements from vec or "NaN" if out of range
            if (i < vec.size()) {
                std::cout << std::left << std::setw(10) << vec[i];
            } else {
                std::cout << std::left << std::setw(10) << "NaN";
            }
        }
        std::cout << std::endl;
    }
}

template <typename T>
void printVectors(const std::vector<Date>& vec1, const std::vector<T>& vec2) {
    size_t maxSize = std::max(vec1.size(), vec2.size());

    // Print headers
    std::cout << std::left << std::setw(30) << "Date1" << "| Date2" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    for (size_t i = 0; i < maxSize; i++) {
        // Print date from vec1 or "NaN" if out of range
        if (i < vec1.size()) {
            std::cout << std::left << std::setw(30) << vec1[i].toString(Date::ISO) << "| ";
        } else {
            std::cout << std::left << std::setw(30) << "NaN" << "| ";
        }

        // Print date from vec2 or "NaN" if out of range
        if (i < vec2.size()) {
            std::cout << vec2[i] << std::endl;
        } else {
            std::cout << "NaN" << std::endl;
        }
    }
}

void printVectors(const std::vector<MeteoData>& vec1, const std::vector<Date>& vec2, const size_t& paramindex);
void printVectors(const std::vector<Date>& vec1, const std::vector<Date>& vec2);
void printVector(const std::vector<MeteoData>& vec1);
void printVectors(const std::vector<double>& vec);

} // namespace mio
#endif //UTILS_H