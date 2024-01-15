#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cassert>
#include <meteoio/MeteoIO.h>
#include <cmath>

std::vector<double> slice(const std::vector<double> &vec, int start, int N) { 
    std::vector<double> vec_sliced(N);
    for (int i = 0; i < N; i++) {
        vec_sliced[i] = vec[start+i];
    }
    return vec_sliced;
};
std::vector<double> slice(const std::vector<double> &vec, int start) { 
    std::vector<double> vec_sliced(vec.size()-start);
    for (int i = 0; i < vec.size()-start; i++) {
        vec_sliced[i] = vec[start+i];
    }
    return vec_sliced;
};

std::vector<double> arange(int start, int N) {
    std::vector<double> vec(N);
    for (int i = 0; i < N; i++) {
        vec[i] = start + i;
    }
    return vec;
};

template <typename T>
T findMinMax(const std::vector<T>& vec, bool findMin = true) {
    assert(!vec.empty()); // Ensure the vector is not empty

    T extremeValue = vec[0];
    for(const auto& value : vec) {
        if(findMin ? value < extremeValue : value > extremeValue) {
            extremeValue = value;
        }
    }
    return extremeValue;
}

//calculate the mean of a vector
double mean(std::vector<double> vec) {
    double sum = 0;
    for (int i = 0; i < vec.size(); i++) {
        sum += vec[i];
    }
    return sum/vec.size();
}

//calculate the standard deviation of a vector
double stdDev(std::vector<double> vec) {
    double mean_vec = mean(vec);
    double sum = 0;
    for (int i = 0; i < vec.size(); i++) {
        sum += (vec[i] - mean_vec)*(vec[i] - mean_vec);
    }
    return std::sqrt(sum/vec.size());
}

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

std::vector<doubl> toVector(std::vector<MeteoData> vecM, const std::string &paramname) {
    int paramindex = vecM[0].getParamIndex(paramname);
    std::vector<double> vec(vecM.size());
    for (int i = 0; i < vecM.size(); i++) {
        vec[i] = vecM[i](paramindex);
    }
    return vec;
}

std::vector<double> decideDirection(std::vector<double> data, std::string direction, bool forward) {
    if (forward) {
        if (direction == "forward") {
            return data;
        } else if (direction == "backward") {
            reverseVector(data);
            return data;
        } else {
            throw mio::IOException("Direction " + direction + " not recognized");
        }
    } else {
        return std::vector<double>();
    }
}
#endif //UTILS_H