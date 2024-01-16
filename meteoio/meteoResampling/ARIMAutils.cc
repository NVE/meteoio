#include "ARIMAutils.h"
#include <cmath>
#include <map>


namespace mio {

// slice a vector from start to start+N
std::vector<double> slice(const std::vector<double> &vec, size_t start, int N) { 
    std::vector<double> vec_sliced(N);
    for (int i = 0; i < N; i++) {
        vec_sliced[i] = vec[start+i];
    }
    return vec_sliced;
};

// slice a vector from start to end
std::vector<double> slice(const std::vector<double> &vec, size_t start) { 
    std::vector<double> vec_sliced(vec.size()-start);
    for (size_t i = 0; i < vec.size()-start; i++) {	
        vec_sliced[i] = vec[start+i];
    }
    return vec_sliced;
};

// np.arange for c++
std::vector<double> arange(size_t start, int N) {
    std::vector<double> vec(N);
    for (size_t i = 0; i < static_cast<size_t>(N); i++) {
        vec[i] = start + i;
    }
    return vec;
};

//calculate the of a vector
double calcVecMean(std::vector<double> vec) {
    double sum = 0;
    for (size_t i = 0; i < vec.size(); i++) {
        sum += vec[i];
    }
    return sum/static_cast<int>(vec.size());
}

//calculate the standard deviation of a vector
double stdDev(std::vector<double> vec) {
	double mean_vec = calcVecMean(vec);
	double sum = 0;
	for (size_t i = 0; i < vec.size(); i++) {
		sum += (vec[i] - mean_vec)*(vec[i] - mean_vec);
	}
	return std::sqrt(sum/static_cast<double>(vec.size()));
}

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(std::vector<MeteoData> vecM, const std::string &paramname) {
    size_t paramindex = vecM[0].getParameterIndex(paramname);
    std::vector<double> vec(vecM.size());
    for (size_t i = 0; i < vecM.size(); i++) {
        vec[i] = vecM[i](paramindex);
    }
    return vec;
}

// converts a vector of MeteoData to a vector of doubles
std::vector<double> toVector(std::vector<MeteoData> vecM, const size_t &paramindex) {
    std::vector<double> vec(vecM.size());
    for (size_t i = 0; i < vecM.size(); i++) {
        vec[i] = vecM[i](paramindex);
    }
    return vec;
}

// helper to parse direction argument for interpolarima
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

//return true if a valid point could be found backward from pos
size_t searchBackward(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size)
{
	const Date windowStart( resampling_date - i_window_size );
	const bool knownGap = (!last_gap.startDate.isUndef() && !last_gap.endDate.isUndef());
	const bool currentInGap = (knownGap)? (resampling_date>=last_gap.startDate && resampling_date<=last_gap.endDate) : false;
	const bool windowStartInGap = (knownGap)? (windowStart>=last_gap.startDate && windowStart<=last_gap.endDate) : false;
	
	//the current point and window start are in a known gap, there is no hope
	if (currentInGap && windowStartInGap) return IOUtils::npos;
	
	//the current point is NOT in a known gap
	if (!currentInGap) { //or !knownGap
		const Date dateStart = (windowStartInGap)? last_gap.endDate : windowStart;
		size_t ii = pos; //because idx will get decremented right away
		for (; ii-- >0; ) {
			if (vecM[ii].date < dateStart) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) {
				last_gap.setStart(ii+1, vecM);
				last_gap.setEnd(pos, vecM);
				return ii;
			}
		}
		
		//no valid point found
		if (windowStartInGap) { //we can extend the current known gap
			last_gap.extend(pos, vecM);
		} else { //this is a new gap
			last_gap.setStart(ii+1, vecM);
			last_gap.setEnd(pos, vecM);
		}
		return IOUtils::npos;
	} else { //what's left: the current point is in a known gap, but there might be some data before
		size_t ii = last_gap.start; //start from the begining of the last known gap (and idx will be decremented right away)
		for (; ii-- >0; ) {
			if (vecM[ii].date < windowStart) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) { //there is some data before the gap!
				last_gap.extend(ii+1, vecM);
				return ii;
			}
		}
		
		last_gap.extend(ii+1, vecM);
		return IOUtils::npos;
	}
}

//return true if a valid point could be found forward from pos
size_t searchForward(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size, const size_t& indexP1)
{
	const Date windowEnd = (indexP1 != IOUtils::npos)? vecM[indexP1].date+i_window_size : resampling_date+i_window_size;
	const bool knownGap = (!last_gap.startDate.isUndef() && !last_gap.endDate.isUndef());
	const bool currentInGap = (knownGap)? (resampling_date>=last_gap.startDate && resampling_date<=last_gap.endDate) : false;
	const bool windowEndInGap = (knownGap)? (windowEnd>=last_gap.startDate && windowEnd<=last_gap.endDate) : false;
	
	//the current point and window start are in a known gap, there is no hope
	if (currentInGap && windowEndInGap) return IOUtils::npos;
	
	//the current point is NOT in a known gap
	if (!currentInGap) { //or !knownGap
		const Date dateEnd = (windowEndInGap)? last_gap.startDate : windowEnd;
		size_t ii = pos;
		for (; ii<vecM.size(); ++ii) {
			if (vecM[ii].date > dateEnd) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) {
				last_gap.setStart(pos, vecM);
				last_gap.setEnd(ii-1, vecM);
				return ii;
			}
		}
		
		if (windowEndInGap) { //we can extend the current known gap
			last_gap.extend(pos, vecM);
		} else { //this is a new gap
			last_gap.setStart(pos, vecM);
			last_gap.setEnd(ii-1, vecM);
		}
		return IOUtils::npos;
	} else { //what's left: the current point is in a known gap, but there might be some data after
		size_t ii = last_gap.end;
		for (; ii<vecM.size(); ++ii) { //start from the end of the last known gap
			if (vecM[ii].date > windowEnd) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) { //there is some data after the gap!
				last_gap.extend(ii-1, vecM);
				return ii;
			}
		}
		
		last_gap.extend(ii-1, vecM);
		return IOUtils::npos;
	}
}

bool requal(Date &date1, Date &date2) {
    double tolerance = 1e-6; // Define your tolerance level
    return std::abs((date1 - date2).getJulian(true)) <= tolerance;
}


void computeGap(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
                                              const double& i_window_size, size_t& indexP1, size_t& indexP2)
{
	indexP1 = searchBackward(last_gap, pos, paramindex, vecM, resampling_date, i_window_size);
	indexP2 = searchForward(last_gap, pos, paramindex, vecM, resampling_date, i_window_size, indexP1);
    if (indexP1 == IOUtils::npos || indexP2 == IOUtils::npos) {
        std::cerr << "Could not find the end of the gap" << std::endl;
        std::cerr << "Gap start: " << indexP1 << std::endl;
        std::cerr << "Gap end: " << indexP2 << std::endl;
    }
}

// returns the most often accuring value in a vector
double mostLikelyValue(const std::vector<double>& vec) {
    if (vec.empty()) {
        throw mio::IOException("Vector of sampling rates is empty");
    }
    std::map<double, int> counts;
    for (double num : vec) {
        counts[num]++;
    }

    return std::max_element(counts.begin(), counts.end(), [](const std::pair<double,int>& pair1, const std::pair<double,int>& pair2) {
        return pair1.second < pair2.second;
    })->first;
}

// compute the most often occuring sampling rate rounded to 1e-6
double computeSamplingRate(Date data_start_date, Date data_end_date, std::vector<MeteoData> vecM) {
    std::vector<double> time_diffs;
    for (size_t i = 0; i < vecM.size(); i++) {
        if (vecM[i].date >= data_start_date && vecM[i].date <= data_end_date) {
            double value = 1/((vecM[i].date - vecM[i-1].date).getJulian(true));
            time_diffs.push_back(std::round(value*100000/100000));
        }
    }
    return mostLikelyValue(time_diffs);
}
} // namespace mio
