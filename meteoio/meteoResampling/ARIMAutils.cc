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
std::vector<double> decideDirection(std::vector<double> data, std::string direction, bool forward, size_t gap_loc) {
    if (forward) {
        if (direction == "forward") {
            return slice(data,0,gap_loc-1);
        } else if (direction == "backward") {
            reverseVector(data);
            return slice(data,0,gap_loc-1);
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
				last_gap.setStart(ii, vecM);
				last_gap.setEnd(pos, vecM);
				return ii;
			}
		}
		
		//no valid point found
		if (windowStartInGap) { //we can extend the current known gap
			last_gap.extend(pos, vecM);
		} else { //this is a new gap
			last_gap.setStart(ii, vecM);
			last_gap.setEnd(pos, vecM);
		}
		return IOUtils::npos;
	} else { //what's left: the current point is in a known gap, but there might be some data before
		size_t ii = last_gap.start; //start from the begining of the last known gap (and idx will be decremented right away)
		for (; ii-- >0; ) {
			if (vecM[ii].date < windowStart) break;
			if (vecM[ii](paramindex) != IOUtils::nodata) { //there is some data before the gap!
				last_gap.extend(ii, vecM);
				return ii;
			}
		}
		
		last_gap.extend(ii, vecM);
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

		if (ii == pos) {
			if (vecM[ii-1].date <= dateEnd && vecM[ii].date >dateEnd) {
				if (windowEndInGap) { //we can extend the current known gap
					last_gap.extend(pos-1, vecM);
				} else { //this is a new gap
					last_gap.setStart(pos-1, vecM);
					last_gap.setEnd(ii, vecM);
				}				
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
    bool is_equal = std::abs((date1 - date2).getJulian(true)) <= tolerance;
	return is_equal;
}


void computeARIMAGap(ARIMA_GAP &last_gap, const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date, size_t& indexP1, size_t& indexP2, double& before_window, double& after_window, double& window_size, Date& data_start_date, Date& data_end_date)
{
	indexP1 = searchBackward(last_gap, pos, paramindex, vecM, resampling_date, window_size);
	indexP2 = searchForward(last_gap, pos, paramindex, vecM, resampling_date, window_size, indexP1);
    if (indexP1 == IOUtils::npos || indexP2 == IOUtils::npos) {
        std::cerr << "Could not find the end of the gap "<< last_gap.toString() << std::endl;
        std::cerr << "Gap start: " << indexP1 << std::endl;
        std::cerr << "Gap end: " << indexP2 << std::endl;
		std::cerr << "Consider increasing the Window Size" << std::endl;
    }
	data_start_date = last_gap.startDate - before_window;
	data_end_date = last_gap.endDate + after_window;	

	// check that window size is not overreached
	if (data_start_date < vecM[0].date) {
		data_start_date = vecM[0].date;
	}
	if (data_end_date > vecM[vecM.size()-1].date) {
		data_end_date = vecM[vecM.size()-1].date;
	}
	if (data_end_date - data_start_date > window_size) {
		throw IOException("The data window needed to interpolate the gap " + last_gap.toString() + " is larger than the resampling window size");
	}
	last_gap.sampling_rate = computeSamplingRate(data_start_date, data_end_date, vecM);
	
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

void printVectors(const std::vector<MeteoData>& vecM, const std::vector<Date>& dates, const size_t& paramindex) {
	size_t maxSize = std::max(vecM.size(), dates.size());

	// Print headers
	std::cout << std::left << std::setw(30) << "MeteoData Date" << "| Date" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;

	for (size_t i = 0; i < maxSize; i++) {
		// Print date from MeteoData or "NaN" if out of range
		if (i < vecM.size()) {
			std::cout << std::left << std::setw(30) << vecM[i].date.toString(Date::ISO)<< ":" << vecM[i](paramindex) << "| ";
		} else {
			std::cout << std::left << std::setw(30) << "NaN" << "| ";
		}

		// Print date from dates vector or "NaN" if out of range
		if (i < dates.size()) {
			std::cout << dates[i].toString(Date::ISO) << std::endl;
		} else {
			std::cout << "NaN" << std::endl;
		}
	}
}

void printVectors(const std::vector<Date>& dates1, const std::vector<Date>& dates2){ 
	size_t maxSize = std::max(dates1.size(), dates2.size());

	// Print headers
	std::cout << std::left << std::setw(30) << "Date1" << "| Date2" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;

	for (size_t i = 0; i < maxSize; i++) {
		// Print date from dates1 vector or "NaN" if out of range
		if (i < dates1.size()) {
			std::cout << std::left << std::setw(30) << dates1[i].toString(Date::ISO) << "| ";
		} else {
			std::cout << std::left << std::setw(30) << "NaN" << "| ";
		}

		// Print date from dates2 vector or "NaN" if out of range
		if (i < dates2.size()) {
			std::cout << dates2[i].toString(Date::ISO) << std::endl;
		} else {
			std::cout << "NaN" << std::endl;
		}
	}
}

void printVector(const std::vector<MeteoData>& vec1) {
	// Print header
	std::cout << std::left << std::setw(30) << "MeteoData Date" << std::endl;
	std::cout << "------------------------------" << std::endl;

	for (size_t i = 0; i < vec1.size(); i++) {
		// Print date from MeteoData
		std::cout << std::left << std::setw(30) << vec1[i].date.toString() << std::endl;
	}
}

void printVectors(const std::vector<double>& vec) {
    // Wrap vec in another vector
    std::vector<std::vector<double>> vecWrapped = {vec};

    // Call the appropriate function
    printVectors(vecWrapped);
}

} // namespace mio
