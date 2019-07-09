#ifndef FILTERKALMAN_H
#define FILTERKALMAN_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <meteoio/Eigen/Dense>

#include <string>
#include <vector>

namespace mio {

class FilterKalman : public ProcessingBlock {
	public:
		FilterKalman(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		        std::vector<MeteoData>& ovec);

	private:
		size_t buildObservationsMatrix(const unsigned int& param, const std::vector<MeteoData>& ivec, Eigen::MatrixXd& zz);
		Eigen::MatrixXd parseMatrix(const std::string& line, const size_t& rows, const size_t& cols);
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);

		std::string matrix_input[3]; //strings from .ini to read matrices from later
		std::vector<std::string> meas_params;

		bool be_verbose; //output warnings/info?
};

} //namespace

#endif //FILTERKALMAN_H
