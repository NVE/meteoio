#ifndef FILTERKALMAN_H
#define FILTERKALMAN_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <Dense>

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
		Eigen::MatrixXd buildControlSignal(const size_t& nx, const size_t& TT, const std::vector<MeteoData>& ivec);
		Eigen::MatrixXd parseMatrix(const std::string& line, const size_t& rows, const size_t& cols);
		Eigen::MatrixXd bloatMatrix(const std::string& line, const size_t& rows, const size_t& cols);
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);

		std::string mat_in_xx, mat_in_AA, mat_in_HH, mat_in_PP, mat_in_QQ, mat_in_RR, mat_in_BB, mat_in_uu;
		std::vector<std::string> meas_params;

		bool be_verbose; //output warnings/info?
};

} //namespace

#endif //FILTERKALMAN_H
