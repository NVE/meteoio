#ifndef FILTERKALMAN_H
#define FILTERKALMAN_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <Core> //<Eigen/Core>

#include <string>
#include <vector>

namespace mio {

class FilterKalman : public ProcessingBlock {
	public:
		FilterKalman(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		        std::vector<MeteoData>& ovec);

	private:
		Eigen::VectorXd buildInitialStates(const std::vector<std::string>& xx_str, std::vector<size_t>& meas_idx,
		        const std::vector<MeteoData>& ivec, const size_t& nr_observations);
		size_t buildObservationsMatrix(const unsigned int& param, const std::vector<MeteoData>& ivec, Eigen::MatrixXd& zz,
				std::vector<size_t>& meas_idx, const size_t& nx) const;
		Eigen::MatrixXd buildControlSignal(const size_t& nx, const size_t& TT, const std::vector<MeteoData>& ivec) const;
		Eigen::MatrixXd parseMatrix(const std::string& line, const size_t& rows, const size_t& cols,
				const std::string& block) const;
		std::vector<std::string> parseSystemMatrix(const std::string& line, const size_t& rows) const;
		Eigen::MatrixXd buildSystemMatrix(const std::vector<std::string>& AA_str, const size_t& sz, const double& dt,
				const std::vector<MeteoData>& ivec, const size_t& kk) const;
		double substitute(const std::string& expr, const double& dt, const std::vector<MeteoData>& ivec, const size_t& kk) const;
		Eigen::MatrixXd bloatMatrix(const std::string& line, const size_t& rows, const size_t& cols, const std::string& block) const;
		Eigen::VectorXd buildTimeVector(const std::vector<MeteoData>& ivec) const;
		bool checkNodata(const Eigen::VectorXd& ivec) const;
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void cleanBrackets(std::string& iline);

		std::string mat_in_xx, mat_in_AA, mat_in_HH, mat_in_PP, mat_in_QQ, mat_in_RR, mat_in_BB, mat_in_uu;
		std::vector<std::string> meas_params;

		bool be_verbose; //output warnings/info?
		std::string unrecognized_key; //to warn about unknown ini keys
};

} //namespace

#endif //FILTERKALMAN_H
