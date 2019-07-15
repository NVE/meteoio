#ifndef FILTERPARTICLE_H
#define FILTERPARTICLE_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <Core> //<Eigen/Core>
#include <meteoio/tinyexpr.h>

#include <inttypes.h> //for RNG int types
#include <string>
#include <vector>

namespace mio {

class FilterParticle : public ProcessingBlock {
	public:
		FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		        std::vector<MeteoData>& ovec);

	private:
		void resamplePaths(Eigen::MatrixXd& xx, Eigen::MatrixXd& ww, const size_t& kk, RandomNumberGenerator& RNU) const;
		bool checkInitialState(const std::vector<MeteoData>& ivec, const size_t& param);
		void initFunctionVars(te_variable* vars, const std::vector<std::string>& names, const std::vector<double>& meteo);
		te_expr* compileExpression(const std::string& expression, const te_variable* te_vars, const size_t& sz) const;
		void parseBracketExpression(std::string& line, std::vector<std::string>& sub_expr,
		        std::vector<std::string>& sub_params) const;
		bool dumpInternalStates(Eigen::MatrixXd& particles, Eigen::MatrixXd& weights) const;
		bool readInternalStates(Eigen::MatrixXd& particles, Eigen::MatrixXd& weights) const;
		bool dumpParticlePaths(Eigen::MatrixXd& particles) const;
		Eigen::VectorXd buildTimeVector(const std::vector<MeteoData>& ivec) const;
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void seedGeneratorsFromIni(RandomNumberGenerator& RNGU, RandomNumberGenerator& RNGV, RandomNumberGenerator& RNG0,
		        RandomNumberGenerator& RNU) const;
		void vecMeteoToEigen(const std::vector<MeteoData>& vec, Eigen::VectorXd& eig, const unsigned int& param) const;
		void readLineToVec(const std::string& line_in, std::vector<uint64_t>& vec_out) const;

		typedef enum PF_FILTER_ALGORITHM {
			SIR, //only SIR is implemented so far
			SIS
		} pf_filter_algorithm;
		typedef enum PF_RESAMPLE_ALGORITHM {
			SYSTEMATIC, //only SYSTEMATIC implemented so far
			EPANECHNIKOV
		} pf_resample_algorithm;

		pf_filter_algorithm filter_alg;
		pf_resample_algorithm resample_alg; //resampling of particle paths

		unsigned int NN; //number of particles
		bool path_resampling; //has nothing to do with temporal or spatial meteo resampling

		std::string model_expression; //model formula
		std::string obs_model_expression;
		double model_x0; //initial state at T=0

		float resample_percentile;

		struct rng_settings {
			RandomNumberGenerator::RNG_TYPE algorithm;
			RandomNumberGenerator::RNG_DISTR distribution;
			std::vector<double> parameters;
			std::vector<uint64_t> seed;
			rng_settings() : algorithm(RandomNumberGenerator::RNG_XOR),
			                 distribution(RandomNumberGenerator::RNG_GAUSS),
			                 parameters(),
			                 seed() {}
		};
		struct rng_settings rng_model; //noise generator for model function
		struct rng_settings rng_obs; //observation noise pdf
		struct rng_settings rng_prior; //generator to sample from prior pdf
		std::vector<uint64_t> resample_seed;

		bool be_verbose; //output warnings/info?
		std::string unrecognized_key; //to warn about unknown ini keys
		std::string dump_particles_file;
		std::string dump_states_file;
		std::string input_states_file;

};

} //namespace

#endif //FILTERPARTICLE_H
