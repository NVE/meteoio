#ifndef FILTERPARTICLE_H
#define FILTERPARTICLE_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <vector>
#include <string>

namespace mio {

class FilterParticle : public ProcessingBlock {
	public:
		FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		bool fill_state(const unsigned int& param, const std::vector<MeteoData>& ivec);

		enum FILTER_ALGORITHM{SIS, SIR}; //only SIR is implemented so far
		enum RESAMPLE_ALGORITHM{SYSTEMATIC, EPANECHNIKOV}; //only SYSTEMATIC implemented so far

		FILTER_ALGORITHM filter_alg;
		RESAMPLE_ALGORITHM resample_alg; //resampling of particle paths

		unsigned int NN; //number of particles
		bool path_resampling; //has nothing to do with temporal or spatial meteo resampling

		std::string model_expression; //model formula (as opposed to file input)
		double model_x0; //initial state at T=0

		std::vector<double> xx; //state / model
		std::vector<double> zz; //observation

		RandomNumberGenerator::RNG_TYPE rng_alg; //all RNGs will use this algorithm

		bool be_verbose; //output warnings/info?

};

} //namespace

#endif //FILTERPARTICLE_H
