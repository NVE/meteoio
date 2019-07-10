#include <meteoio/meteoFilters/FilterParticle.h>
#include <meteoio/tinyexpr.h>

#include <limits>
#include <fstream> //for dump files
#include <sstream> //for readLineToVec

namespace mio {

FilterParticle::FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name), filter_alg(SIR), resample_alg(SYSTEMATIC), NN(2000), path_resampling(false),
          model_expression(""), obs_model_expression(""), model_x0(IOUtils::nodata), resample_percentile(0.5),
          rng_model(), rng_obs(), rng_prior(), resample_seed(),
          be_verbose(true), dump_particles_file(""), dump_states_file(""), input_states_file("")
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}


void FilterParticle::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{

	/* INITIALIZATION */

	const size_t TT = ivec.size(); //number of time steps

	if (obs_model_expression.empty()) {
		std::cerr << "[W] Model to relate observations to state is missing for particle filter. Picking obs(x)=state(x).\n";
		obs_model_expression = "xx";
	}

	//prepare system model expression
	double te_kk, te_x_km1; //2 substitutions available: index and state value at previous time step
	te_variable te_vars[] = {{"kk", &te_kk, 0, 0}, {"x_km1", &te_x_km1, 0, 0}}; //read: "x_(k-1)"

	int te_err;
	te_expr *expr_model = te_compile(model_expression.c_str(), te_vars, 2, &te_err); //ready the lazy equation with variables including syntax check
	if (!expr_model)
		throw InvalidFormatException("Arithmetic expression '" + model_expression +
		        "'could not be evaluated for particle filter; parse error at " + IOUtils::toString(te_err), AT);
	if (model_x0 == IOUtils::nodata) {
		if (be_verbose) std::cerr << "[W] No initial state value x_0 provided for particle filter; using 1st measurement.\n";
		model_x0 = ivec.front()(param); //TODO: nodata checks
	}

	//prepare observation model expression
	double te_xx;
	te_variable te_vars_obs[] = {{"xx", &te_xx, 0, 0}};
	int te_err_obs;
	te_expr *expr_obs = te_compile(obs_model_expression.c_str(), te_vars_obs, 1, &te_err_obs);
	if (!expr_obs)
		throw InvalidFormatException("Arithmetic expression '" + obs_model_expression +
		        "'could not be evaluated for particle filter; parse error at " + IOUtils::toString(te_err_obs), AT);

	//init random number generators:
	RandomNumberGenerator RNGU(rng_model.algorithm, rng_model.distribution, rng_model.parameters); //process noise
	RandomNumberGenerator RNGV(rng_obs.algorithm, rng_obs.distribution, rng_obs.parameters); //observation pdf
	RandomNumberGenerator RNG0(rng_prior.algorithm, rng_prior.distribution, rng_prior.parameters); //prior pdf
	RandomNumberGenerator RNU; //uniforms for resampling
	seedGeneratorsFromIni(RNGU, RNGV, RNG0, RNU);

	//init states and noise:
	Eigen::MatrixXd xx(NN, TT); //particles
	Eigen::VectorXd zz(TT); //observations
	vecMeteoToEigen(ivec, zz, param);
	Eigen::MatrixXd ww(NN, TT); //weights of particles

	bool instates_success(false);
	if (!input_states_file.empty()) { //there is data saved from a previous run
		instates_success = readInternalStates(xx, ww);  //online data aggregation
	}
	if (!instates_success) { //start from the initial value
		xx(0, 0) = model_x0;
		ww(0, 0) = 1. / NN;
		for (size_t nn = 1; nn < NN; ++nn) { //draw from prior pdf for initial state of particles at T=0
			xx(nn, 0) = xx(0, 0) + RNG0.doub();
			ww(nn, 0) = 1. / NN; //starting up, all particles have the same weight
		}
	}

	/* PARTICLE FILTER */

	for (size_t kk = 1; kk < TT; ++kk) { //for each TIME STEP (starting at 2nd)...

		//SIR algorithm
		for (size_t nn = 0; nn < NN; ++nn) { //for each PARTICLE...
			te_kk = (double)kk;
			te_x_km1 = xx(nn, kk-1);

			const double res = te_eval(expr_model); //evaluate expression with current substitution values
			xx(nn, kk) = res + RNGU.doub();
			ww(nn, kk) = ww(nn, kk-1) * RNGV.pdf( zz(kk) - xx(nn, kk) );
		} //endfor nn

		ww.col(kk) /= ww.col(kk).sum(); //normalize weights to sum=1 per timestep

	    if (path_resampling)
	    	resample_path(xx, ww, kk, RNU);

	} //endfor kk

	Eigen::VectorXd xx_mean = (xx.array() * ww.array()).colwise().sum(); //average over NN weighted particles at each time step

	ovec = ivec; //copy with all special parameters etc.
	for (size_t kk = 0; kk < TT; ++kk) {
		te_xx = xx_mean(kk);
		ovec[kk](param) = te_eval(expr_obs); //filtered observation (model function of mean state [= estimated likely state])
	}

	te_free(expr_model);
	te_free(expr_obs);

	if (!dump_states_file.empty())
		(void) dumpInternalStates(xx, ww);
	if (!dump_particles_file.empty())
		(void) dumpParticlePaths(xx);

}

void FilterParticle::resample_path(Eigen::MatrixXd& xx, Eigen::MatrixXd& ww, const size_t& kk, RandomNumberGenerator& RNU)
{ //if a lot of computational power is devoted to particles with low contribution (low weight), resample the paths
	if (resample_alg == SYSTEMATIC) { //choose resampling algorithm

		double N_eff = 0.;
		for (size_t nn = 0; nn < NN; ++nn)
			N_eff += ww(nn, kk)*ww(nn, kk); //effective sample size
		N_eff = 1. / N_eff;

		if (N_eff < resample_percentile * (float)NN)
		{
			std::vector<double> cdf(NN);
			cdf.front() = 0.;
			for (size_t nn = 1; nn < NN; ++nn)
				cdf[nn] = cdf[nn-1] + ww(nn, kk); //construct cumulative density function
			cdf.back() = 1.0; //round-off protection

			double rr = RNU.doub() / NN;

			for (size_t nn = 0; nn < NN; ++nn) //for each PARTICLE...
			{
				size_t jj = 0;
				while (rr > cdf[jj])
					++jj; //check which range in the cdf the random number belongs to...

				xx(nn, kk) = xx(jj, kk); //... and use that index

				ww(nn, kk) = 1. / NN; //all resampled particles have the same weight
				rr += 1. / NN; //move along cdf
			}

		} //endif N_eff

	} else {
		throw InvalidArgumentException("Resampling strategy for particle filter not implemented.", AT);
	} //end switch resampling
}

void FilterParticle::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where("Filters::" + block_name);

	bool has_prior(false), has_obs(false);

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {

		/*** FILTER settings ***/
		if (vecArgs[ii].first == "NO_OF_PARTICLES") {
			IOUtils::parseArg(vecArgs[ii], where, NN);
		} else if (vecArgs[ii].first == "PATH_RESAMPLING") {
			IOUtils::parseArg(vecArgs[ii], where, path_resampling);
		}

		/*** MODEL FUNCTION settings ***/
		else if (vecArgs[ii].first == "MODEL_FUNCTION") {
			model_expression = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "MODEL_X0") {
			IOUtils::parseArg(vecArgs[ii], where, model_x0);
		} else if (vecArgs[ii].first == "OBS_MODEL_FUNCTION") {
			obs_model_expression = vecArgs[ii].second;
		}

		/*** MODEL RNG settings ***/
		else if (vecArgs[ii].first == "MODEL_RNG_ALGORITHM") { //everything RNG will be defaulted if not provided - cf. RNG doc
			rng_model.algorithm = RandomNumberGenerator::strToRngtype(vecArgs[ii].second);
		} else if (vecArgs[ii].first == "MODEL_RNG_DISTRIBUTION") {
			rng_model.distribution = RandomNumberGenerator::strToRngdistr(vecArgs[ii].second); //convert from int to enum RNG_DISTR
		} else if (vecArgs[ii].first == "MODEL_RNG_PARAMETERS") {
			IOUtils::readLineToVec(vecArgs[ii].second, rng_model.parameters);
		} else if (vecArgs[ii].first == "MODEL_RNG_SEED") {
			readLineToVec(vecArgs[ii].second, rng_model.seed);
		}

		/*** PRIOR PDF RNG settings ***/
		else if (vecArgs[ii].first == "PRIOR_RNG_ALGORITHM") {
			rng_prior.algorithm = RandomNumberGenerator::strToRngtype(vecArgs[ii].second);
			has_prior = true;
		} else if (vecArgs[ii].first == "PRIOR_RNG_DISTRIBUTION") {
			rng_prior.distribution = RandomNumberGenerator::strToRngdistr(vecArgs[ii].second); //convert from int to enum RNG_DISTR
			has_prior = true;
		} else if (vecArgs[ii].first == "PRIOR_RNG_PARAMETERS") {
			IOUtils::readLineToVec(vecArgs[ii].second, rng_prior.parameters);
			has_prior = true;
		} else if (vecArgs[ii].first == "PRIOR_RNG_SEED") {
			readLineToVec(vecArgs[ii].second, rng_prior.seed);
			has_prior = true;
		}

		/*** OBSERVATION PDF RNG settings ***/
		else if (vecArgs[ii].first == "OBS_RNG_ALGORITHM") {
			rng_obs.algorithm = RandomNumberGenerator::strToRngtype(vecArgs[ii].second);
			has_obs = true;
		} else if (vecArgs[ii].first == "OBS_RNG_DISTRIBUTION") {
			rng_obs.distribution = RandomNumberGenerator::strToRngdistr(vecArgs[ii].second); //convert from int to enum RNG_DISTR
			has_obs = true;
		} else if (vecArgs[ii].first == "OBS_RNG_PARAMETERS") {
			IOUtils::readLineToVec(vecArgs[ii].second, rng_obs.parameters);
			has_obs = true;
		} else if (vecArgs[ii].first == "OBS_RNG_SEED") {
			readLineToVec(vecArgs[ii].second, rng_obs.seed);
			has_obs = true;
		}

		/*** MISC settings ***/
		else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);
		} else if (vecArgs[ii].first == "DUMP_PARTICLES_FILE") {
			dump_particles_file = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "DUMP_INTERNAL_STATES_FILE") {
			dump_states_file = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "INPUT_INTERNAL_STATES_FILE") {
			input_states_file = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "RESAMPLE_PERCENTILE") {
			IOUtils::parseArg(vecArgs[ii], where, resample_percentile);
		} else if (vecArgs[ii].first == "RESAMPLE_RNG_SEED") {
			readLineToVec(vecArgs[ii].second, resample_seed);
		}

	} //endfor vecArgs

	if (!has_prior) //not one prior pdf RNG setting -> pick importance density
		rng_prior = rng_model;
	if (!has_obs) //no dedicated observation noise pdf - use prior for importance sampling
		rng_obs = rng_prior;

	if (model_expression.empty())
		throw NoDataException("No model function supplied for the particle filter.", AT);
}

bool FilterParticle::dumpInternalStates(Eigen::MatrixXd& particles, Eigen::MatrixXd& weights)
{ //using this, we are able to resume our filter without having to recalculate the past if new data arrives
	std::ofstream oss(dump_states_file.c_str(), std::ofstream::out);
	if (oss.fail()) {
		std::ostringstream ss;
		ss << "Particle filter could not dump internal states to \"" << dump_states_file;
		ss << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	oss << "# This file was generated by MeteoIO's particle filter and holds the particles and weights at the last time step." << std::endl;
	oss << "# [particles]   [weights]" << std::endl;
	const int digits = std::numeric_limits<double>::digits10;
	oss.precision(digits);
	oss.setf(std::ios::fixed);
	for (int ii = 0; ii < particles.rows(); ++ii) {
		oss << std::setw(digits) << particles.rightCols(1)(ii) << "   ";
		oss << std::setw(digits) << weights.rightCols(1)(ii) << std::endl;
	}
	oss.close();
	return true;
}

bool FilterParticle::readInternalStates(Eigen::MatrixXd& particles, Eigen::MatrixXd& weights)
{
	std::vector<double> xx, ww;
	readCorrections(block_name, input_states_file, xx, ww);
	if ( (unsigned int)particles.rows() != xx.size() || (unsigned int)weights.rows() != ww.size() ) {
		if (be_verbose) std::cerr << "[W] Particle filter file input via INPUT_STATES_FILE does not match the number of particles. Using MODEL_X0.\n";
		return false;
	}

	Eigen::Map<Eigen::VectorXd> tmp_p(&xx.data()[0], xx.size()); //memcopy read vectors to Eigen types
	particles.col(0) = tmp_p;
	Eigen::Map<Eigen::VectorXd> tmp_w(&xx.data()[0], ww.size());
	weights.col(0) = tmp_w;
	return true;
}

bool FilterParticle::dumpParticlePaths(Eigen::MatrixXd& particles)
{ //to plot paths and kernel density outside of MeteoIO
	std::ofstream oss(dump_particles_file.c_str(), std::ofstream::out);
	if (oss.fail()) {
		std::ostringstream ss;
		ss << "Particle filter could not dump particle paths states to \"" << dump_states_file;
		ss << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	oss << "# This file was generated by MeteoIO's particle filter and holds the paths of all particles." << std::endl;
	oss << "# Rows are the particles, columns the time steps." << std::endl;
	const int digits = std::numeric_limits<double>::digits10;
	oss.precision(digits);
	oss.setf(std::ios::fixed);
	oss << particles;
	oss.close();
	return true;
}

void FilterParticle::seedGeneratorsFromIni(RandomNumberGenerator& RNGU, RandomNumberGenerator& RNGV, RandomNumberGenerator& RNG0,
        RandomNumberGenerator& RNU)
{ //to keep process(...) more readable
	if (!rng_model.seed.empty()) //seed integers given in ini file
		RNGU.setState(rng_model.seed);
	if (!rng_obs.seed.empty())
		RNGV.setState(rng_obs.seed);
	if (!rng_prior.seed.empty())
		RNG0.setState(rng_prior.seed);
	if (!resample_seed.empty())
		RNU.setState(resample_seed);
}

void FilterParticle::vecMeteoToEigen(const std::vector<MeteoData>& vec, Eigen::VectorXd& eig, const unsigned int& param)
{ //naive copy of 1 meteo parameter in STL vector to Eigen library vector
	eig.resize(vec.size());
	for (size_t i = 0; i < vec.size(); ++i)
		eig[i] = vec[i](param);
}

void FilterParticle::readLineToVec(const std::string& line_in, std::vector<uint64_t>& vec_out)
{ //uint64 type version
	vec_out.clear();
	std::istringstream iss(line_in);
	uint64_t val;
	while (!iss.eof()) {
		iss >> std::skipws >> val;
		if (iss.fail())
			throw InvalidFormatException("Unable to parse process noise seed integers for particle filter.", AT);
		vec_out.push_back(val);
	}
}

} //namespace
