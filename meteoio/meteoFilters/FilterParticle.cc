#include <meteoio/meteoFilters/FilterParticle.h>
#include <meteoio/tinyexpr.h>

#include <sstream> //for readLineToVec

namespace mio {

FilterParticle::FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name), filter_alg(SIR), resample_alg(SYSTEMATIC), NN(2000), path_resampling(false),
          model_expression(""), obs_model_expression(""), model_x0(IOUtils::nodata),
          rng_model(), rng_prior(),
          be_verbose(true)
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
	seedGeneratorsFromIni(RNGU, RNGV, RNG0);

	//init states and noise:
	Eigen::MatrixXd xx(NN, TT); //particles
	Eigen::VectorXd zz(TT); //observation
	vecMeteoToEigen(ivec, zz, param);

	Eigen::MatrixXd ww(NN, TT); //weights of particles

	xx(0, 0) = model_x0;
	ww(0, 0) = 1. / NN;

	for (size_t nn = 1; nn < NN; ++nn) { //draw from prior pdf for initial state of particles at T=0
		xx(nn, 0) = xx(0, 0) + RNG0.doub();
		ww(nn, 0) = 1. / NN; //starting up, all particles have the same weight
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
	    	resample_path(xx, ww, kk);  // This manipulates the original arrays by reference

	} //endfor kk


	Eigen::VectorXd xx_mean = (xx.array() * ww.array()).colwise().sum(); //average over NN weighted particles at each time step

	ovec = ivec; //copy with all special parameters etc.


	for (size_t kk = 0; kk < TT; ++kk) {
		te_xx = xx_mean(kk);
		ovec[kk](param) = te_eval(expr_obs); //filtered observation (model function of mean state [= estimated likely state])
	}

	te_free(expr_model);
	te_free(expr_obs);

}

void FilterParticle::resample_path(Eigen::MatrixXd& xx, Eigen::MatrixXd& ww, const int& kk)
{ //if a lot of computational power is devoted to particles with low contribution (low weight), resample the paths
	switch(resample_alg) //choose resampling algorithm
	{
	case SYSTEMATIC:

		double N_eff = 0.;
		for (int nn = 0; nn < NN; ++nn)
			N_eff += ww(nn, kk)*ww(nn, kk); //effective sample size
		N_eff = 1. / N_eff;

		static const double rc = 0.5;

		RandomNumberGenerator RNU; //uniform random numbers

		if (N_eff < rc * NN)
		{

			double cdf[NN];
			cdf[0] = 0.;
			for (int nn = 1; nn < NN; ++nn)
				cdf[nn] = cdf[nn-1] + ww(nn, kk); //construct cumulative density function
			cdf[NN-1] = 1.0; //round-off protection

			double rr = RNU.doub() / NN;

			for (int nn = 0; nn < NN; ++nn) //for each PARTICLE...
			{
				size_t jj = 0;
				while (rr > cdf[jj])
					++jj; //check which range in the cdf the random number belongs to...

				xx(nn, kk) = xx(jj, kk); //... and use that index

				ww(nn, kk) = 1. / NN; //all resampled particles have the same weight
				rr += 1. / NN; //move along cdf
			}

		} //endif N_eff
		break;

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
		}

		if (!has_prior) //not one prior pdf RNG setting -> pick importance density
			rng_prior = rng_model;
		if (!has_obs) //no dedicated observation noise pdf - use prior for importance sampling
			rng_obs = rng_prior;

		if (model_expression.empty())
			throw NoDataException("No model function supplied for the particle filter.", AT);

	} //endfor ii
}


void FilterParticle::seedGeneratorsFromIni(RandomNumberGenerator& RNGU, RandomNumberGenerator& RNGV, RandomNumberGenerator& RNG0)
{ //to keep process(...) more readable
	if (!rng_model.seed.empty()) //seed integers given in ini file
		RNGU.setState(rng_model.seed);
	if (!rng_obs.seed.empty())
		RNGV.setState(rng_obs.seed);
	if (!rng_prior.seed.empty())
		RNG0.setState(rng_prior.seed);
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
