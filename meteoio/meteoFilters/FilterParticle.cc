#include <meteoio/meteoFilters/FilterParticle.h>

#include <cmath> //for isnan()
#include <limits>
#include <fstream> //for dump files
#include <sstream> //for readLineToVec

namespace mio {

FilterParticle::FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name), filter_alg(SIR), resample_alg(SYSTEMATIC), NN(2000), path_resampling(false),
          model_expression(""), obs_model_expression(""), model_x0(IOUtils::nodata), resample_percentile(0.5),
          rng_model(), rng_obs(), rng_prior(), resample_seed(),
          be_verbose(true), unrecognized_key(""), dump_particles_file(""), dump_states_file(""), input_states_file("")
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}

void FilterParticle::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{

	/* INITIALIZATION */

	if (!unrecognized_key.empty() && be_verbose)
		std::cerr << "[W] Unrecognized ini key(s) ignored for particle filter, one of them is: \"" + unrecognized_key + "\"\n";
	const size_t TT = ivec.size(); //number of time steps
	ovec = ivec; //copy with all special parameters etc.

	const bool nodata_check = checkInitialState(ivec, param); //if not provided, find first valuable meteo data point for x0
	if (!nodata_check)
		return; //nothing to do, keep nodata values

	//init random number generators:
	RandomNumberGenerator RNGU(rng_model.algorithm, rng_model.distribution, rng_model.parameters); //process noise
	RandomNumberGenerator RNGV(rng_obs.algorithm, rng_obs.distribution, rng_obs.parameters); //observation pdf
	RandomNumberGenerator RNG0(rng_prior.algorithm, rng_prior.distribution, rng_prior.parameters); //prior pdf
	RandomNumberGenerator RNU; //uniforms for resampling
	seedGeneratorsFromIni(RNGU, RNGV, RNG0, RNU);

	//init states:
	Eigen::MatrixXd xx(NN, TT); //particles
	Eigen::VectorXd zz(TT); //observations
	vecMeteoToEigen(ivec, zz, param);
	Eigen::MatrixXd ww(NN, TT); //weights of particles
	const Eigen::VectorXd tVec = buildTimeVector(ivec);

	bool instates_success(false);
	if (!input_states_file.empty()) //there is data saved from a previous run
		instates_success = readInternalStates(xx, ww);  //online data aggregation
	if (!instates_success) { //start from the initial value
		xx(0, 0) = model_x0;
		ww(0, 0) = 1. / NN;
		for (size_t nn = 1; nn < NN; ++nn) { //draw from prior pdf for initial state of particles at T=0
			xx(nn, 0) = xx(0, 0) + RNG0.doub();
			ww(nn, 0) = 1. / NN; //starting up, all particles have the same weight
		}
	}

	//prepare system model and observation model expressions
	std::vector<std::string> sub_expr, sub_params;
	parseBracketExpression(model_expression, sub_expr, sub_params); //get substitution strings and index map for the meteo parameters
	std::vector<double> sub_values(sub_expr.size()); //empty so far but with reserved memory to point to

	te_variable *te_vars = new te_variable[sub_expr.size()];
	initFunctionVars(te_vars, sub_expr, sub_values); //build te_variables from substitution vectors

	te_expr *expr_model = compileExpression(model_expression, te_vars, sub_expr.size()); //ready the lazy equation
	te_expr *expr_obs = compileExpression(obs_model_expression, te_vars, sub_expr.size()); //(with syntax check)

	/*
	 * SUBSTITUTIONS:
	 *     [0]: index (k)
	 *     [1]: time index (t)
	 *     [2]: value (x_k) - only in observations equation as x_mean
	 *     [3]: previous value (x_km1) - arbitrarily x0 for k=0 in observations equation
	 */
	static const size_t nr_hardcoded_sub = 4; //in order not to forget to keep this in sync

	/* PARTICLE FILTER */

	bool saw_nodata(false);
	try {
		for (size_t kk = 1; kk < TT; ++kk) { //for each TIME STEP (starting at 2nd)...

			sub_values[0] = (double)kk; //this vector is linked to tinyexpr expression memory
			sub_values[1] = tVec(kk);

			for (size_t jj = 0; jj < sub_params.size(); ++jj) //fill current meteo parameters
				sub_values[jj+nr_hardcoded_sub] = ivec[kk](IOUtils::strToUpper(sub_params[jj]));

			if (ivec[kk](param) == IOUtils::nodata) {
				xx.col(kk) = xx.col(kk-1); //repeat particles for nodata values
				ww.col(kk) = ww.col(kk-1); //the mean gets skewed a little
				saw_nodata = true;
			} else {
				//SIR algorithm
				for (size_t nn = 0; nn < NN; ++nn) { //for each PARTICLE...
					sub_values[3] = xx(nn, kk-1);
					sub_values[2] = sub_values[3]; //we don't know x but put something for protection (meant for obs model)
					const double res = te_eval(expr_model); //evaluate expression with current substitution values
					xx(nn, kk) = res + RNGU.doub(); //generate system noise
					ww(nn, kk) = ww(nn, kk-1) * RNGV.pdf( zz(kk) - xx(nn, kk) );
				} //endfor nn
			} //endif nodata

			ww.col(kk) /= ww.col(kk).sum(); //normalize weights to sum=1 per timestep

			if (path_resampling)
				resamplePaths(xx, ww, kk, RNU);

		} //endfor kk
	} catch (...) { //we could get a "nonexistent meteo parameter" error above
		te_free(expr_model);
		te_free(expr_obs);
		delete[] te_vars;
		throw;
	}

	Eigen::VectorXd xx_mean = (xx.array() * ww.array()).colwise().sum(); //average over NN weighted particles at each time step
	for (size_t kk = 0; kk < TT; ++kk) {
		sub_values[0] = (double)kk;
		sub_values[1] = tVec(kk);
		sub_values[2] = xx_mean(kk);
		sub_values[3] = (kk == 0)? model_x0 : xx_mean(kk-1); //somewhat arbitrary - this substitution is meant for the system model
		const double res = te_eval(expr_obs); //filtered observation (model function of mean state [= estimated likely state])
		ovec[kk](param) = isnan(res)? IOUtils::nodata : res; //NaN to nodata
	}

	te_free(expr_model);
	te_free(expr_obs);
	delete[] te_vars;

	if (be_verbose && saw_nodata) std::cerr << "[W] Nodata value(s) encountered in particle filter. For this, the previous particle was repeated. You should probably resample beforehand.\n";

	if (!dump_states_file.empty())
		(void) dumpInternalStates(xx, ww);
	if (!dump_particles_file.empty())
		(void) dumpParticlePaths(xx);

}

void FilterParticle::resamplePaths(Eigen::MatrixXd& xx, Eigen::MatrixXd& ww, const size_t& kk, RandomNumberGenerator& RNU) const
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

bool FilterParticle::checkInitialState(const std::vector<MeteoData>& ivec, const size_t& param)
{ //check for nodata values in the input meteo data
	if (model_x0 == IOUtils::nodata) {
		for (size_t kk = 0; kk < ivec.size(); ++kk) { //find 1st data element for starting point
			if (ivec[kk](param) != IOUtils::nodata) {
				model_x0 = ivec[kk](param); //if there are many the filter will run through, but the user should really resample
				break;
			}
		}
		if (be_verbose) std::cerr << "[W] No initial state value x_0 provided for particle filter; using 1st available measurement.\n";
		if (model_x0 == IOUtils::nodata) //all values are nodata
			return false;
	} //endif nodata
	return true;
}

void FilterParticle::initFunctionVars(te_variable* vars, const std::vector<std::string>& expr, const std::vector<double>& values)
{ //build a substitutions expression for tinyexpr
	if (values.size() != expr.size())
		throw InvalidArgumentException("Particle filter: error mapping meteo(param) fields to meteo values. Are all fields available?\n", AT);

	for (size_t ii = 0; ii < expr.size(); ++ii) {
		vars[ii].name = expr[ii].c_str();
		vars[ii].address = &values.data()[ii];
		vars[ii].type = 0;
		vars[ii].context = 0;
	}
}

te_expr* FilterParticle::compileExpression(const std::string& expression, const te_variable* te_vars, const size_t& sz) const
{
	int te_err;
	te_expr *expr = te_compile(expression.c_str(), te_vars, (int)sz, &te_err); //ready the lazy equation with variables including syntax check
	if (!expr)
		throw InvalidFormatException("Arithmetic expression \"" + model_expression +
		        "\" could not be evaluated for particle filter; parse error at " + IOUtils::toString(te_err), AT);
	return expr;
}

void FilterParticle::parseBracketExpression(std::string& line, std::vector<std::string>& sub_expr,
        std::vector<std::string>& sub_params) const
{ //list names of meto parameters given by meteo(...) in input matrix
	sub_expr.clear();
	sub_params.clear();
	sub_expr.push_back("kk"); //hard-coded parameters
	sub_expr.push_back("tt");
	sub_expr.push_back("xx");
	sub_expr.push_back("x_km1");

	const std::string prefix("meteo");
	const size_t len = prefix.length() + 1;
	size_t pos1 = 0, pos2;
	while (true) {
		pos1 = line.find(prefix, pos1);
		if (pos1 == std::string::npos)
			break; //done
		pos2 = line.find(")", pos1+len);
		if (pos2 == std::string::npos || pos2-pos1-len == 0) { //no closing bracket
			throw InvalidArgumentException("Missing closing bracket in meteo(...) part of particle filter's system model.", AT);
		}

		const std::string pname = IOUtils::strToLower(line.substr(pos1+len, pos2-pos1-len));
		sub_params.push_back(pname); //meteo name
		line.replace(pos1, pos2-pos1+1, prefix + pname + "  "); //to make parsable with tinyexpr: 'meteo(RH)' --> 'meteorh  '
		sub_expr.push_back(prefix + pname); //full expression, lower case and without brackets
		pos1 += len;
	}
}

bool FilterParticle::dumpInternalStates(Eigen::MatrixXd& particles, Eigen::MatrixXd& weights) const
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
	static const int digits = std::numeric_limits<double>::digits10;
	oss.precision(digits);
	oss.setf(std::ios::fixed);
	for (int ii = 0; ii < particles.rows(); ++ii) {
		oss << std::setw(digits) << particles.rightCols(1)(ii) << "   ";
		oss << std::setw(digits) << weights.rightCols(1)(ii) << std::endl;
	}
	oss.close();
	return true;
}

bool FilterParticle::readInternalStates(Eigen::MatrixXd& particles, Eigen::MatrixXd& weights) const
{
	std::vector<double> xx, ww;
	try {
		readCorrections(block_name, input_states_file, xx, ww); //file not available yet - will become active on next run
	} catch (...) {
		return false;
	}

	if ( (unsigned int)particles.rows() != xx.size() || (unsigned int)weights.rows() != ww.size() ) {
		if (be_verbose) std::cerr << "[W] Particle filter file input via INPUT_STATES_FILE does not match the number of particles. Using MODEL_X0.\n";
		return false;
	}

	Eigen::Map<Eigen::VectorXd> tmp_p(&xx.data()[0], xx.size()); //memcopy read vectors to Eigen types
	particles.col(0) = tmp_p;
	Eigen::Map<Eigen::VectorXd> tmp_w(&ww.data()[0], ww.size());
	weights.col(0) = tmp_w;
	return true;
}

bool FilterParticle::dumpParticlePaths(Eigen::MatrixXd& particles) const
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

Eigen::VectorXd FilterParticle::buildTimeVector(const std::vector<MeteoData>& ivec) const
{ //time representation of the index: shift and normalize to t(0)=0, t(1)=1
	double base_time = ivec.front().date.getJulian();
	double dt(1.);

	if (ivec.size() > 1)
		dt = ivec[1].date.getJulian() - base_time;

	Eigen::VectorXd vecRet(ivec.size());
	for (size_t ii = 0; ii < ivec.size(); ++ii) {
		const double tt = (ivec[ii].date.getJulian() - base_time) / dt;
		vecRet(ii) = tt;
	}

	return vecRet;
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
		} else if (vecArgs[ii].first == "INITIAL_STATE") {
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
		} else {
			unrecognized_key = vecArgs[ii].first; //constructor is always called twice - show only when processing
		}

	} //endfor vecArgs

	if (model_expression.empty())
		throw NoDataException("No model function supplied for the particle filter.", AT);

	if (obs_model_expression.empty()) {
		if (be_verbose) std::cerr << "[W] Model to relate observations to state is missing for particle filter. Picking obs(k)=state(k).\n";
		obs_model_expression = "xx";
	}
	if (!has_prior) //not one prior pdf RNG setting -> pick importance density
		rng_prior = rng_model;
	if (!has_obs) //no dedicated observation noise pdf - use prior for importance sampling
		rng_obs = rng_prior;
}

void FilterParticle::seedGeneratorsFromIni(RandomNumberGenerator& RNGU, RandomNumberGenerator& RNGV, RandomNumberGenerator& RNG0,
        RandomNumberGenerator& RNU) const
{ //to keep process(...) less cluttered
	if (!rng_model.seed.empty()) //seed integers given in ini file
		RNGU.setState(rng_model.seed);
	if (!rng_obs.seed.empty())
		RNGV.setState(rng_obs.seed);
	if (!rng_prior.seed.empty())
		RNG0.setState(rng_prior.seed);
	if (!resample_seed.empty())
		RNU.setState(resample_seed);
}

void FilterParticle::vecMeteoToEigen(const std::vector<MeteoData>& vec, Eigen::VectorXd& eig, const unsigned int& param) const
{ //naive copy of 1 meteo parameter in STL vector to Eigen library vector
	eig.resize(vec.size());
	for (size_t i = 0; i < vec.size(); ++i)
		eig[i] = vec[i](param);
}

void FilterParticle::readLineToVec(const std::string& line_in, std::vector<uint64_t>& vec_out) const
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
