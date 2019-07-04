#include <meteoio/meteoFilters/FilterParticle.h>
#include <meteoio/tinyexpr.h>

#include <sstream> //for readLineToVec

namespace mio {

static const size_t nx = 1; //dimension of the state vector (how many observables?), fixed so far
static const size_t nu = 1; //dimension of process noise, fixed so far
static const size_t nv = 1; //dimension of observation noise, fixed so far

FilterParticle::FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
          : ProcessingBlock(vecArgs, name), filter_alg(SIR), resample_alg(SYSTEMATIC), NN(2000), path_resampling(false),
            model_appendix("_MOD"), model_expression(""), model_x0(IOUtils::nodata), xx(), zz(),
            model_rng(), add_model_noise(true),
            be_verbose(true)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}

void FilterParticle::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{
	const size_t TT = ivec.size(); //number of time steps

	const bool init_success = fill_state(param, ivec);
	if (!init_success)
    	throw NoDataException("The particle filter could not initialize the state vector. Neither good model input data nor a usable model function are available.", AT);

	ovec = ivec;
	for (size_t ii = 0; ii < TT; ++ii){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		if (tmp < 269) {
			tmp = IOUtils::nodata;
		}
	}
}

bool FilterParticle::fill_state(const unsigned int& param, const std::vector<MeteoData>& ivec)
{
	const size_t TT = ivec.size(); //number of time steps
	xx.clear();

	const bool has_model_expression = !model_expression.empty(); //user has entered a model equation, valid or not

	//check if model data exists in meteo data as per our naming convention:
	const std::string param_name = ivec.front().getNameForParameter(param);
	const size_t param_mod = ivec.front().getParameterIndex(param_name + model_appendix);
	const bool has_param_mod = (param_mod != IOUtils::npos); //user has supplied model data alongside the meteo data

	if (has_model_expression && has_param_mod)
		if (be_verbose) std::cerr << "[W] Both model data and model function are supplied for the particle filter. Ignoring model function.\n";

	if (has_param_mod) {

		for (size_t kk = 1; kk < TT; ++kk) {
			const double val_model = ivec[kk](param_mod);
			if (val_model != IOUtils::nodata && ivec[kk](param) != IOUtils::nodata)
				xx.push_back(val_model);
			else
				throw NoDataException("No model data at " + ivec[kk].date.toString(Date::ISO, false) + ". Please refine your model output to match the timestamps or enable meteo resampling.");
		}

	} else if (has_model_expression) {

		double te_kk, te_x_km1; //2 substitutions available: index and state value at previous time step
		te_variable te_vars[] = {{"kk", &te_kk, 0, 0}, {"x_km1", &te_x_km1, 0, 0}}; //read: "x_(k-1)"

		int te_err;
		te_expr *expr = te_compile(model_expression.c_str(), te_vars, 2, &te_err); //ready the lazy equation with variables including syntax check

		if (expr) {
			if (model_x0 == IOUtils::nodata) {
				if (be_verbose) std::cerr << "[W] No initial state value x_0 provided for particle filter; using 1st measurement.\n";
				model_x0 = ivec.front()(param);
			}
			xx.push_back(model_x0);

			for (size_t kk = 1; kk < TT; ++kk) { //fill rest of vector
				te_kk = (double)kk;
				te_x_km1 = xx[kk-1];
				const double res = te_eval(expr); //evaluate expression with current substitution values
				xx.push_back(res);
			}
			te_free(expr);

			if (add_model_noise) { //user asks to add noise to model function
				RandomNumberGenerator RNU(model_rng.algorithm); //init the generator for just this scope
				RNU.setDistribution(model_rng.distribution, model_rng.parameters);
				if (!model_rng.seed.empty()) //seed integers given in ini file
					RNU.setState(model_rng.seed);

				for (size_t kk = 0; kk < TT; ++kk)
					xx[kk] += RNU.draw();

			} //endif noise


	    } else {
	    	throw InvalidFormatException("Arithmetic expression '" + model_expression +
	    	    "'could not be evaluated for particle filter; parse error at " + IOUtils::toString(te_err), AT);
	    } //endif expr

	} else {
		throw InvalidArgumentException("No modeled data found in meteo data (expected '"
		    + param_name + model_appendix + "' for the particle filter). No model function found either.");
	} //endif model_expression available

	return true;

}

void FilterParticle::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where("Filters::"+block_name);

	std::string model_rng_algorithm, model_rng_distribution, model_rng_parameters;

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {

		/*** FILTER settings ***/
		if (vecArgs[ii].first == "NO_OF_PARTICLES") {
			IOUtils::parseArg(vecArgs[ii], where, NN);
		} else if (vecArgs[ii].first == "PATH_RESAMPLING") {
			IOUtils::parseArg(vecArgs[ii], where, path_resampling);

		/*** MODEL FUNCTION settings ***/
		} else if (vecArgs[ii].first == "MODEL_APPENDIX") {
			model_appendix = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "MODEL_FUNCTION") {
			model_expression = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "MODEL_X0") {
			IOUtils::parseArg(vecArgs[ii], where, model_x0);
		} else if (vecArgs[ii].first == "MODEL_ADD_NOISE") {
			IOUtils::parseArg(vecArgs[ii], where, add_model_noise);

		/*** MODEL RNG settings ***/
		} else if (vecArgs[ii].first == "MODEL_RNG_ALGORITHM") {
			std::string tmp; //everything RNG will be defaulted if not provided - cf. RNG doc
			IOUtils::parseArg(vecArgs[ii], where, tmp);
			model_rng.algorithm = RandomNumberGenerator::strToRngtype(tmp);
		} else if (vecArgs[ii].first == "MODEL_RNG_DISTRIBUTION") {
			std::string tmp;
			IOUtils::parseArg(vecArgs[ii], where, tmp);
			model_rng.distribution = RandomNumberGenerator::strToRngdistr(tmp); //convert from int to enum RNG_DISTR
		} else if (vecArgs[ii].first == "MODEL_RNG_PARAMETERS") {
			std::string tmp;
			IOUtils::parseArg(vecArgs[ii], where, tmp);
			IOUtils::readLineToVec(tmp, model_rng.parameters);
		} else if (vecArgs[ii].first == "MODEL_RNG_SEED") {
			std::string tmp;
			IOUtils::parseArg(vecArgs[ii], where, tmp);
			readLineToVec(tmp, model_rng.seed);

		/*** MISC settings ***/
		} else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);
		}

	} //endfor ii

}

void FilterParticle::readLineToVec(const std::string& line_in, std::vector<uint64_t>& vec_out)
{ //uint64 type version
	vec_out.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input
	uint64_t val;
	while (!iss.eof()) {
		iss >> std::skipws >> val;
		if (iss.fail())
			throw InvalidFormatException("Unable to parse process noise seed integers for particle filter.", AT);
		vec_out.push_back(val);
	}
}

} //namespace
