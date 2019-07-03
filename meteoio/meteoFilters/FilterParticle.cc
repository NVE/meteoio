#include <meteoio/meteoFilters/FilterParticle.h>
#include <meteoio/tinyexpr.h>

using namespace std;

namespace mio {

static const std::string model_appendix("_MOD");

static const size_t nx = 1; //dimension of the state vector (how many observables?), fixed so far
static const size_t nu = 1; //dimension of process noise, fixed so far
static const size_t nv = 1; //dimension of observation noise, fixed so far

FilterParticle::FilterParticle(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
          : ProcessingBlock(vecArgs, name), filter_alg(SIR), resample_alg(SYSTEMATIC), NN(2000), path_resampling(false),
            model_expression(""), model_x0(IOUtils::nodata), xx(), zz(),
			rng_alg(RandomNumberGenerator::RNG_XOR),
			be_verbose(true)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::second;
}

void FilterParticle::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	const size_t TT = ivec.size(); //number of time steps

	const bool init_success = fill_state(param, ivec);
	if (!init_success)
    	throw NoDataException("The Particle Filter could not initialize the state vector. There is neither good model input data nor a usable model function available.", AT);

	ovec = ivec;
	for (size_t ii = 0; ii < TT; ++ii){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		if (tmp < 0.8){
			tmp = IOUtils::nodata;
		}
	}
}

bool FilterParticle::fill_state(const unsigned int& param, const std::vector<MeteoData>& ivec)
{
	const size_t TT = ivec.size(); //number of time steps
	xx.clear();

	const bool has_model_expression = model_expression.size() != 0; //user has entered a model equation, valid or not

	//check if model data exists in meteo data as per our naming convention:
	const std::string param_name = ivec.front().getNameForParameter(param);
	const size_t param_mod = ivec.front().getParameterIndex(param_name + model_appendix);
	const bool has_param_mod = param_mod != IOUtils::npos; //user has supplied model data alongside the meteo data

	if (has_model_expression && has_param_mod)
		if (be_verbose) std::cerr << "[W] Both model data and model function are supplied for the Particle Filter. Ignoring model function.\n";

	if (has_param_mod) {

		for (size_t kk = 1; kk < TT; ++kk) {
			const double val_model = ivec[kk](param_mod);
			if (val_model != IOUtils::nodata && ivec[kk](param) != IOUtils::nodata)
				xx.push_back(val_model);
			//else
				//throw NoDataException("No model data at " + ivec[kk].date.toString(Date::ISO, false) + ". Please refine your model output to match the timestamps or enable meteo resampling.");
			//std::cout << kk << ": " << ivec[kk](param_mod) << std::endl;
		}

	} else if (has_model_expression) {

		double te_kk, te_x_km1; //2 substitutions available: index and state value at last time step
		te_variable te_vars[] = {{"kk", &te_kk, 0, 0}, {"x_km1", &te_x_km1, 0, 0}};

		int te_err;
		te_expr *expr = te_compile(model_expression.c_str(), te_vars, 2, &te_err); //ready equation with variables including syntax check

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
				//std::cout << kk << ": " << res << std::endl;
			}
			te_free(expr);
	    } else {
	    	throw InvalidNameException("Arithmetic expression '" + model_expression +
	    	    "'could not be evaluated for particle filter; parse error at " + IOUtils::toString(te_err), AT);
	    } //endif expr

	} else {
		throw InvalidArgumentException("No modeled data found in meteo data (expected '"
		    + param_name + model_appendix + "' for the Particle Filter). No model function found either.");
	} //endif model_expression available

	return true;

}

void FilterParticle::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where("Filters::"+block_name);

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="NO_OF_PARTICLES") {
			IOUtils::parseArg(vecArgs[ii], where, NN);
		} else if (vecArgs[ii].first=="PATH_RESAMPLING") {
			IOUtils::parseArg(vecArgs[ii], where, path_resampling);
		} else if (vecArgs[ii].first=="MODEL") {
			model_expression = vecArgs[ii].second;
		} else if (vecArgs[ii].first=="MODEL_X0") {
			IOUtils::parseArg(vecArgs[ii], where, model_x0);
		} else if (vecArgs[ii].first=="VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);

		} else if (vecArgs[ii].first=="RNG_ALGORITHM") {
			if (vecArgs[ii].second == "XOR")
				rng_alg = RandomNumberGenerator::RNG_XOR;
			else if (vecArgs[ii].second == "MTW")
				rng_alg = RandomNumberGenerator::RNG_MTW;
			else if (vecArgs[ii].second == "PCG")
				rng_alg = RandomNumberGenerator::RNG_PCG;
			else
				throw InvalidArgumentException("Unknown RNG algorithm '" + vecArgs[ii].second + "' provided for the Particle Filter.", AT);
		}
	} //endfor ii

}



} //namespace
