/* DO NOT USE this file until you see the usual header here. */

/* - doc piece -
 * RANDOM NUMBER GENERATION
 * Random numbers outside of the insidious standard library.
 * We offer two inherently 32 bit generators, and an inherently 64 bit generator,
 * (although all three need 64 bits space) aswell as some convenience methods.
 *
 * The goal is to have a generator suite that satisfies all needs for statistical filters / Monte Carlo methods (and not
 * more), especially when working within meteoIO. In a way, statistical filters are what ultimately justify this class,
 * and therefore it is meant to be tailored to their needs (and be ANSI C).
 *
 * So, if you are currently using this[1]:
 *    srand( time(NULL) );
 *    return rand() % range;
 *    return rand() / double(RAND_MAX + 1);
 * then switch to meteoIO's RNG. If however you rely heavily on the best quality random numbers, maybe even crypto-secure,
 * there are some links to dedicated libraries in the bibliography.
 * Apart from the generators and distributions, this class aims to take away all the small steps that are often quickly
 * deemed good enough, i.e. generator choice, seeding, saving states, range calculations, ...
 *
 * What it can already do:
 *  - produce quality 64 bit, 32 bit and double random numbers with one simple call
 *  - doubles with uniform and Gaussian distribution
 *  - probability density functions and cumulative distribution functions
 *  - make use of quality hardware and time seeds
 *  - fast downscaling of random numbers to a range
 *  - true floating point random numbers without rounding
 *  - can be resumed from a saved state
 *  - sidesteps some widespread misuse of quick & dirty solutions
 *  - offers a ready-to-use interface for implementing new distributions (or even generators)
 *  - passes statistical tests
 *  - very good benchmarks for the generator cores, memory check passed
 *
 * What's left to do:
 *  - some distributions
 *  - Monte Carlo sampling template for arbitrary distribution functions
 *  - write the doc
 *
 * For developers of statistical filters it may be important to be able to implement custom probability distributions,
 * for example for an empirical nonlinear sensor response. This class tries to be easy to expand in that regard.
 * There are comment markers in the header file and in here leading with "CUSTOM_DIST step #: ..." in the 6 places
 * you need to register your custom distribution functions at. These 6 steps are:
 *  1) Give your distribution a name within meteoIO
 *  2) Put your functions' prototypes in the header
 *  3) Point to your distribution function in the generic setDistribution() function,
 *     and use the interface to the caller to set your distribution parameters
 *  4) Give a small output info string
 *  5) Write your distribution function, its pdf and cdf (if only to throw a not-implemented error)
 *  6) If you want, you can map your parameters to names in the get- and setDistributionParameter() functions.
 *
 * REFERENCES
 * [AS73] Abramowitz, Stegun.
 *        Handbook of Mathematical Functions.
 *        Applied Mathematics Series 55, 10th edition, 1973.
 * [DK81] Donald E. Knuth.
 *        The art of computer programming 2
 *        Addison-Wesley series in computer science and information processing, 2nd edition, 1981.
 * [GM03] George Marsaglia.
 *        Xorshift RNGs.
 *        Journal of Statistical Software, Articles, 8/14, 2003.
 * [MN98] Makoto Matsumoto and Takuji Nishimura.
 *        Mersenne Twister: A 623-dimensionally equidistributed uniform pseudo-random number generator.
 *        ACM Transactions on Modeling and Computer Simulation, 8/1, 1998.
 *        http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 * [MO14] Melissa O'Neill.
 *        PCG: A family of simple fast space-efficient statistically good algorithms
 *        for random number generation.
 *        Harvey Mudd College, 2014.
 *        http://www.pcg-random.org
 * [MT00] G. Marsaglia and W. Tsang.
 *        A simple method for generating gamma variables.
 *        ACM Transactions on Mathematical Software, 26/3, 2000.
 *  [NR3] Press, Teukolsky, Vetterling, Flannery.
 *        Numerical Recipes. The Art of Scientific Computing.
 *        Cambridge University Press, 3rd edition, 2007.
 * [PE97] Pierre L'Ecuyer.
 *        Distribution properties of multiply-with-carry random number generators.
 *        Mathematics of Computation, 66/218i, 1997.
 * [PE99] Pierr L'Ecuyer.
 *        Tables of linear congruential generators of different sizes and good
 *        lattice structure.
 *        Mathematics of Computation, 68/225, 1999.
 *        (Errata for the paper at
 *        https://www.iro.umontreal.ca/~lecuyer/myftp/papers/latrules99Errata.pdf, read on 18-10-23) 
 * [TC14] Taylor R. Campbell.
 *        http://mumble.net/~campbell/tmp/random_real.c (read on 18-10-23)
 *
 * APPENDIX
 * [1] Why is
 *      srand( time(NULL) );
 *      return rand() % range;
 *      return rand() / double(RAND_MAX + 1);
 * bad?
 * - A purely linear congruential RNG has purely bad statistical qualities
 * - Quiet type collision between time() and srand()
 * - (% range) distorts the distribution at the borders and (range + 1) should be used
 * - Careful not to hit RAND_MAX = INT_MAX, maybe ((double)RAND_MAX) + 1.
 * Why is
 *      std::random_device RNG;
 *      std::seed_seq seed{RNG()};
 *      std::mt19937 RNG_MT(seed);
 * not good?
 * - A 624 state Mersenne Twister is seeded with a single 32 bit value, not a sequence
 * - Leads to statistical flaws; some numbers are never drawn
 * - std::seed_seq isn't a bijection like it's supposed to be
 * - Can produce zero-state
 *
 * Author: Michael Reisecker, 2018-10
*/

#include <cmath>
#include <fstream> //for hardware seed
#include <limits> //for numeric_limits
#include <sstream> //for toString()

#include <meteoio/IOUtils.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

namespace mio { //the holy land

///////////////////////////////////////////////////////////////////////////////
//    RANDOM NUMBER GENERATOR class                                          //
///////////////////////////////////////////////////////////////////////////////

/* CONSTRUCTOR */

/* -doc piece-
 * Each time a RNG is constructed, it auto-seeds from hardware noise, or if that fails the system time.
 * Manually seeding the generator is done after the fact with RNG.setState(), for example, to resume
 * experiments after the state was saved via RNG.getState().
 * Finally, we offer the getUniqueSeed() function to the outside, so if someone has set up their
 * calculations with a grandfathered in, better, faster, ... RNG we can at least help with the seeding.
 */
RandomNumberGenerator::RandomNumberGenerator(const RNG_TYPE& type, const RNG_DISTR& distribution,
    const std::vector<double>& distribution_params) :
    rng_core(RngFactory::getCore(type)),
    rng_type(type),
    rng_distribution(distribution),
    DistributionParameters(distribution_params),
    rng_muller_generate(false),
    rng_muller_z1(0.),
    doubFunc(&RandomNumberGenerator::doubUniform),
    pdfFunc(&RandomNumberGenerator::pdfUniform),
    cdfFunc(&RandomNumberGenerator::cdfUniform)
{
	setDistribution(distribution, distribution_params); //checks if the parameters fit the distribution
}

RandomNumberGenerator::RandomNumberGenerator(const RandomNumberGenerator& source) : 
    rng_core(RngFactory::getCore(source.rng_type)),
    rng_type(source.rng_type),
    rng_distribution(source.rng_distribution),
    DistributionParameters(source.DistributionParameters),
    rng_muller_generate(source.rng_muller_generate),
    rng_muller_z1(source.rng_muller_z1),
    doubFunc(source.doubFunc),
    pdfFunc(source.pdfFunc),
    cdfFunc(source.cdfFunc)
{
	std::vector<uint64_t> transfer_states;
	source.getState(transfer_states);
	setState(transfer_states); //a little different per generator
}

RandomNumberGenerator::~RandomNumberGenerator()
{
	delete rng_core; //TODO: good practice to guard this further?
}

RandomNumberGenerator& RandomNumberGenerator::operator=(const RandomNumberGenerator& source)
{
	if (this != &source) {
		rng_core = RngFactory::getCore(source.rng_type);
		rng_type = source.rng_type;
		rng_distribution = source.rng_distribution;
		DistributionParameters = source.DistributionParameters;
		rng_muller_generate = source.rng_muller_generate;
		rng_muller_z1 = source.rng_muller_z1;
		doubFunc = source.doubFunc;
		pdfFunc = source.pdfFunc;
		cdfFunc = source.cdfFunc;

		std::vector<uint64_t> transfer_states;
		source.getState(transfer_states);
		setState(transfer_states);
	}
	return *this;
}

/* PUBLIC FUNCTIONS */
uint64_t RandomNumberGenerator::int64()
{
	return rng_core->int64();
}

uint32_t RandomNumberGenerator::int32()
{
	return rng_core->int32();
}

/* - doc piece -
 * The RNG.doub() function returns a double within [0, 1] that is rounded to the nearest 1/2^64th.
 * If you absolutely need to control the properties further, look into the RNG.doub(RNG_BOUND bound, bool true_double)
 * function call:
 * You can call the doub() function with an RNG_BOUND argument and choose from RNG_AINCBINC [0, 1], RNG_AINCBEXC [0, 1),
 * RNG_AEXBINC (0, 1] and RNG_AEXCBEXC (0, 1). This can be only done for the uniform distribution, where it's clear
 * what the borders are.
 * You can also set true_double to use an algorithm that calculates doubles in [0, 1] without the usual limitation of
 * floating point randoms being on a grid (but then you must use RNG_AINCBINC and guard that in your own code).
 * Example 1:
 *     double rr = RNG.doub(RNG_AEXCBINC); //make sure it's not 0
 *     rr = log(rr);
 * Example 2:
 *     double rr;
 *     do {
 *         rr = RNG.doub(RNG_AINCBINC, true); //get a random float on continuous axis
 *     } while (rr == 0.); //make sure it's not 0
 *     rr = log(rr)
 */
double RandomNumberGenerator::doub()
{
	return (this->*doubFunc)();
}

double RandomNumberGenerator::doub(const RNG_BOUND& bounds, const bool& true_double)
{
	if (rng_distribution != RNG_UNIFORM)
		throw InvalidArgumentException("RNG: Boundary exclusion not implemented for non-uniform distributions. Please use doub() and check against the limits in your own code.", AT);
	
	double rr = true_double? RngCore::trueDoub() : doubUniform();
	if (bounds != RNG_AINCBINC) {
		while ( (bounds == RNG_AINCBEXC && rr == 1.) ||
		        (bounds == RNG_AEXCBINC && rr == 0.) ||
		        (bounds == RNG_AEXCBEXC && (rr == 0. || rr == 1.)) )
			rr = true_double? RngCore::trueDoub() : doubUniform();	
	}

	return rr;

	////some ways to convert integers to doubles with borders (within limitations):	
	//return 5.42101086242752217E * int64(); //divide by 2^64-1 [0, 1]
	////(boundaries run into precision problems for 64 bits but aren't realistically happening anyway if the rng performs)
	//return int32() * (1. / 4294967295.); //divide by 2^32-1 for [0, 1]
	//return int32() * (1. / 4294967296.); //divide by 2^32 for [0, 1)
	//return ((int64() >> 12) + 0.5) * (1./4503599627370496.); //for (0, 1)
	//return int32() / (double)(UINT32_MAX + 1); //[0, 1]
	//uint32_t a = int32() >> 5, b = int32() >> 6; 
	//return (a * 67108864. + b) * (1. / 9007199254740992.); //53 bits resolution
}

double RandomNumberGenerator::draw() //convenience call with no gimmicks
{
	return doubUniform();
}

double RandomNumberGenerator::pdf(const double& xx) //probability density function of selected distribution
{
	return (this->*pdfFunc)(xx);
}

double RandomNumberGenerator::cdf(const double& xx) //cumulative distribution function of selected distribution
{
	return (this->*cdfFunc)(xx);
}

/* - doc piece -
 * Note that whatever you do, for an arbitrary count of random numbers you cannot downscale them and keep the distribution intact
 * (although "non-trivial" methods are under investigation) due to the Pigeonhole principle (https://en.wikipedia.org/wiki/Pigeonhole_principle).
 * The only way not to distort the (uniform) distribution is to generate lots of numbers and reject
 * out of boundary values. This is done by the trueRange32(...) function with a default 1e6 tries before
 * resorting to downscaling. You can crank this up, but to state the obvious if the range gets small this gets
 * costly quickly. The bitshift-methods here only avoid the slow modulo and its possible inherent bias.
 */
uint64_t RandomNumberGenerator::range64(const uint64_t& aa, const uint64_t& bb) 
{ //needs 64 bits space
	const uint64_t diff = &bb + 1 - &aa; //[a, b]
	if (diff < 0xffffffff) //avoid UINT32_MAX (platform and compiler dependent whether it's there)
		return (uint64_t)(range32( (uint32_t)aa, (uint32_t)bb ));

	//multiply that keeps only the top 64 bits of the resulting 128 bit number:
	const uint64_t rlow  = int32();
	const uint64_t rhigh = int32();
	const uint64_t vlow  = diff & 0xffffffff;
	const uint64_t vhigh = diff >> 32;

	uint64_t rn = ((rhigh * vlow) >> 32 );
		rn += ((rlow * vhigh) >> 32 ); 
		rn += (rhigh * vhigh);
		rn += aa;
	return rn;
}

uint32_t RandomNumberGenerator::range32(const uint32_t& aa, const uint32_t& bb)
{ //also needs 64 bits space
	uint64_t rn = int32();
	rn *= (bb + 1 - aa); //[a, b]
	uint64_t tmp = (( rn >> 32 ) + aa); //keep only the bits above 32
	return (uint32_t)tmp;
}

bool RandomNumberGenerator::trueRange32(const uint32_t& aa, const uint32_t& bb, uint32_t& result, const unsigned int& nmax)
{ //only 32 bits to have more chances
	uint64_t it(0);
	uint32_t rn(0);
	do {
		++it;
		if (it > nmax) {
			result = range32(aa, bb); //fallback
			return false;
		}
		rn = int32();
	} while ( (rn < aa) || (rn > bb) ); //[a, b]

	result = rn;
	return true;
}

void RandomNumberGenerator::getState(std::vector<uint64_t>& ovec_seed) const
{
	rng_core->getState(ovec_seed);
}

void RandomNumberGenerator::setState(const std::vector<uint64_t>& ivec_seed)
{
	rng_core->setState(ivec_seed);
}

/* -doc piece -
 * - DISTRIBUTION -
 * When you change the distribution, you switch to a completely new one.
 * The distribution parameters _have to be_ provided each time, or they will be defaulted.
 *
 * - DISTRIBUTION PARAMETERS -
 * You can set distribution parameters in two ways:
 *  1) By setting (getting) them one by one after a distribution has been set:
 *        RNG.setDistribution(mio::RandomNumberGenerator::RNG_GAUSS);
 *        RNG.setDistributionParameter("mean", 5.);
 *        RNG.setDistributionParameter("sigma", 2.);
 *        double mean_out = RNG.getDistributionParameter("mean");
 *  2) Via accessing the DistributionParameters vector directly. This is a generic interface between the
 *    implementation of the distribution algorithm and the end user allowing for arbitrary parameter passing.
 *    A std::vector<double> is provided by the user with input parameters, or it has the output stored to it.
 *    This vector is given to the set- or getDistribution() call (the latter also returns the distribution type).
 *    Set distribution and parameters:
 *        std::vector<double> distribution_params;
 *        double gamma_in = 0.7;
 *        distribution_params.push_back(gamma_in);
 *        RNG.setDistribution(mio::RandomNumberGenerator::RNG_GAMMA, distribution_params); //set distribution with params
 *    Be aware that if you forget this, default parameters may be set without a warning.
 *    Get distribution and parameters:
 *        distribution_params.clear();
 *        const mio::RandomNumberGenerator::RNG_DISTR dist_t = RNG.getDistribution(distribution_params);
 *        const double gamma_out = distribution_params[0]; //check doc for indices
 *
 * RNG_UNIFORM:
 *     no parameters
 * RNG_GAUSS = RNG_NORMAL:
 *     #1:  "mean" ... center of curve (default: 0)
 *     #2: "sigma" ... standard deviation (default: 1)
 * RNG_GAMMA:
 *     #1: "alpha" ... shape parameter 1 (default: 1)
 *     #2:  "beta" ... shape parameter 2 (default: 1)
 * RNG_CHISQUARE:
 *     #1:    "nu" ... number of degrees of freedom (default: 1)
 * RNG_STUDENTT:
 *     #1:    "nu" ... shape parameter (default: 1)
 *     #2:  "mean" ... center of curve (default: 0)
 *     #3: "sigma" ... standard deviation (default: 1)
 * RNG_BETA:
 *     #1: "alpha" ... shape parameter 1 (default: 1)
 *     #2:  "beta" ... shape parameter 2 (default: 1)
 * RNG_F:
 *     #1:   "nu1" ... degrees of freedom in numerator (default: 1)
 *     #2:   "nu2" ... degrees of freedom in denominaator (default: 1)
*/
RandomNumberGenerator::RNG_DISTR RandomNumberGenerator::getDistribution(std::vector<double>& vec_params) const
{
	vec_params = DistributionParameters;
	return rng_distribution;
}

//CUSTOM_DIST step 3/6: Add a case for your distribution here and set your functions from step 2, aswell as defaults
//for all parameters the distribution needs via the vector DistributionParameters. Please make sure all mandatory ones
//are described properly. Cf. notes in doc setDistribution().
void RandomNumberGenerator::setDistribution(const RNG_DISTR& distribution, const std::vector<double>& vec_params)
{
	rng_distribution = distribution;
	DistributionParameters.clear(); //we enforce a complete reset to make things clear
	switch (distribution) {
	case RNG_UNIFORM:
		this->doubFunc = &RandomNumberGenerator::doubUniform;
		this->pdfFunc = &RandomNumberGenerator::pdfUniform;
		this->cdfFunc = &RandomNumberGenerator::cdfUniform;
		break;
	case RNG_GAUSS: case RNG_NORMAL:
		this->doubFunc = &RandomNumberGenerator::doubGauss;
		this->pdfFunc = &RandomNumberGenerator::pdfGauss;
		this->cdfFunc = &RandomNumberGenerator::cdfGauss;
		if (vec_params.size() == 0) { //no input -> choose defaults (or throw error)
			DistributionParameters.push_back(0.); //def. mean = 0
			DistributionParameters.push_back(1.); //def. standard deviation = 1
		} else if (vec_params.size() == 2) { //user inits with right number of parameters
			DistributionParameters = vec_params;
		} else {
			throw InvalidArgumentException("RNG: Incorrect number of input parameters for distribution (length of input vector). Expected: 2 (mean, sigma).", AT);
		}
		break;
	case RNG_GAMMA:
		this->doubFunc = &RandomNumberGenerator::doubGamma;
		this->pdfFunc = &RandomNumberGenerator::pdfNotImplemented; //this is some 2000 lines of code in GSL
		this->cdfFunc = &RandomNumberGenerator::cdfNotImplemented;
		if (vec_params.size() == 0) {
			DistributionParameters.push_back(1.); //def. alpha = 1
			DistributionParameters.push_back(1.); //def. beta = 1
		} else if (vec_params.size() == 2) {
			DistributionParameters = vec_params;
		} else {
			throw InvalidArgumentException("RNG: Incorrect number of input parameters for distribution (length of input vector). Expected: 2 (alpha, beta).", AT);
		}
		break;
	case RNG_CHISQUARE:
		this->doubFunc = &RandomNumberGenerator::doubChiSquare;
		this->pdfFunc = &RandomNumberGenerator::pdfNotImplemented;
		this->cdfFunc = &RandomNumberGenerator::cdfNotImplemented;
		if (vec_params.size() == 0) {
			DistributionParameters.push_back(1.); //def. nu = 1
		} else if (vec_params.size() == 1) {
			DistributionParameters = vec_params;
		} else {
			throw InvalidArgumentException("RNG: Incorrect number of input parameters for distribution (length of input vector). Expected: 1 (nu).", AT);
		}
		break;
	case RNG_STUDENTT:
		this->doubFunc = &RandomNumberGenerator::doubStudentT;
		this->pdfFunc = &RandomNumberGenerator::pdfNotImplemented;
		this->cdfFunc = &RandomNumberGenerator::cdfNotImplemented;
		if (vec_params.size() == 0) {
			DistributionParameters.push_back(1.); //def. nu = 1
			DistributionParameters.push_back(0.); //def. mean = 0
			DistributionParameters.push_back(1.); //def. sigma = 1
		} else if (vec_params.size() == 3) {
			DistributionParameters = vec_params;
		} else {
			throw InvalidArgumentException("RNG: Incorrect number of input parameters for distribution (length of input vector). Expected: 3 (nu, mean, sigma).", AT);
		}
		break;
	case RNG_BETA:
		this->doubFunc = &RandomNumberGenerator::doubBeta;
		this->pdfFunc = &RandomNumberGenerator::pdfNotImplemented;
		this->cdfFunc = &RandomNumberGenerator::cdfNotImplemented;
		if (vec_params.size() == 0) {
			DistributionParameters.push_back(1.); //def. alpha = 1
			DistributionParameters.push_back(1.); //def. beta = 1
		} else if (vec_params.size() == 2) {
			DistributionParameters = vec_params;
		} else {
			throw InvalidArgumentException("RNG: Incorrect number of input parameters for distribution (length of input vector). Expected: 2 (alpha, beta).", AT);
		}
		break;
	case RNG_F:
		this->doubFunc = &RandomNumberGenerator::doubF;
		this->pdfFunc = &RandomNumberGenerator::pdfNotImplemented;
		this->cdfFunc = &RandomNumberGenerator::cdfNotImplemented;
		if (vec_params.size() == 0) {
			DistributionParameters.push_back(1.); //def. nu1 = 1
			DistributionParameters.push_back(1.); //def. nu2 = 1
		} else if (vec_params.size() == 2) {
			DistributionParameters = vec_params;
		} else {
			throw InvalidArgumentException("RNG: Incorrect number of input parameters for distribution (length of input vector). Expected: 2 (nu1, nu2).", AT);
		}
		break;
	default:
		throw InvalidArgumentException("RNG: This distribution is not implemented. Check your RNG_DISTR.", AT);
	}
	rng_muller_generate = false; //clear cached Box-Muller value 
}

/*
 * Hardcoded mapping of parameters to values, since giving access by name slows down distribution functions that look them up
 * by more than 300%. So, we list the indices that lead to certain distribution parameters per hand here.
 * In these two functions, a distribution is meant to provide a nice interface to the user, but it's not mandatory to do so.
 * (About 1% speed loss per value in vector as opposed to global constants in uniform distribution.)
 */
double RandomNumberGenerator::getDistributionParameter(const std::string& param_name) const
{ //convenience call taking away the vector from the user
	const std::string str_param_error("RNG: Distribution parameter " + param_name +
	    " not set when trying to read. If you are sure you set the correct distribution before encountering this error please notify the developers.");

	//we have supposedly made sure in setDistribution() that we have correctly mapped vector indices
	switch (rng_distribution) {
	case RNG_UNIFORM:
		throw InvalidArgumentException("RNG: No parameters associated with uniform distribution. Set another one via RNG.setDistribution(RNG_GAUSS) for example.", AT);
		break;
	case RNG_GAUSS: case RNG_NORMAL:
		if (IOUtils::strToLower(param_name) == "mean")
			return DistributionParameters.at(0);
		else if (IOUtils::strToLower(param_name) == "sigma")
			return DistributionParameters.at(1);
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_GAMMA:
		if (IOUtils::strToLower(param_name) == "alpha")
			return DistributionParameters.at(0);
		else if (IOUtils::strToLower(param_name) == "beta")
			return DistributionParameters.at(1);
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_CHISQUARE:
		if (IOUtils::strToLower(param_name) == "nu")
			return DistributionParameters.at(0);
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_STUDENTT:
		if (IOUtils::strToLower(param_name) == "nu")
			return DistributionParameters.at(0);
		else if (IOUtils::strToLower(param_name) == "mean")
			return DistributionParameters.at(1);
		else if (IOUtils::strToLower(param_name) == "sigma")
			return DistributionParameters.at(2);
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_BETA:
		if (IOUtils::strToLower(param_name) == "alpha")
			return DistributionParameters.at(0);
		else if (IOUtils::strToLower(param_name) == "beta")
			return DistributionParameters.at(1);
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_F:
		if (IOUtils::strToLower(param_name) == "nu1")
			return DistributionParameters.at(0);
		else if (IOUtils::strToLower(param_name) == "nu2")
			return DistributionParameters.at(1);
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	default:	
		throw InvalidArgumentException("RNG: Friendly parameter fetching not implemented yet for this distribution. Use RNG_TYPE type = RNG.getDistribution(vector<double> output) and check the doc for the indices.", AT);
	} //end switch
	return 0.;
}

//CUSTOM_DIST step 6/6: Convenience mapping of your distribution parameters to names (not mandatory if you have many) 
void RandomNumberGenerator::setDistributionParameter(const std::string& param_name, const double& param_val)
{ //convenience
	const std::string str_param_error("RNG: Distribution parameter " + param_name + " not available when trying to set. If you are sure you set the correct distribution before encountering this error please notify the developers.");

	switch (rng_distribution) {
	case RNG_UNIFORM:
		throw InvalidArgumentException("RNG: No parameters to set for uniform distribution. Switch to another one via RNG.setDistribution(RNG_GAUSS).", AT);
		break;
	case RNG_GAUSS: case RNG_NORMAL:
		if (IOUtils::strToLower(param_name) == "mean")
			DistributionParameters.at(0) = param_val;
		else if (IOUtils::strToLower(param_name) == "sigma")
			DistributionParameters.at(1) = param_val;
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_GAMMA:
		if (IOUtils::strToLower(param_name) == "alpha")
			DistributionParameters.at(0) = param_val;
		else if (IOUtils::strToLower(param_name) == "beta")
			DistributionParameters.at(1) = param_val;
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_CHISQUARE:
		if (IOUtils::strToLower(param_name) == "nu")
			DistributionParameters.at(0) = param_val;
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_STUDENTT:
		if (IOUtils::strToLower(param_name) == "nu")
			DistributionParameters.at(0) = param_val;
		else if (IOUtils::strToLower(param_name) == "mean")
			DistributionParameters.at(1) = param_val;
		else if (IOUtils::strToLower(param_name) == "sigma")
			DistributionParameters.at(2) = param_val;
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_BETA:
		if (IOUtils::strToLower(param_name) == "alpha")
			DistributionParameters.at(0) = param_val;
		else if (IOUtils::strToLower(param_name) == "beta")
			DistributionParameters.at(1) = param_val;
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	case RNG_F:
		if (IOUtils::strToLower(param_name) == "nu1")
			DistributionParameters.at(0) = param_val;
		else if (IOUtils::strToLower(param_name) == "nu2")
			DistributionParameters.at(1) = param_val;
		else
			throw InvalidArgumentException(str_param_error, AT);
		break;
	default:
		throw InvalidArgumentException("RNG: Friendly parameter setting not implemented yet for this distribution. Use RNG.setDistribution(vector<double> input) and check the doc for the indices.", AT);
	} //end switch
}

bool RandomNumberGenerator::getHardwareSeedSuccess() const
{ //has the RNG successfully seeded from hardware noise?
	return rng_core->hardware_seed_success;
}

bool RandomNumberGenerator::getUniqueSeed(uint64_t& store) const
{ //allow for outside calls to the seeding function if someone wants only that for their own RNG
	return rng_core->getUniqueSeed(store);
}

std::string RandomNumberGenerator::toString()
{
	std::stringstream ss;
	switch (rng_type) {
	case RNG_XOR:
		ss << "Name: RNG_XOR\n";
		ss << "Family: Xor, shift, multiply\n";
		ss << "Size: 64 bit\n";
		ss << "Period: ~3.138*10^57\n";
		break;
	case RNG_PCG:
		ss << "Name: RNG_PCG\n";
		ss << "Family: Permuted linear congruential\n";
		ss << "Size: 32 bit\n";
		ss << "Period: ~2^64 ~ 1.8e19\n";
		break;
	case RNG_MTW:
		ss << "Name: Mersenne Twister\n";
		ss << "Family: Twisted feedback shift register\n";
		ss << "Size: 32 bit\n";
		ss << "Period: 2^19937-1 ~ 4.3e6001\n";
		break;
	}
	ss << "Hardware seeded: " << (getHardwareSeedSuccess()? "yes" : "no") << "\n";

	switch (rng_distribution) {
	case RNG_UNIFORM: 
		ss << "Distribution: Uniform\n";
		break;
	case RNG_GAUSS: case RNG_NORMAL:
		ss << "Distribution: Gauss\n";
		ss << "Mean: " << DistributionParameters.at(0)
		   << ", sigma: " << DistributionParameters.at(1) << "\n";
		break;
	case RNG_GAMMA:
		ss << "Distribution: Gamma\n";
		ss << "Alpha: " << DistributionParameters.at(0)
		   << ", beta: " << DistributionParameters.at(1) << "\n";
		break;
	case RNG_CHISQUARE:
		ss << "Distribution: Chi-Square\n";
		ss << "Nu: " << DistributionParameters.at(0) << "\n";
		break;
	case RNG_STUDENTT:
		ss << "Distribution: Student-T\n";
		ss << "Nu: " << DistributionParameters.at(0)
		   << ", mean: " << DistributionParameters.at(1)
		   << ", sigma: " << DistributionParameters.at(2) << "\n";
		break;
	case RNG_BETA:
		ss << "Distribution: Beta\n";
		ss << "Alpha: " << DistributionParameters.at(0)
		   << ", beta: " << DistributionParameters.at(1) << "\n";
		break;
	case RNG_F:
		ss << "Distribution: F (Fisher)\n";
		ss << "Nu1: " << DistributionParameters.at(0)
		   << ", nu2: " << DistributionParameters.at(1) << "\n";
		break;
//CUSTOM_DIST step 4/6: Give a small info string in this function. 
	default:
		ss << "Distribution: custom\n";	
	}

	return ss.str();
}

/* PRIVATE FUNCTIONS */
double RandomNumberGenerator::doubUniform()
{
	const uint64_t rn = rng_core->int64();
	return RngCore::doubFromInt(rn);
}

double RandomNumberGenerator::pdfUniform(const double& xx) const
{
	//in a given interval, it is 1/(b-a) or 0 outside; for [0, 1] this is 1
	return xx*0.+ 1.; //guard against "unused parameter" compiler warning
}

double RandomNumberGenerator::cdfUniform(const double& xx) const
{
	return xx; //in [0, 1]
}

double RandomNumberGenerator::doubGauss() //Gauss double with user-set parameters
{ //Gauss double with user set distribution parameters
	const double mean = DistributionParameters[0];
	const double sigma = DistributionParameters[1];
	
	return doubGaussKernel(mean, sigma);
}

/* The "Kernel" versions are separated for deviates that are needed with fixed parameters (i. e. other than
 * stored in DistributionParameters) by other deviates, and for recursive calls for transformations. */
double RandomNumberGenerator::doubGaussKernel(const double& mean, const double& sigma) //Box-Muller
{ //generate 2 uniform doubles and transform to 2 Gaussians
	const double eps = std::numeric_limits<double>::min();
	
	//2 independent numbers are generated at once -> new calculation every 2nd call
	rng_muller_generate = !rng_muller_generate;
	if (!rng_muller_generate)
		return rng_muller_z1 * sigma + mean;

	double x1, x2;
	do {
		x1 = doubUniform();
		x2 = doubUniform();
	} while (x1 <= eps);

	double z0;
	z0 = sqrt(-2. * log(x1)) * cos(2.*M_PI * x2); //TODO: make this a little faster
	rng_muller_z1 = sqrt(-2. * log(x1)) * sin(2.*M_PI * x2);
	return z0 * sigma + mean;
} //http://mathworld.wolfram.com/Box-MullerTransformation.html


double RandomNumberGenerator::pdfGauss(const double& xx) const
{ //Gauss curve around mean and with standard deviation at point xx (probability density function)
	const double mean = DistributionParameters[0];
	const double sigma = DistributionParameters[1];

	return 1. / (sigma * sqrt(2. * M_PI)) * exp( -pow((xx-mean), 2.) / (2. * sigma*sigma) );
}

double RandomNumberGenerator::cdfGauss(const double& xx) const
{ //cumulative distribution function for a Gauss curve at point xx
	const double mean = DistributionParameters[0];
	const double sigma = DistributionParameters[1];
	const double xabs = fabs(xx - mean) / sigma; //formula is for mu=0 and sigma=1 --> transform
	
	const double bb[5] = { //|error| < 7.5e-8
	    0.319381530, 
	   -0.356563782,
	    1.781477937,
	   -1.821255978,
	    1.330274429
	};
	const double pp = 0.2316419;

	//Ref. [AS72] formula 26.2.17:
	const double tt = 1. / (1. + pp * xabs);
	//use the Gauss-pdf without mu and sigma as the transformation is nonlinear:
	double yy = 1. - ( ((((bb[4]*tt + bb[3])*tt) + bb[2])*tt + bb[1])*tt + bb[0] )*tt *
	    1. / (sqrt(2. * M_PI)) * exp(-(xabs*xabs) / 2.); //Horner's method to evaluate polynom

	//formula is for mu=0 and x > 0, but it is symmetric around mu and can be shifted trivially:
	const int sign = (xx > mean) - (xx < mean); 
	return 0.5 + sign * (yy - 0.5);
}

/* - doc piece - 
 * The Gamma distribution as proposed by Ref. [MT00].
 * The GNU Scientific Library also implements its Gamma-distribution like this but chooses a slightly different pdf.
 */
double RandomNumberGenerator::doubGamma() //Ref. [MT00]
{
	const double alpha = DistributionParameters[0];
	const double beta = DistributionParameters[1];

	return doubGammaKernel(alpha, beta); //needs to be able to change a, b
}

double RandomNumberGenerator::doubGammaKernel(const double& alpha, const double& beta)
{ //Gamma deviates
	if (alpha < 1.) { //Gamma(a, b) ~ Gamma(a+1, b) * U^(1/a), U ~ (0, 1)
		double ru;
		do {
			ru = doubUniform();
		} while (ru == 0. || ru == 1.);
		return doubGammaKernel(1. + alpha, beta) * pow(ru, 1. / alpha); //cf. GSL/randist/gamma.c
	}

	const double dd = alpha - 1. / 3.;
	const double cc = 1. / sqrt(9. * dd);
	double rn, vv;
	while (true) {
		do {
			rn = doubGaussKernel(0., 1.); //normal double
			vv = 1. + cc * rn;
		} while (vv <= 0.);
		vv = vv * vv * vv;
		double uu = doubUniform(); //uniform double
		if ( uu < 1. - 0.0331 * rn*rn * rn*rn )
			break;
		if ( log(uu) < 0.5 * rn * rn + dd * (1. - vv + log(vv)) )
			break;
	}
	return dd * vv / beta; //Gamma(a, b) ~ Gamma(a, 1) / b
}

double RandomNumberGenerator::doubChiSquare()
{ //ChiSquare(nu) ~ Gamma(nu/2, 1/2)
	const double nu = DistributionParameters[0];
	return doubGammaKernel(nu/2., 0.5);
}

double RandomNumberGenerator::doubStudentT()
{ //Student-t deviates
	const double nu = DistributionParameters[0];
	const double mean = DistributionParameters[1];
	const double sigma = DistributionParameters[2];

	//There's a formula similar to Box-Muller, but it doesn't generate 2 values at once,
	//diminishing this advantage. So, we use what we already have:
	//x ~ N(0, 1), y ~ Gamma(nu/2, 1/2) --> x*sqrt(nu/y) = StudentT(nu, 0, 1)
	const double xx = doubGaussKernel(0., 1.);
	const double yy = doubGammaKernel(nu/2., 0.5);
	const double st = xx * sqrt(nu / yy);

	//s0 ~ StudentT(nu, 0, 1) --> mean + sigma*s0 ~ StudentT(nu, mean, sigma)
	return mean + sigma * st;
}


double RandomNumberGenerator::doubBeta()
{ //Beta deviates
	const double alpha = DistributionParameters[0];
	const double beta = DistributionParameters[1];
	return doubBetaKernel(alpha, beta);	
}

double RandomNumberGenerator::doubBetaKernel(const double& alpha, const double& beta)
{ //TODO: maybe some limits assertions
	//x ~ Gamma(a, 1), y ~ Gamma(b, 1) --> x/(x+y) ~ Beta(a, b)
	const double xx = doubGammaKernel(alpha, 1.);
	const double yy = doubGammaKernel(beta, 2.);
	return xx / (xx + yy);
}

double RandomNumberGenerator::doubF()
{ //Fisher deviates
	const double nu1 = DistributionParameters[0];
	const double nu2 = DistributionParameters[1];

	//x ~ Beta(1/2 nu1, 1/2 nu2) --> nu2*x/(nu1*(1-x)) ~ F(nu1, nu2)
	const double xx = doubBetaKernel(0.5*nu1, 0.5*nu2);
	return nu2 * xx / (nu1 * (1. - xx));
}

//CUSTOM_DIST step 5/6: Implement your 3 functions for the distribution, its pdf and cdf here, corresponding to,
//for example, doubGauss(), pdfGauss() and cdfGauss(). Implement all of them and throw an appropriate error
//if it's not actually there. Please also properly document them here.

double RandomNumberGenerator::pdfNotImplemented(const double& xx) const
{ //pdfs and cdfs are often very hard - we only implement them as needed
	throw InvalidArgumentException("RNG: Probability density function (pdf) not implemented for this distribution.", AT);
	return xx;
}

double RandomNumberGenerator::cdfNotImplemented(const double& xx) const
{
	throw InvalidArgumentException("RNG: Cumulative distribution function (cdf) not implemented for this distribution.", AT);
	return xx;
}

///////////////////////////////////////////////////////////////////////////////
//    XOR GENERATOR class                                                    //
///////////////////////////////////////////////////////////////////////////////

/* CONSTRUCTOR */
RngXor::RngXor() : state(0), uu(0), vv(0), ww(0)
{
	hardware_seed_success = initAllStates();
}

/* PUBLIC FUNCTIONS */

/* - doc piece -
 * Random Number Generator 1 - XOR
 * - Generator with xor, shift and multiplication
 *   This is a fast combined generator that should be suitable for all but very special Monte Carlo applications.
 *   Since more than one internal states are being propagated and combined to the output, this makes it
 *   somewhat less predictable than similar generators.
 * - Seed with any value except vv.
 * - Facts:
 *   size: 64 bit, period: ~3.138e57
 */
uint64_t RngXor::int64()
{
	//first, a linear congruential generator with good figures of merit:
	uu = uu * 7664345821815920749LL + 6204829405619482337LL; //Ref. [PE99] + arb. odd value
	//64 bit xorshift, using one of the empirical triplets preserving order 2^64-1:
	vv ^= (vv << 13); vv ^= (vv >> 7); vv ^= (vv << 17); //Ref. [GM03]
	//multiply with carry:
	ww = 3874257210U * (ww & 0xffffffff) + (ww >> 32); //Ref. [PE97]
	//xorshift on the other states:
	uint64_t xx = uu ^ (uu << 21); xx ^= xx >> 35; xx ^= xx << 4; //Ref. [NR3]
	return (xx + vv) ^ ww;
}

uint32_t RngXor::int32()
{
	return (uint32_t)int64();
}

void RngXor::getState(std::vector<uint64_t>& ovec_seed) const
{
	ovec_seed.clear();
	ovec_seed.push_back(state);
	ovec_seed.push_back(uu);
	ovec_seed.push_back(vv);
	ovec_seed.push_back(ww);
}

void RngXor::setState(const std::vector<uint64_t>& ivec_seed)
{
	if (ivec_seed.size() != 4)
		throw InvalidArgumentException("RNG: Unexpected number of seeds for this generator (needed: 4)", AT);

	state = ivec_seed[0];
	uu = ivec_seed[1];
	vv = ivec_seed[2];
	ww = ivec_seed[3];
}

/* PRIVATE FUNCTIONS */
bool RngXor::initAllStates() //initial XOR-generator states
{
	state = 4; //if everything else fails, seed with random die throw as per a fair experiment
	const bool hardware_success = getUniqueSeed(state);
	vv = 4101842887655102017LL;
	ww = 1;
	uu = state ^ vv; int64();
	vv = uu; int64();
	ww = vv; int64();
	return hardware_success; //pass through whether we had hardware entropy to the constructor
}

///////////////////////////////////////////////////////////////////////////////
//    PCG GENERATOR class                                                    //
///////////////////////////////////////////////////////////////////////////////

/* CONSTRUCTOR */
RngPcg::RngPcg() : state(0), inc(0)
{
	hardware_seed_success = initAllStates();
}

/* PUBLIC FUNCTIONS */
uint64_t RngPcg::int64() //for PCG, draw two 32 bit numbers and combine them to one 64 bit nr
{
	const uint32_t lowpart = int32();
	const uint32_t highpart = int32();
	return RngCore::combine32to64(lowpart, highpart);
}

/* - doc piece -
 * Random Number Generator 2 - PCG
 * - Permuted linear congruential generator by Prof. Melissa O'Neill (Ref. [MO14])
 *   Range is overestimated, and this generator performs very well in statistical tests, i. e. it is less
 *   predictable than related generators. Even smaller versions with only 32 bit entropy pass SmallCrunch,
 *   which is only barely theoretically possible.
 *   The key element is the hashing function from the internal states to the random number.
 *   states. The algorithm author describes this RNG family in her paper (Ref. [MO14]) and offers a huge
 *   sophisticated C-library for free download with from tiny to 128 bit generators (https://github.com/imneme/pcg-c).
 * - You should seed true 64 bit values or discard the first numbers.
 * - If drawing 64 bit naturally is slow on your machine, try this one.
 * - Facts:
 *   size: 32 bit, period: ~2^64 ~ 1.8e19
 */
//---------- The following code is under the Apache license (do what you want and include license)
//https://www.apache.org/licenses/LICENSE-2.0
uint32_t RngPcg::int32()
{
	//linear congruential state transition function:
	const uint64_t oldstate = state; 
	state = oldstate * 6364136223846793005ULL + (inc|1);
	//permutation function of a tuple as output function:
	const uint32_t xorshifted = (uint32_t)( ((oldstate >> 18u) ^ oldstate) >> 27u );
	const uint32_t rot = (uint32_t)(oldstate >> 59u);
	return  (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
//---------- End Apache license

void RngPcg::getState(std::vector<uint64_t>& ovec_seed) const
{
	ovec_seed.clear();
	ovec_seed.push_back(state);
	ovec_seed.push_back(inc);
}

void RngPcg::setState(const std::vector<uint64_t>& ivec_seed)
{
	if (ivec_seed.size() != 2)
		throw InvalidArgumentException("RNG: Unexpected number of seeds for this generator (needed: 2)", AT);
	
	state = ivec_seed[0];
	inc = ivec_seed[1];
}

/* PRIVATE FUNCTIONS */
bool RngPcg::initAllStates() //initial PCG-generator states
{
	state = 0;
	inc = 0;
	const bool hardware_success_1 = getUniqueSeed(state);
	const bool hardware_success_2 = getUniqueSeed(inc);
	return (hardware_success_1 && hardware_success_2);
}

///////////////////////////////////////////////////////////////////////////////
//    MERSENNE TWISTER GENERATOR class                                       //
///////////////////////////////////////////////////////////////////////////////

/* CONSTRUCTOR */
RngMtw::RngMtw() : MT_NN(624), MT_MM(397), current_mt_index(0), vec_states(std::vector<uint32_t>()) 
{
	/* MT_NN denotes the number of states the MT uses and is not really meant to change, but
	 * doing so doesn't break anything. It is hardcoded for simplicity, standard choices are 314 and 624.
	 * If w is the word size and r the separation point (bits in lower bitmask),
	 * then 2^(NN*w-r)-1 must be a Mersenne Prime.
	 */
	hardware_seed_success = initAllStates();
}

/* PUBLIC FUNCTIONS */
uint64_t RngMtw::int64()
{
	const uint32_t lowpart = int32();
	const uint32_t highpart = int32();
	return RngCore::combine32to64(lowpart, highpart);
}

void RngMtw::getState(std::vector<uint64_t>& ovec_seed) const
{
	ovec_seed.clear();
	ovec_seed.push_back(current_mt_index); //1st element is the currently used index
	for (size_t i = 0; i < vec_states.size(); ++i)
		ovec_seed.push_back( (uint64_t)vec_states[i] );
}

void RngMtw::setState(const std::vector<uint64_t>& ivec_seed)
{
	//assert that we have NN seeds or NN+1 seeds with the 1st one usable as the current index:
	if ( (ivec_seed.size() != MT_NN) && ((ivec_seed.size() != MT_NN + 1) || (ivec_seed[0] >= MT_NN)) ) {
		std::stringstream ss;
		ss << "RNG: Unexpected number of seeds for this generator (needed: " << MT_NN << ")";
		throw InvalidArgumentException(ss.str(), AT);
	}
	unsigned int offset = 0;
	if (ivec_seed.size() == MT_NN + 1) { //initializing from a previous run including the last index
		current_mt_index = (unsigned int)ivec_seed[0];
		offset = 1; //actual states start after the stored index
	} else { //initializing from scratch with NN numbers the user provides
		current_mt_index = 0;
	}

	for (size_t i = 0; i < ivec_seed.size(); ++i)
		vec_states[i] = (uint32_t)ivec_seed[i+offset];
}

/* - doc piece -
 * Random Number Generator 3 - Mersenne Twister
 * - Implementation of the wide-spread Mersenne Twister algorithm by M. Matsumoto and T. Nishimura (Ref. [MN98]).
 * - By using 624 internal states, the period is extremely long.
 * - This does not make it crypto-secure (the state can be derived from 624 random numbers), but it
 *   passes many statistical tests and is the standard RNG in numerous well-known software packages.
 * - It needs a few kB buffer size, which is relatively large compared to the other generators.
 * - Facts:
 *   size: 32 bit, period: 2^19937-1 (Mersenne prime) ~ 4.3e6001
 */

//---------- The following code is adapted from copyrighted but completely free-to-use material by M. Matsumoto and T. Nishimura, Ref. [MN98]
uint32_t RngMtw::int32() //[0, 2^32-1] 
{
	const uint32_t UPPER_MASK = 0x80000000UL; //most significant w-r bits
	const uint32_t LOWER_MASK = 0x7fffffffUL; //least significant r bits
	const uint32_t magic[2] = {0x0UL, 0x9908b0dfUL}; //Twister-matrix

	uint32_t xx(0);

	if (current_mt_index >= MT_NN) { //all NN words were used --> calculate next set
		for (size_t i = 0; i < MT_NN - MT_MM; ++i) { //first bit of current, last 32 bits of next state
			xx = (vec_states[i] & UPPER_MASK) | (vec_states[i+1] & LOWER_MASK);
			vec_states[i] = (uint32_t)( vec_states[i+MT_MM] ^ (xx >> 1) ^ magic[(int)(xx & 0x1UL)] );
		}
		for (size_t i = MT_NN - MT_MM; i < MT_NN-1; ++i) {
			xx = (vec_states[i] & UPPER_MASK) | (vec_states[i+1] & LOWER_MASK);
			vec_states[i] = (uint32_t)( vec_states[i+MT_MM-MT_NN] ^ (xx >> 1) ^ magic[(int)(xx & 0x1UL)] );
		}
		//if odd, xor with magic number
		xx = (vec_states[MT_NN-1] & UPPER_MASK) | (vec_states[0] & LOWER_MASK);
		vec_states[MT_NN-1] = (uint32_t)( vec_states[MT_MM-1] ^ (xx >> 1) ^ magic[(int)(xx & 0x1UL)] );
		current_mt_index = 0; //new set of numbers - can run through them again
	}

	xx = vec_states[current_mt_index];
	current_mt_index++;

	xx ^= (xx >> 11); //output tempering
	xx ^= (xx << 7) & 0x9d2c5680UL;
	xx ^= (xx << 15) & 0xefc60000UL;
	xx ^= (xx >> 18);

	return xx;
}

/* PRIVATE FUNCTIONS */
bool RngMtw::initAllStates() //init all states with a mix of "true" and "pseudo" entropy
{
	//first, the full state array is initialized with values generated from a single seed:
	uint64_t store;
	bool hardware_success = getUniqueSeed(store);
	uint32_t seed = (uint32_t)store; //keeps lower bits
	vec_states.clear();
	vec_states.push_back(seed);
	for (size_t i = 1; i < MT_NN; ++i) { //Ref. [DK81] for multiplier
		vec_states.push_back( (1812433253UL * (vec_states[i-1] ^ (vec_states[i-1] >> 30)) + i) );
		vec_states[i] &= 0xffffffffUL;
	}

	//next, we generate a second seeding array independent of the first one:
	std::vector<uint32_t> seed_states; //this vector could be user input at the cost of another seeding function
	const unsigned int n_additional_seeds = 64; //our choice
	for (size_t i = 0; i < n_additional_seeds; i+=2) {
		const bool hw = getUniqueSeed(store);
		if (!hw) hardware_success = false; //report failed hardware seed
		//decompose the 64 bit hardware seed into two 32 bit states:
		seed_states.push_back( (uint32_t)(store >> 32) ); //high bits
		seed_states.push_back( (uint32_t)(store & 0xffffffffUL) ); //low bits
	}

	//then, the initially generated states are mixed with the additional entropy:
	const uint32_t sz = MT_NN > seed_states.size()? MT_NN : seed_states.size();
	uint32_t i = 1, j = 0;
	for (uint32_t k = sz; k > 0; --k) { //1st step with arbitrarily sized 2nd array
		vec_states[i] = (vec_states[i] ^ ((vec_states[i-1] ^ (vec_states[i-1] >> 30))
		    * 1664525UL)) + seed_states[j] + j; //non-linear
		vec_states[i] &= 0xffffffffUL;
		i++; j++;
		if (i >= MT_NN) {
			vec_states[0] = vec_states[MT_NN-1];
			i = 1;
		}
		if (j>=seed_states.size())
			j = 0;
	}
	for (uint32_t k = MT_NN - 1; k > 0; --k) { //2nd step over all states
		vec_states[i] = (vec_states[i] ^ ((vec_states[i-1] ^ (vec_states[i-1] >> 30))
		    * 1566083941UL)) - i;
		vec_states[i] &= 0xffffffffUL;
		i++;
		if (i >= MT_NN) {
			vec_states[0] = vec_states[MT_NN-1];
			i = 1;
		}
	}
	vec_states[0] = 0x80000000UL; //assure non-zero initial array (most significant bit is 1)
	current_mt_index = MT_NN + 1; //make sure int32() inits on the first run

	return hardware_success;
}
//---------- End algorithm Matsumoto/Nishimura

///////////////////////////////////////////////////////////////////////////////
//    The CORE class                                                         //
///////////////////////////////////////////////////////////////////////////////

/* CONSTRUCTOR */
RngCore::RngCore() : hardware_seed_success(0)
{
	//do nothing
}


RngCore::~RngCore() //we need this declared virtual to be able to delete rng_core
{
	//do nothing
}

/* PUBLIC FUNCTIONS */
bool RngCore::getUniqueSeed(uint64_t& store) const 
{
	#if !defined(_WIN32) && ( defined(__unix__)  || defined(__unix) \
	    || defined(__linux) || (defined(__APPLE__) && defined(__MACH__)) \
	    || defined(ANDROID) )
		uint64_t entropy; //TODO: confirm it's there on Android
		const bool success = getEntropy(entropy);
		if (!success) {
			entropy = timeMixer(time(NULL), clock()); //clock() prevents static
		}
		store = entropy;
		return success;
	#else
		//in the future, there could be a case for the Windows crypto box
		store = timeMixer(time(NULL), clock());
		return false;
	#endif
}

/* PROTECTED FUNCTIONS */
uint64_t RngCore::combine32to64(const uint32_t& low, const uint32_t& high) const 
{
	const uint64_t temp_low = (uint64_t)low;
	const uint64_t temp_high = (uint64_t)high;
	return temp_low | (temp_high << 32); //shift one up 32 bits and or together
}

double RngCore::doubFromInt(const uint64_t& rn) const
{
	//For speed, we don't check boundaries in the standard method (if this leads to problems then probably the generator
	//is broken anyway) and rely on the built in implementation:
	return ldexp((double)rn, -64); //rn * 2^(-64)
	//Return values in [2^-11, 1] are overrepresentated, low precision for small exponents
}

/* -doc piece -
 * Uniform random double values are quite hard to generate. The code example at Ref. [TC14] provides a method to do it,
 * which is to interpret a random stream of bits as fractional part of the binary expansion of a number in [0, 1].
 * The file also goes into details about why other methods are troublesome if we rely on quality, e. g. sensitive random searches
 * on a plane due to the gap size of 1/2^(bits).
 * Another generator capable of producing doubles directly is a lagged Fibonacci generator (quick but maybe poor).
 */
//---------- The following code is provided for free use by Taylor Campbell -- see Ref. [TC14] for full comments
double RngCore::trueDoub() //[0, 1]
{
//#define countLeadingZeros __builtin_clz //safe? That would be nice
	int exponent = -64;
	uint64_t significand;
	unsigned int shift;
	
	while ((significand = int64()) == 0) { //read 0 into exp until 1 is hit, rest goes to significant 
		exponent -= 64;
		if (exponent < -1074) return 0.; //won't realistically happen
	}
	shift = countLeadingZeros(significand);
	if (shift != 0) {
		exponent -= shift; //shift leading zeros into exponent
		significand <<= shift; //refill less-significant bits
		significand |= (int64() >> (64 - shift));
	}
	significand |= 1; //sticky bit
	return ldexp((double)significand, exponent);
}
//---------- End 

/* PRIVATE FUNCTIONS */
bool RngCore::getEntropy(uint64_t& store) const
{ //read 64 bits from hardware entropy device
	//this shouldn't be done too often but for isolated generator seeds it's quite suitable:
	std::ifstream urandom("/dev/urandom", std::ios::in|std::ios::binary); //"should be safe too"
	if (urandom) {
		char *memblock;
		const size_t sz = sizeof(uint64_t);
		memblock = new char[sz];
		urandom.read(memblock, sz);
		store = *reinterpret_cast<uint64_t*>(memblock);
		delete[] memblock;
		urandom.close();
		return true;
	} else {
		return false; //stream failed or depleted
	}
}

uint64_t RngCore::timeMixer(const time_t& tt, const clock_t& cc) const
{ //mix time() and clock() together and hash that to get a 64 bit seed
	static uint32_t diff = 0;
	uint32_t ht = 0;
	unsigned char *pt = (unsigned char*)& tt; //Mersenne Twister-hash
	for (size_t i = 0; i < sizeof(tt); ++i) {
		ht *= std::numeric_limits<unsigned char>::max() + 2U; //avoid UCHAR_MAX
		ht += pt[i];
	}
	uint32_t hc = 0;
	pt = (unsigned char*)& cc;
	for (size_t i = 0; i < sizeof(cc); ++i) {
		hc *= std::numeric_limits<unsigned char>::max() + 2U;
		hc += pt[i];
	}
	ht += (++diff); //ensure different seeds for same time and all running generators
	uint32_t res = 0xca01f9dd * ht - 0x4973f715 * hc; //PCG mixer
	res ^= res >> 16;

	return combine32to64(res, hash(res)); //get 64 bits out of 32 by hashing
}

uint32_t RngCore::hash(const uint32_t& nn) const
{ //fast random hashing function recommended by Ref. [NR3]
	uint64_t v = nn * 3935559000370003845LL + 2691343689449507681LL;
	v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
	v *= 4768777513237032717LL;
	v ^= v << 20; v ^= v >> 41; v ^= v << 5;
	return (v & 0xffffffff); //down to 32 bits
}

size_t RngCore::countLeadingZeros(const uint64_t& nn) const //our own poor man's clz-algorithm
{ //HACK: even without __builtin_clz there are much better ways to do this
	size_t clz = 0;

	const unsigned short bit_char = std::numeric_limits<unsigned char>::digits; //avoid CHAR_BIT
	for (size_t i = 0; i < bit_char * sizeof(nn); ++i) 
	{
		if ((nn & (1 << i)) == 0)
			clz++;
		else //first non-zero character hit
			break;
	}
	return clz;
}

///////////////////////////////////////////////////////////////////////////////
//    OBJECT FACTORY                                                         //
///////////////////////////////////////////////////////////////////////////////

RngCore* RngFactory::getCore(const RandomNumberGenerator::RNG_TYPE& algorithm)
{
	if (algorithm == RandomNumberGenerator::RNG_XOR)
		return new RngXor;
	else if (algorithm == RandomNumberGenerator::RNG_PCG)
		return new RngPcg;
	else if (algorithm == RandomNumberGenerator::RNG_MTW)
		return new RngMtw;
	else
		throw InvalidArgumentException("RNG: This random number generator algorithm is not implemented. Check your RNG_TYPE.", AT);
}

} //namespace
