#include <meteoio/meteoFilters/FilterKalman.h>
#include <meteoio/IOUtils.h>


#include <sstream> //for readLineToVec

namespace mio {

FilterKalman::FilterKalman(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name), matrix_input(), meas_params(),
          be_verbose(true)
{

	matrix_input[0] = matrix_input[1] = matrix_input[2] = "";

	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}


void FilterKalman::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{

	//const size_t TT = ivec.size(); //number of time steps

	//TODO: allow input with brackets (no effect)
	Eigen::MatrixXd xx = parseMatrix(matrix_input[0], IOUtils::npos, 1); //this ini line fixes our number of states
	const size_t nx = xx.size();

	//knowing the number of states, we can init state matrix sizes
	Eigen::MatrixXd AA = parseMatrix(matrix_input[1], nx, nx);

	Eigen::MatrixXd zz; //check observations dimension
	const size_t nz = buildObservationsMatrix(param, ivec, zz);

	//now that we know how many states and observables there are, we can size the rest of the matrices:
	Eigen::MatrixXd HH(nz, nx);
	if (matrix_input[2].length() == 0) { //use identity matrix
		HH.setIdentity();
		if (be_verbose) std::cerr << "[W] No model found to relate observations to states. Setting to 'identity', but this may not make much sense.\n";
	} else {
		HH = parseMatrix(matrix_input[2], nz, nx);
	}

	ovec = ivec;
}

Eigen::MatrixXd FilterKalman::parseMatrix(const std::string& line, const size_t& rows, const size_t& cols)
{
	std::vector<double> vecRet;
	const size_t nr_elements = IOUtils::readLineToVec(line, vecRet, ',');

	size_t xrows(rows), xcols(cols); //row or column vector - unambiguous size from input line
	if (rows == IOUtils::unodata && cols == 1) {
		xrows = nr_elements;
	} else if (cols == IOUtils::unodata && rows == 1) {
		xcols = nr_elements;
	} else {
		if (rows*cols != nr_elements)
			throw InvalidArgumentException("The Kalman filter encountered an unfit matrix size (expected: " +
			        IOUtils::toString(rows) + "x" + IOUtils::toString(cols) + ").", AT);
	}

	//we have to copy to be able to let vecRet go out of scope, but we don't need to loop:
	const Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ET(&vecRet.data()[0], xrows, xcols);
	return ET; //RowMajor is the storage order and has the effect of row-wise input in the ini
}

size_t FilterKalman::buildObservationsMatrix(const unsigned int& param, const std::vector<MeteoData>& ivec, Eigen::MatrixXd& zz) {
	std::vector<size_t> meas_idx; //index map of observables
	meas_idx.push_back(param);
	for (size_t j = 0; j < meas_params.size(); ++j)
		meas_idx.push_back( ivec.front().getParameterIndex(meas_params[j]) );

	zz.resize(meas_idx.size(), ivec.size()); //fill observation matrix with desired values from meteo set
	for (size_t kk = 0; kk < ivec.size(); ++kk) {
		for (size_t jj = 0; jj < meas_idx.size(); ++jj)
			zz(jj, kk) = ivec[kk](meas_idx[jj]);
	}

	return meas_idx.size();
}

void FilterKalman::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where("Filters::" + block_name);

	bool has_initial(false), has_dynamics(false);

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {
		if (vecArgs[ii].first == "INITIAL_STATE") {
			matrix_input[0] = vecArgs[ii].second;
			has_initial = true;
		} else if (vecArgs[ii].first == "STATE_DYNAMICS") {
			matrix_input[1] = vecArgs[ii].second;
			has_dynamics = true;
		} else if (vecArgs[ii].first == "ADD_OBSERVATIONS") {
			(void) IOUtils::readLineToVec(vecArgs[ii].second, meas_params);
		} else if (vecArgs[ii].first == "OBSERVATION_RELATION") {
			matrix_input[2] = vecArgs[ii].second;
		}

		else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);
		}
	}


	if (!has_initial)
		throw InvalidArgumentException("The Kalman filter needs an initial state, i. e. a vector given in INITIAL_STATE.", AT);
	if (!has_dynamics)
		throw InvalidArgumentException("The Kalman filter needs a linear model, i. e. a matrix given in STATE_DYNAMICS.", AT);

}

} //namespace
