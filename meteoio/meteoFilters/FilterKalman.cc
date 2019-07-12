#include <LU> //<Eigen/LU>

#include <meteoio/meteoFilters/FilterKalman.h>
#include <meteoio/IOUtils.h>

namespace mio {

FilterKalman::FilterKalman(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name),
		  mat_in_xx(""), mat_in_AA(""), mat_in_HH(""), mat_in_PP("1"), mat_in_QQ("0"), mat_in_RR("0"), mat_in_BB("0"), mat_in_uu("0"),
          meas_params(),
          be_verbose(true)
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;
}


void FilterKalman::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{

	/* INITIALIZATION */

	const size_t TT = ivec.size(); //number of time steps
	const Eigen::VectorXd vecT = buildTimeVector(ivec);

	//TODO: allow input with brackets (no effect)
	Eigen::MatrixXd xx = parseMatrix(mat_in_xx, IOUtils::unodata, 1); //this ini line fixes our number of states
	const size_t nx = xx.size();

	//knowing the number of states, we can init state matrix sizes
	std::vector<std::string> AA_str = parseSystemMatrix(mat_in_AA, nx);
	//TODO: adaptive time stepping

	Eigen::MatrixXd zz;
	std::vector<size_t> meas_idx;
	const size_t nz = buildObservationsMatrix(param, ivec, zz, meas_idx); //check observations dimension

	//now that we know how many states and observables there are, we can size the rest of the matrices:
	Eigen::MatrixXd HH(nz, nx);
	if (!mat_in_HH.empty()) { //use identity matrix or something similar if the matrix isn't square
		HH.setIdentity();
		if (be_verbose) std::cerr << "[W] No Kalman filter model found to relate observations to states. Setting to 'identity', but this may not make much sense.\n";
	} else {
		HH = parseMatrix(mat_in_HH, nz, nx);
	}

	Eigen::MatrixXd PP = bloatMatrix(mat_in_PP, nx, nx); //trust in initial state
	Eigen::MatrixXd QQ = bloatMatrix(mat_in_QQ, nx, nx); //process noise covariance
	Eigen::MatrixXd RR = bloatMatrix(mat_in_RR, nz, nz); //observation noise covariance

	//at last, we input an optional control signal, either as scalar, matrix, or "meteo" data
	Eigen::MatrixXd BB = bloatMatrix(mat_in_BB, nx, nx); //relates control input to state
	Eigen::MatrixXd uu = buildControlSignal(nx, TT, ivec);

    /* KALMAN FILTER */

    Eigen::MatrixXd KK; //Kalman gain
    Eigen::MatrixXd II = Eigen::MatrixXd::Identity(nx, nx);

    ovec = ivec; //copy with all special parameters etc.

    bool saw_nodata(false);
    for (size_t kk = 0; kk < TT; ++kk) //for each time step...
    {

    	const double dt = (kk == 0)? 0 : vecT(kk) - vecT(kk-1); //initial state is at time of 1st measurement
    	Eigen::MatrixXd AA = buildSystemMatrix(AA_str, nx, dt);

    	const bool has_nodata = checkNodata(zz.col(kk));
    	if (!has_nodata) { //we can calculate
    		//prediction
    		xx = ( AA * xx + BB * uu.col(kk) ).eval(); //guard against aliasing
    		PP = ( AA * PP * AA.transpose() + QQ ).eval();

    		//update
    		KK = PP * HH.transpose() * (HH * PP * HH.transpose() + RR).inverse();
    		xx = ( xx + KK * (zz.col(kk) - HH * xx) ).eval();
    		PP = ( (II - KK * HH) * PP ).eval();

    		//for (size_t ii = 0; ii < nz; ++ii)
    			//ovec[kk](meas_idx[ii]) = xx(ii); //Note! This will filter all specified parameters, not just param!
    		ovec[kk](param) = xx(0); //TODO: correlate states with observables

    		//last_nodata_idx = ... for dt

    	} else { //all states stay as if this time step wasn't encountered
    		saw_nodata = true;
    		for (size_t ii = 0; ii < nz; ++ii)
    			ovec[kk](meas_idx[ii]) = ivec[kk](meas_idx[ii]); //transport all parameters as-is
    	} //endif nodata


    } // endfor kk

    if (be_verbose && saw_nodata) std::cerr << "[W] Nodata value(s) encountered in Kalman filter; some values were ignored. You should probably resample beforehand.\n";
}

size_t FilterKalman::buildObservationsMatrix(const unsigned int& param, const std::vector<MeteoData>& ivec, Eigen::MatrixXd& zz,
		std::vector<size_t>& meas_idx) const
{
	meas_idx.clear(); //index map of observables
	meas_idx.push_back(param);
	for (size_t jj = 0; jj < meas_params.size(); ++jj)
		meas_idx.push_back( ivec.front().getParameterIndex(meas_params[jj]) );

	zz.resize(meas_idx.size(), ivec.size()); //fill observation matrix with desired values from meteo set
	for (size_t kk = 0; kk < ivec.size(); ++kk) {
		for (size_t jj = 0; jj < meas_idx.size(); ++jj)
			zz(jj, kk) = ivec[kk](meas_idx[jj]);
	}
	return meas_idx.size();
}

Eigen::MatrixXd FilterKalman::buildControlSignal(const size_t& nx, const size_t& TT, const std::vector<MeteoData>& ivec) const
{
	Eigen::MatrixXd uu(nx, TT);

	std::vector<std::string> vecU;
	const size_t nr_uu = IOUtils::readLineToVec(mat_in_uu, vecU, ',');

	double aa;
	std::istringstream iss(vecU.front());
	iss >> aa; //test-read double

	if (!iss.fail()) { //there's at least a numerical value in the beginning, so it should not be a meteo param name
		if (nr_uu == 1) { //case 1: single value --> replicate to whole matrix
			uu.setOnes();
			uu *= aa;
		} else if (nr_uu == nx) { //case 2: nx values --> repeat column vector
			for (size_t jj = 0; jj < nx; ++jj) {
				iss.clear(); iss.str(vecU[jj]);
				iss >> aa;
				uu(jj, 0) = aa;
			}
			uu.colwise() = uu.col(0);
		} else {
			throw InvalidArgumentException("Control signal vector for Kalman filter ill-formatted (expected 1 or " +
			        IOUtils::toString(nx) + " elements).", AT);
		}
	} else { //case 3: parameter names given --> separate vector each timestep
		if (nr_uu != nx)
			throw InvalidArgumentException("Control signal vector for Kalman filter cannot be constructed (expected: " +
			        IOUtils::toString(nx) + " parameter names).", AT);
		for (size_t jj = 0; jj < TT; ++jj)
			for (size_t ii = 0; ii < nx; ++ii) {
				IOUtils::trim(vecU[ii]);
				uu(ii, jj) =  ivec[jj](vecU[ii]);
			}
	}

	return uu;
}

Eigen::MatrixXd FilterKalman::parseMatrix(const std::string& line, const size_t& rows, const size_t& cols) const
{
	std::vector<double> vecRet; //handle growing vector with stl
	const size_t nr_elements = IOUtils::readLineToVec(line, vecRet, ',');

	size_t xrows(rows), xcols(cols);
	if (rows == IOUtils::unodata && cols == 1) { //row or column vector - unambiguous size from input line
		xrows = nr_elements;
	} else if (cols == IOUtils::unodata && rows == 1) {
		xcols = nr_elements;
	} else {
		if (rows*cols != nr_elements)
			throw InvalidArgumentException("The Kalman filter encountered an unfit input matrix size (expected: " +
			        IOUtils::toString(rows) + "x" + IOUtils::toString(cols) + ").", AT);
	}

	//we have to copy to be able to let vecRet go out of scope, but we don't need to loop:
	const Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ET(&vecRet.data()[0], xrows, xcols);
	return ET; //RowMajor is the storage order and has the effect of row-wise input in the ini
}

std::vector<std::string> FilterKalman::parseSystemMatrix(const std::string& line, const size_t& sz) const
{
	std::vector<std::string> vecRet; //keep as strings to be able to do substitutions
	const size_t nr_elements = IOUtils::readLineToVec(line, vecRet, ',');

	if (sz*sz != nr_elements) //it's a square matrix
		throw InvalidArgumentException("The Kalman filter encountered an unfit state input matrix size (expected: " +
		        IOUtils::toString(sz) + "x" + IOUtils::toString(sz) + ").", AT);

	return vecRet;
}

Eigen::MatrixXd FilterKalman::buildSystemMatrix(const std::vector<std::string>& AA_str, const size_t& sz, const double& dt) const
{ //substitutions in system matrix
	Eigen::MatrixXd AA(sz, sz);
	for (size_t ii = 0; ii < sz; ++ii) {
		for (size_t jj = 0; jj < sz; ++jj) {
			AA(ii, jj) = sub(AA_str[ii*sz + jj], dt);
		}
	}
	return AA;
}

double FilterKalman::sub(const std::string& expr, const double& dt) const
{
	if (expr == "dt") { //time step
		return dt;
	} else {
		std::istringstream ss(expr);
		double aa;
		ss >> aa;
		if (ss.fail())
			throw InvalidArgumentException("Unrecognized value in the Kalman filter's system matrix: \"" + expr + "\".", AT);
		return aa;
	}
}

Eigen::MatrixXd FilterKalman::bloatMatrix(const std::string& line, const size_t& rows, const size_t& cols) const
{ //if a scalar is given put it on the diagonal, else read the matrix as usual
	Eigen::MatrixXd ret;
	const size_t nr = IOUtils::count(line, ",");
	if (nr > 0) { //not an isolated double
		ret = parseMatrix(line, rows, cols);
	} else { //single value: put on diagonal
		std::istringstream ss(line);
		double aa;
		ss >> aa;
		ret.setIdentity(rows, cols);
		ret *= aa;
	}
	return ret;
}

Eigen::VectorXd FilterKalman::buildTimeVector(const std::vector<MeteoData>& ivec) const
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

	return vecRet; //TODO: save last time step to disk to be able to resume?
}

bool FilterKalman::checkNodata(const Eigen::VectorXd& ivec) const
{ //is any vector element nodata?
	for (int ii = 0; ii < ivec.size(); ++ii)
		if (ivec(ii) == IOUtils::nodata)
			return true;
	return false;
}

void FilterKalman::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where("Filters::" + block_name);

	bool has_initial(false), has_dynamics(false);

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {
		if (vecArgs[ii].first == "INITIAL_STATE") {
			mat_in_xx = vecArgs[ii].second;
			has_initial = true;
		} else if (vecArgs[ii].first == "INITIAL_TRUST") {
			mat_in_PP = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "STATE_DYNAMICS") {
			mat_in_AA = vecArgs[ii].second;
			has_dynamics = true;
		} else if (vecArgs[ii].first == "PROCESS_COVARIANCE") {
			mat_in_QQ = vecArgs[ii].second;

		}

		else if (vecArgs[ii].first == "ADD_OBSERVATIONS") {
			(void) IOUtils::readLineToVec(vecArgs[ii].second, meas_params);
		} else if (vecArgs[ii].first == "OBSERVATION_RELATION") {
			mat_in_HH = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "OBSERVATION_COVARIANCE") { //TODO: make mandatory
			mat_in_RR = vecArgs[ii].second;
		}

		else if (vecArgs[ii].first == "CONTROL_SIGNAL") {
			mat_in_uu = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "CONTROL_RELATION") {
			mat_in_BB = vecArgs[ii].second;
		}

		else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);
		}
	} //endfor vecArgs

	if (!has_initial) //TODO: 1st measurement?
		throw InvalidArgumentException("The Kalman filter needs an initial state, i. e. a vector given in INITIAL_STATE.", AT);
	if (!has_dynamics)
		throw InvalidArgumentException("The Kalman filter needs a linear model, i. e. a matrix given in STATE_DYNAMICS.", AT);

}

} //namespace