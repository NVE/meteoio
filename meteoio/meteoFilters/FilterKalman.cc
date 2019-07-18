#include <LU> //<Eigen/LU>

#include <meteoio/meteoFilters/FilterKalman.h>
#include <meteoio/IOUtils.h>

namespace mio {

FilterKalman::FilterKalman(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
        : ProcessingBlock(vecArgs, name),
		  mat_in_xx(""), mat_in_AA(""), mat_in_HH(""), mat_in_PP("1"), mat_in_QQ("0"), mat_in_RR(""), mat_in_BB("0"), mat_in_uu("0"),
          meas_params(), error_params(), QQ_params(), RR_params(),
          be_verbose(true), unrecognized_key("")
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first;

	cleanBrackets(mat_in_xx); cleanBrackets(mat_in_AA); cleanBrackets(mat_in_HH); cleanBrackets(mat_in_PP); //remove all brackets
	cleanBrackets(mat_in_QQ); cleanBrackets(mat_in_RR); cleanBrackets(mat_in_BB); cleanBrackets(mat_in_uu); //(they are just visual help)
}


void FilterKalman::process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec)
{

	/* INITIALIZATION */

	if (!unrecognized_key.empty() && be_verbose) //because we have so many input keys
		std::cerr << "[W] Unrecognized ini key(s) ignored for Kalman filter, one of them is: \"" + unrecognized_key + "\"\n";

	const size_t TT = ivec.size(); //number of time steps
	const Eigen::VectorXd vecT = buildTimeVector(ivec);

	std::vector<std::string> xx_str;
	size_t nx = IOUtils::readLineToVec(mat_in_xx, xx_str, ','); //fix number of internal states

	Eigen::MatrixXd zz;
	std::vector<size_t> meas_idx; //indices of given observables in meteo set
	const size_t nz = buildObservationsMatrix(param, ivec, zz, meas_idx, nx); //checks observations dimension - must match nx

	if (nx == 0) //special case: no initial state given
		nx = nz; //--> number of observations fixes the states

	Eigen::MatrixXd xx = buildInitialStates(xx_str, meas_idx, ivec, nz); //once the observables are known, go back to parsing the initial states
	std::vector<std::string> AA_str = parseSystemMatrix(mat_in_AA, nx); //keep as strings so we can do substitutions

	//Now that we know how many states and observables there are, we can size the rest of the matrices.
	//We keep nz for clarity, but this implementation forces nz=nx (but not all have to be actually used).
	Eigen::MatrixXd HH(nz, nx);
	if (mat_in_HH.empty()) { //use identity matrix or something similar if the matrix isn't square
		HH.setIdentity();
		if (be_verbose) std::cerr << "[W] No Kalman filter model found to relate observations to states. Setting to \"identity\", but this may not make much sense.\n";
	} else {
		HH = bloatMatrix(mat_in_HH, nz, nx, "observation model");
	}

	Eigen::MatrixXd PP = bloatMatrix(mat_in_PP, nx, nx, "initial trust"); //trust in initial state
	Eigen::MatrixXd QQ = bloatMatrix(mat_in_QQ, nx, nx, "process noise"); //process noise covariance
	Eigen::MatrixXd RR = bloatMatrix(mat_in_RR, nz, nz, "observation noise"); //observation noise covariance
	bool has_RR_params, has_QQ_params;
	assertInputCovariances(ivec, has_RR_params, has_QQ_params, nx);

	//at last, we input an optional control signal, either as scalar, matrix, or "meteo" data:
	Eigen::MatrixXd BB = bloatMatrix(mat_in_BB, nx, nx, "control relation"); //relates control input to state
	Eigen::MatrixXd uu = buildControlSignal(nx, TT, ivec);

    /* KALMAN FILTER */

    Eigen::MatrixXd KK; //Kalman gain
    Eigen::MatrixXd II = Eigen::MatrixXd::Identity(nx, nx);

    ovec = ivec; //copy with all special parameters etc.


    for (size_t pp = 0; pp < error_params.size(); ++pp) //create parameters to output PP if they don't exist
    	for (size_t ii = 0; ii < ovec.size(); ++ii)
    		ovec[ii].addParameter(error_params[pp]);

    bool saw_nodata(false);
    size_t last_valid_idx(0);
    for (size_t kk = 0; kk < TT; ++kk) //for each time step...
    {
    	const double dt = (kk == 0)? 0 : vecT(kk) - vecT(last_valid_idx); //initial state is at time of 1st measurement
    	const Eigen::MatrixXd AA = buildSystemMatrix(AA_str, nx, dt, ivec, kk);
    	if (has_RR_params) //if desired, get RR and QQ from input data at each time step
    		RR = extractInputMatrix(RR_params, ivec, kk);
    	if (has_QQ_params)
    		QQ = extractInputMatrix(QQ_params, ivec, kk);
    	const bool has_nodata = checkNodata(zz.col(kk));

    	if (!has_nodata) { //we can calculate

    		//prediction:
    		xx = ( AA * xx + BB * uu.col(kk) ).eval(); //guard against aliasing
    		PP = ( AA * PP * AA.transpose() + QQ ).eval();

    		//update:
    		KK = PP * HH.transpose() * (HH * PP * HH.transpose() + RR).inverse();
    		xx = ( xx + KK * (zz.col(kk) - HH * xx) ).eval();
    		PP = ( (II - KK * HH) * PP ).eval();

    		for (size_t ii = 0; ii < nz; ++ii) //nx = nz! The user is responsible to provide enough output fields.
    			ovec[kk](meas_idx[ii]) = xx(ii); //Note! This will filter all specified parameters, not just param!

    		last_valid_idx = kk; //dt is calculated between valid data values

    		if (error_params.size() == nx) //output estimated error
    			for (size_t mm = 0; mm < nx; ++mm)
    				ovec[kk](error_params[mm]) = PP(mm, mm); //save diagonal elements only

    	} else { //all states stay as if this time step wasn't encountered
    		saw_nodata = true;
    	} //endif nodata
    } // endfor kk

    if (be_verbose && saw_nodata) std::cerr << "[W] Nodata value(s) or missing parameter encountered in Kalman filter; some values were ignored. You should probably resample beforehand.\n";
}

Eigen::VectorXd FilterKalman::buildInitialStates(const std::vector<std::string>& xx_str_in, std::vector<size_t>& meas_idx,
        const std::vector<MeteoData>& ivec, const size_t& nz)
{
	bool saw_nodata(false);
	Eigen::VectorXd vecRet(nz);
	if (xx_str_in.empty()) { //rely on number of observations given
		for (size_t jj = 0; jj < nz; ++jj) { //guaranteed to be at least 1 - the one the filter runs on
			vecRet[jj] = ivec.front()(meas_idx[jj]); //initial_state[t] = observation[t]
			if (vecRet[jj] == IOUtils::nodata)
				saw_nodata = true;
		}
	} else { //not empty
		for (size_t ii = 0; ii < nz; ++ii) {
			const std::string xx_str = IOUtils::trim(xx_str_in[ii]);
			if (xx_str == "") { //some values are given, but not this one --> fill with 1st observation
				vecRet[ii] = ivec.front()(meas_idx[ii]);
				if (vecRet[ii] == IOUtils::nodata)
					saw_nodata = true;
			} else if (xx_str == "average") { //take the average of all specified observables at t=0
				double av_sum = 0.;
				for (size_t jj = 0; jj < nz; ++jj) {
					av_sum += ivec.front()(meas_idx[jj]);
					if (ivec.front()(meas_idx[jj]) == IOUtils::nodata)
						saw_nodata = true;
				}
				vecRet[ii] = av_sum / (double)nz;
			} else if ( xx_str.size() > 6 && (xx_str.compare(0, 6, "meteo(") == 0) ) { //meteo parameters, e. g. "meteo(RH)"
				vecRet[ii] = ivec.front()(xx_str.substr(6, xx_str.length()-7)); //e. g. "RH"
				if (vecRet[ii] == IOUtils::nodata)
					saw_nodata = true;
			} else { //double value expected if we arrive here
				std::istringstream iss(xx_str);
				iss >> vecRet[ii];
				if (iss.fail())
					throw InvalidArgumentException("The Kalman filter could not read initial state \"" + xx_str	+ "\".", AT);
			}
		} //endfor ii
	} //endif empty

	if (saw_nodata && be_verbose) std::cerr << "[W] Nodata values found (and kept) in the Kalman filter's automatic state initialization!\n";
	return vecRet;
}

size_t FilterKalman::buildObservationsMatrix(const unsigned int& param, const std::vector<MeteoData>& ivec, Eigen::MatrixXd& zz,
		std::vector<size_t>& meas_idx, const size_t& nx) const
{ //extract observables from the meteo set to a dedicated vector
	meas_idx.clear(); //index map of observables
	meas_idx.push_back(param); //1st one is always the one the filter runs on
	for (size_t jj = 0; jj < meas_params.size(); ++jj) //then come the ones specified in ADD_OBSERVABLES
		meas_idx.push_back( ivec.front().getParameterIndex(meas_params[jj]) );

	if ( (meas_idx.size() != nx) && (nx != 0) ) //nx = 0: pick input states --> no warning
		throw InvalidArgumentException("You need to provide as many observation variables as you have states variables for the Kalman filter ("
		        + IOUtils::toString(nx) + " including the one the filter runs on).", AT);

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
				iss.clear();
				iss.str(vecU[jj]);
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

Eigen::MatrixXd FilterKalman::parseMatrix(const std::string& line, const size_t& rows, const size_t& cols,
        const std::string& blockname) const
{ //matrix parsing where no substitutions are done - i. e. double values only
	std::vector<double> vecRet; //handle growing vector with stl
	const size_t nr_elements = IOUtils::readLineToVec(line, vecRet, ',');

	if (rows*cols != nr_elements)
		throw InvalidArgumentException("The Kalman filter encountered an unfit " + blockname + " input matrix size (expected: " +
		        IOUtils::toString(rows) + "x" + IOUtils::toString(cols) + ").", AT);

	//we have to copy to be able to let vecRet go out of scope, but we don't need to loop:
	const Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ET(&vecRet.data()[0], rows, cols);
	return ET; //RowMajor is the storage order and has the effect of row-wise input in the ini
}

std::vector<std::string> FilterKalman::parseSystemMatrix(const std::string& line, const size_t& sz) const
{
	std::vector<std::string> vecRet; //keep as strings to be able to do substitutions
	const size_t nr_elements = IOUtils::readLineToVec(line, vecRet, ',');

	if (nr_elements == 1) { //single value --> put on diagonal, rest is 0
		vecRet.resize(sz*sz);
		for (size_t ii = 0; ii < sz; ++ii)
			for (size_t jj = 0; jj < sz; ++jj)
				vecRet[ii*sz + jj] = (ii == jj)? vecRet.front() : "0";
	} else if (nr_elements != sz*sz) { //it's a square matrix
		throw InvalidArgumentException("The Kalman filter encountered an unfit state input matrix size (expected: " +
		        IOUtils::toString(sz) + "x" + IOUtils::toString(sz) + ").", AT);
	}
	return vecRet;
}

Eigen::MatrixXd FilterKalman::buildSystemMatrix(const std::vector<std::string>& AA_str, const size_t& sz, const double& dt,
        const std::vector<MeteoData>& ivec, const size_t& kk) const
{ //substitutions in system matrix
	Eigen::MatrixXd AA(sz, sz);
	for (size_t ii = 0; ii < sz; ++ii)
		for (size_t jj = 0; jj < sz; ++jj)
			AA(ii, jj) = substitute(AA_str[ii*sz + jj], dt, ivec, kk);
	return AA;
}

Eigen::MatrixXd FilterKalman::extractInputMatrix(const std::vector<std::string>& vecParams, const std::vector<MeteoData>& mvec,
        const size_t& kk) const
{ //read meteo parameters (given in a vector) to a diagnoal matrix
	Eigen::VectorXd dia(vecParams.size());
	for (size_t ii = 0; ii < vecParams.size(); ++ii)
		dia(ii) = mvec[kk](vecParams[ii]);
	return dia.asDiagonal();
}

void FilterKalman::assertInputCovariances(const std::vector<MeteoData>& ivec, bool& has_RR_params, bool& has_QQ_params,
        const size_t& nx) const
{ //check if all given parameters to input RR and QQ exist
	if (RR_params.size() == nx) {
		has_RR_params = true;
		for (size_t rr = 0; rr < RR_params.size(); ++rr) {
			if ( ivec.front().getParameterIndex(RR_params[rr]) == IOUtils::npos ) { //check 1st one only...
				has_RR_params = false;
				break;
			}
		}
	} else  {
		has_RR_params = false;
	}
	if (QQ_params.size() == nx) {
		has_QQ_params = true;
		for (size_t qq = 0; qq < QQ_params.size(); ++qq) {
			if ( ivec.front().getParameterIndex(QQ_params[qq]) == IOUtils::npos ) {
				has_QQ_params = false;
				break;
			}
		}
	} else  {
		has_QQ_params = false;
	}
	if (mat_in_RR.empty() && !has_RR_params) //no data here, and no fixed input matrix
		throw InvalidArgumentException("The Kalman filter needs an OBSERVATION_COVARIANCE matrix. No valid parameters found in meteo data either.", AT);
	if (has_RR_params && !mat_in_RR.empty())
		if (be_verbose) std::cerr << "[W] The Kalman filter's observation covariance matrix is provided as single matrix and in the meteo data. Using the latter.\n";
	if (has_QQ_params && !mat_in_QQ.empty())
		if (be_verbose) std::cerr << "[W] The Kalman filter's process covariance matrix is provided as single matrix and in the meteo data. Using the latter.\n";
}

double FilterKalman::substitute(const std::string& expr, const double& dt, const std::vector<MeteoData>& ivec, const size_t& kk) const
{
	const std::string texp = IOUtils::trim(expr);
	if (texp == "dt") { //current time step (time between measurements)
		return dt;
	} else if ( texp.size() > 6 && (texp.compare(0, 6, "meteo(") == 0) ) { //meteo parameters
		return ivec[kk](texp.substr(6, texp.length()-7));
	} else { //double values
		std::istringstream ss(texp);
		double aa;
		ss >> aa;
		if (ss.fail())
			throw InvalidArgumentException("Unrecognized value in the Kalman filter's system matrix: \"" + texp + "\".", AT);
		return aa;
	}
}

Eigen::MatrixXd FilterKalman::bloatMatrix(const std::string& line, const size_t& rows, const size_t& cols, const std::string& blockname) const
{ //if a scalar is given put it on the diagonal, same with a vector, else read the complete matrix
	Eigen::MatrixXd ret;
	if (line.length() == 0)
		return ret;
	const size_t nr = IOUtils::count(line, ",");
	if ( (nr == rows - 1) && (rows == cols) ) { //parse as row vector and put that on diagonal
		Eigen::MatrixXd tmp = parseMatrix(line, rows, 1, blockname);
		ret = tmp.asDiagonal();
	} else if (nr > 0) { //not an isolated double
		ret = parseMatrix(line, rows, cols, blockname);
	} else { //single value: put on diagonal
		std::istringstream ss(line);
		double aa;
		ss >> aa;
		if (ss.fail())
			throw InvalidArgumentException("Kalman filter could not parse value \"" + line + "\" for " + blockname + " matrix.", AT);
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
	for (size_t ii = 0; ii < ivec.size(); ++ii)
		vecRet(ii) = (ivec[ii].date.getJulian() - base_time) / dt;

	return vecRet;
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

	for (size_t ii = 0; ii < vecArgs.size(); ii++) {
		if (vecArgs[ii].first == "INITIAL_STATE") {
			mat_in_xx = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "INITIAL_TRUST") {
			mat_in_PP = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "STATE_DYNAMICS") {
			mat_in_AA = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "PROCESS_COVARIANCE") {
			mat_in_QQ = vecArgs[ii].second;

		}

		else if (vecArgs[ii].first == "ADD_OBSERVABLES") {
			(void) IOUtils::readLineToVec(vecArgs[ii].second, meas_params);
		} else if (vecArgs[ii].first == "OBSERVATION_RELATION") {
			mat_in_HH = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "OBSERVATION_COVARIANCE") {
			mat_in_RR = vecArgs[ii].second;
		}

		else if (vecArgs[ii].first == "CONTROL_SIGNAL") {
			mat_in_uu = vecArgs[ii].second;
		} else if (vecArgs[ii].first == "CONTROL_RELATION") {
			mat_in_BB = vecArgs[ii].second;
		}

		else if (vecArgs[ii].first == "VERBOSE") {
			IOUtils::parseArg(vecArgs[ii], where, be_verbose);
		} else if (vecArgs[ii].first == "OUT_ESTIMATED_ERROR") {
			(void) IOUtils::readLineToVec(vecArgs[ii].second, error_params);
		} else if (vecArgs[ii].first == "PROCESS_COVARIANCE_PARAMS") {
			(void) IOUtils::readLineToVec(vecArgs[ii].second, QQ_params);
		} else if (vecArgs[ii].first == "OBSERVATION_COVARIANCE_PARAMS") {
			(void) IOUtils::readLineToVec(vecArgs[ii].second, RR_params);
		} else {
			unrecognized_key = vecArgs[ii].first; //constructor is always called twice - show only when processing
		}
	} //endfor vecArgs

	if (mat_in_AA.empty())
		throw InvalidArgumentException("The Kalman filter needs a linear model, i. e. a matrix given in STATE_DYNAMICS.", AT);
}

void FilterKalman::cleanBrackets(std::string& iline)
{ //allow input with brackets to no effekt
	IOUtils::replace_all(iline, "][", ", ");
	IOUtils::replace_all(iline, "] [", ", ");
	IOUtils::replace_all(iline, "[", "");
	IOUtils::replace_all(iline, "]", "");
}

} //namespace
