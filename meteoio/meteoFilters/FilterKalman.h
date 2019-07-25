#ifndef FILTERKALMAN_H
#define FILTERKALMAN_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoStats/RandomNumberGenerator.h>

#include <Core> //<Eigen/Core>

#include <string>
#include <vector>

namespace mio {

/**
 * @class FilterKalman
 * @ingroup processing
 * @author Michael Reisecker
 * @date   2019-07
 * @brief A statistical filter for state likelihood estimation: the <a href="https://en.wikipedia.org/wiki/Kalman_filter">Kalman filter</a>.
 * @details
 *
 * @tableofcontents
 *
 * \ref introduction <br>
 * \ref kalmanfilter : \ref kalmanoverview, \ref kalmanalgorithm, \ref kalmanexample <br>
 * \ref kalmanfeatures : \ref kalmancontrol, \ref kalmanerror, \ref kalmanvariances, \ref kalmannotes <br>
 * \ref kalmankeylist
 *
 * @section introduction Introduction
 *
 * @note After you have read this, some concepts are further explained in FilterParticle; references are also found there.
 * The documentation is meant to be read in conjunction.
 * @note It is advised that you have a TeX installation in place to compile this part of the documentation yourself.
 *
 * This filter suite implements a particle and a Kalman Filter.

 * __Disclaimer__: The Kalman and particle filters are <b>not</b> meant for operational use, rather they aim to be an accessible framework to
 * scientific research backed by the convenience of MeteoIO, to quickly try a few things.
 *
 * What do they do?
 *
 * Both filters follow a statistical approach. The particle filter is a _classical Monte Carlo method_ and samples states from a
 * distribution function. The Kalman filter solves matrix equations of _Bayesian statistics_ theory.
 * They solve the problem of combining model data with noisy measurements in an optimal way.

 * Who are they for?
 *
 * a) Suppose you have your toolchain set up with MeteoIO and are interested in Bayesian filters. With this software, you can
 * (relatively) easily start to explore in this direction and hopefully decide if it is an idea worth pursuing.
 *
 * b) You try a filter, play with the settings, and for whatever reason it works.
 *
 * The difference to other, easier to use, filters is that you, the user, must provide the Physics in form of analytical models and
 * covariance matrices. Almost all will be defaulted, but this may not make much sense in your use case.
 *
 * @section kalmanfilter The Kalman filter
 *
 * A Kalman filter works on a _linear state model_ which propagates an initial state forward in time. At each time step, the filter will
 * adjust this model taking into account _incoming measurements_ tainted with white (Gaussian) noise. For this, we can tune how much trust
 * we put in our measurements (initial and online), and thus how strongly the filter will react to observations.
 * @note In this scope, "online" means that new measurements can be incorporated on the fly as they come in without having to recalculate
 * the past.
 *
 * @subsection kalmanoverview Overview
 *
 * You need your <b>linear model function</b> in matrix form, initial values at \f$T=0\f$ and a second <b>model that relates
 * observations</b> to the internal states. If the model describes the measured parameters directly, this is the identity matrix.
 * You supply an <b>initial trust</b> in the initial states, and this matrix will be propagated as error estimation.
 * Furthermore, you assign a level of <b>trust in the model function</b>, and a level of <b>trust in the measurements</b>.
 *
 * @subsection kalmanalgorithm Algorithm
 *
 * The Kalman filter solves the following set of recursive equations:
 *
 * _Prediction_:
 * \f[
 * \hat{x}_k=A\hat{x}_{k-1} + B u_k
 * \f]
 * \f[
 * P_k = A P_{k-1}A^T+Q
 * \f]
 * _Update_:
 * \f[
 * K_k = P_k H^T(H P_k H^T + R)^{-1}
 * \f]
 * \f[
 * \hat{x}_k = \hat{x}_k + K_k(z_k - H\hat{x}_k)
 * \f]
 * \f[
 * P_k = (I - K_k H)P_k
 * \f]
 *
 * where \f$\hat{x}_k\f$ is the model value at time step \f$k\f$, \f$A\f$ is the model itself, i. e. the system transformation matrix,
 * \f$u_k\f$ is a control signal that is related to the states via \f$B\f$, \f$P\f$ is an error estimation, and \f$Q\f$ is the covariance
 * matrix of the process noise. Furthermore, \f$K\f$ denotes the Kalman gain, \f$H\f$ is the observation model relating observations to
 * the states, \f$R\f$ is the observation covariance matrix, and \f$I\f$ is the identity. For details see Ref. [AM+02] equations (8)-(16).
 *
 * @note These statistical methods stem mostly from robotics, where the state tracking is used to narrow down the position of a vehicle
 * given noisy secondary measurements. See for example Ref. [SC06] for an assessment of snow data assimilation through an
 * (Ensemble) Kalman filter.
 *
 * @subsection kalmanexample Example: Sensor fusion
 *
 * Let's start with a real-world example and work through the settings.
 * @note This example is exercised in the `/doc/examples/statistical_filters.cc` program to follow along.
 *
 * At Innsbruck's Seegrube, there is a spot were two stations from different owners are situated right next to each other.
 * Of course they both measure the air temperature and we want to use the Kalman filter to find the _most likely real temperature_.
 * @note The Kalman filter can only be used for linear models (that can be expressed in matrix form). It is only valid for Gaussian
 * measurement noise.
 *
 * @code
 * TA::FILTER1 = KALMAN
 * @endcode
 *
 * First, we need the heart of the problem, our <b>system model</b>. It describes the evolution of the states through time and is read from
 * the `STATE_DYNAMICS` key. This is the matrix A from above.
 * Matrices can be input in three different ways:
 * 1. A scalar that the matrix is filled with
 * 2. A vector that is put on the diagonal
 * 3. A complete matrix.
 * Matrices and vectors are input by a list of comma-separated elements. To init both a vector that is 4 elements long or a matrix that
 * is 2 by 2 you would write: `1, 0, 0, 1`. You can put brackets around your numbers for readability: `[1, 0][0, 1]`. Be careful to abide
 * to this syntax strictly.
 *
 * We will use a most easy _model_:
 * @code
 * STATE_DYNAMICS = 1
 * @endcode
 * This means that we assume a constant temperature throughout our data window.
 * @note You can use the following substitutions in the system matrix:
 * 1. "dt" will be replaced with the normalized and scaled time delta between the current and the last measurement. At T=0, dt=0. At time
 * of the first observation, dt=1.
 * 2. "meteo(XX)" where "XX" stands for any meteo parameter present in the data set at this time, e. g. "meteo(TA)".
 *
 * Next, we need an <b>initial state</b>. We have two internal states (both of which are energy states), namely the two temperatures.
 * For both filters, you can provide the values as a list (e. g. `268, 271`). You can also leave all of them or some of them empty to make
 * MeteoIO pick the initial state(s) from the 1st available observation(s). You can use the token "average" (case sensitive)
 * to let MeteoIO average over all 1st observations. We want the average between the first measurement of both stations so we write:
 * @code
 * INITIAL_STATE = [average, average]
 * @endcode
 * @note Other valid input formats would be to not use this key (to pick earliest observation for all states) or `[268, ]` (to set
 * the first state to 268 and pick the earliest measurement of the second observation variable for the 2nd state).
 * @note This is usually the line that fixes the number of states and thus the matrix sizes. MeteoIO will complain if they don't match,
 * but input errors could still lead to surprising results.
 * @note MeteoIO will handle cases where the user requests to choose initial values automatically, but nodata values are encountered. It will
 * search for the first valid dataset, and use this for the beginning of the measurement. If any one initial state encounters
 * nodata values only, an error is thrown. If the filter itself encounters nodata values a warning is displayed.
 * This makes it run but could have undesired side effects. You should resample beforehand via `[1DINTERPOLATIONS]`.
 * All filters are called after resampling.
 *
 * We have our states initialized and our state transition matrix ready. Now we need some <b>statistical quantities</b>.
 * We provide the trust in the initial state via the `INITIAL_TRUST` key. This inputs the matrix `P` from above. <i>Set low values if you
 * don't trust \f$\hat{x}_0\f$, and very high ones if you do</i>. Let's not worry about it too much and pick a trusty 1 at the diagonal:
 * @code
 * INITIAL_TRUST = 1
 * @endcode
 * @note Other valid input formats that do exactly the same would be `[1, 0] [0, 1]` (full matrix), `[1, 1]` (diagonal vector), or
 * `1, 1` (you do not need brackets).
 *
 * Next, we need the <b>covariance of the process noise</b> (noise in the model function by external forces acting on the state variables).
 * A poor process model might be mitigated by injecting enough uncertainty here. We don't have knowledge of any specific process
 * disturbances but we know that our model is quite crude, so we _inject a bit of uncertainty in the model_:
 * @code
 * PROCESS_COVARIANCE = 0.05
 * @endcode
 * @note The above line will construct this 2x2-matrix: `[0.05, 0][0, 0.05]`.
 *
 * We now have everything regarding the system model. On to the _observations_.
 *
 * The matrix \f$H\f$ from above is the <b>observation model</b> relating our measurements to the internal states. In our case, we measure
 * the temperatures which are also our internal states, so there is nothing to do. We can however decide how much each station contributes
 * to the estimate. Let's trust both stations the same:
 * @code
 * OBSERVATION_RELATION = [1, 1]
 * @endcode
 * Next, the matrix \f$R\f$ which is the observation covariance. Let's assume a standard deviation sigma of 0.8, i. e. the uncertainty of the
 * sensor is somewhere around 0.8 °C and we put the variance (sigma squared) on the diagonal:
 * @code
 * OBSERVATION_COVARIANCE = 0.64
 * @endcode
 *
 * @note If \f$Q\f$ and \f$R\f$ are time invariant then \f$P\f$ and \f$K\f$ will stabilize quickly, cf. Fig. 4 (they could even be computed offline).
 *
 * We don't want to include any _control signal_ \f$u\f$, so we can neglect to provide `CONTROL_SIGNAL` and `CONTROL_RELATION`.
 * Only thing missing is a list of parameters that stand for the states. The filter runs on `TA` which is the first state, and we want
 * the second one to be `TA_HYD`. By default, MeteoIO will only filter `TA` in the output file. This is to comply with all other filters and avoid
 * unexpected results. We can set however that all observations that are filtered internally are also filtered in the meteo set:
 * @code
 * ADD_OBSERVABLES       = TA_HYD
 * FILTER_ALL_PARAMETERS = TRUE
 * @endcode
 * Note that we require the number of observables to be equal to the number of states, simply to avoid having to set how
 * many, and which, observables to output. You must supply them, but you can easily discard them (e. g. by leaving
 * `FILTER_ALL_PARAMETERS = FALSE` or by designing your matrices in a way that they zero out).
 *
 * @note If you want to output all filtered parameters but also keep the originals you can simply use MeteoIO's `COPY` command.
 *
 * @note Even though MeteoIO can handle nodata cases, it is advised that you _resample your meteo data first_ via `[1DINTERPOLATIONS]`.
 * As the number of nodata values increases, the results will become more and more skewed. A warning will be displayed.
 *
 * The complete ini section now reads:
 * @code
 * [FILTERS]
 * TA::FILTER1                          =   KALMAN
 * TA::ARG1::STATE_DYNAMICS             =   1
 * TA::ARG1::INITIAL_STATE              =   [average, average]
 * TA::ARG1::INITIAL_TRUST              =   1
 * TA::ARG1::PROCESS_COVARIANCE         =   0.05
 * TA::ARG1::ADD_OBSERVABLES            =   TA_HYD
 * TA::ARG1::FILTER_ALL_PARAMETERS      =   TRUE
 * TA::ARG1::OBSERVATION_RELATION       =   1
 * TA::ARG1::OBSERVATION_COVARIANCE     =   0.6
 * @endcode
 *
 * \image html kalmanfilter_trust_observation.png "Fig. 1: Kalman filtering result with a process covariance of 0.05 and an observation covariance of 0.6, i. e. there is a little uncertainty injected in the system model and the observation is expected to be noisy with approximately sqrt(0.6) degrees uncertainty."
 *
 * \image html kalmanfilter_trust_model.png "Fig. 2: Kalman filtering result with a process covariance of 0 and an observation covariance of 0.6, i. e. the model is trusted completely."
 *
 * \image html kalmanfilter_trust_medium.png "Fig. 3: Kalman filtering result with a process covariance of 0.2 and an observation covariance of 0.3, i. e. a middle ground where spikes are followed, but not to an extreme extent."
 *
 * @subsection kalmanfeatures Other features
 *
 * @subsubsection kalmancontrol Control signal input
 *
 * With `CONTROL_SIGNAL` you can input an _external control signal acting on the states_ in three different ways:
 * 1. A scalar (e. g. `0.2`) that an appropriate vector is filled with
 * 2. A vector (e. g. `[0.1, 0.2]`) which will be applied for all time steps
 * 3. A list of meteo parameters holding the value for each timestep (e. g. `UU1 UU2`), i. e. a time-variant vector
 * The matrix \f$H\f$ is applied to this vector and can be input via `CONTROL_RELATION` as usual.
 *
 * @subsubsection kalmanerror Error estimation
 *
 * Say you were in a car entering a tunnel. Then you could still measure your speed, but not your position. The standard deviation
 * of the position would grow (less reliable), and the standard deviation of the speed would decrease (more reliable with new
 * measurements). That would show in the matrix \f$P\f$.
 *
 * `OUT_ESTIMATED_ERROR` lets you pick one or more field names that the filter will output the _evolution of the diagonal elements of
 * the error estimation_ \f$P\f$ to.
 * @code
 * OUT_ESTIMATED_ERROR = ERR1 ERR2
 * @endcode
 * If you have two states, you need two parameter names (that will be created).
 *
 * \image html kalmanfilter_errorestimation.png "Fig. 4: Diagonal components of the error estimation matrix from the sensor fusion example."
 *
 * @subsubsection kalmanvariances Time-variant covariance matrices
 *
 * With `PROCESS_COVARIANCE_PARAMS` and `OBSERVATION_COVARIANCE_PARAMS` you can supply a list of parameters that \f$Q\f$ and \f$R\f$
 * are read from. This has priority over a fixed matrix that is found in `PROCESS_COVARIANCE` and `OBSERVATION_COVARIANCE`.
 * They can only be _input as vectors_ that will be put on a matrix diagonal with the rest being 0. For this you need one meteo parameter
 * per vector element (state) and you need to name the parameters with `PROCESS_COVARIANCE_PARAMS` and `OBSERVATION_COVARIANCE_PARAMS`.
 * @code
 * PROCESS_COVARIANCE_PARAMS = COV_P1 COV_P2 //if either is missing, this key will be ignored
 * OBSERVATION_COVARIANCE_PARAMS = COV_O1 COV_O2
 * @endcode
 *
 * @subsubsection kalmanverbose Suppression of warnings
 *
 * The `VERBOSE` keyword can be used to enable/disable some warnings the filter shows in the console. This will normally concern missing
 * input that was defaulted or something similar. _Usually you will want to get rid of all of them_, but if the filter must be quiet no
 * matter what then you can set `VERBOSE = FALSE`.
 *
 * @subsubsection kalmannotes Further notes
 *
 * - While matrices are input as comma-separated values, other lists are read from a space delimited line (e. g. `PROCESS_COVARIANCE_PARAMS`).
 * - You may find MeteoIO's `[DATASOURCEXXX]` command to read multiple input files and stream them to the same container
 * useful to combine meteo data and model input.
 * - <b>Important:</b> For technical reasons the Kalman and particle filters can run on a larger time frame than it is
 * requested in the getMeteoData call. This allowes windowed filters to have enough data even at the beginning and end to do
 * their work. This is controlled by the `BUFF_BEFORE` and `BUFFER_SIZE` keys. If you set both to 0 then the filters receive
 * exactly the time frame covered by the input dates. However, if you are using windowed filters also, then these will lack
 * data. This is important in understanding how the initial states are picked: if MeteoIO inits the states automatically with
 * the first available data point then this is the first point MeteoIO sees. Thus, with a large enough buffer the model will
 * start at a much earlier date and initial states are also picked a while before the input start date. Additionally, if there
 * are nodata elements in there then even if they would be ignored for the they would make the filters display warnings.
 * You have to be careful to either start your analytical model at the date `data window - buffer`, or ideally provide exactly
 * the same amount of data on the file system as you request MeteoIO to filter.
 *
 * @subsection kalmankeylist List of ini keys
 *  <table>
 *  <tr><th>Keyword</th><th>Meaning</th><th>Optional</th><th>Default Value</th></tr>
 *  <tr><td>STATE_DYNAMICS</td><td>State transition matrix (\f$A\f$).</td><td>no</td><td>-</td></tr>
 *  <tr><td>INITIAL_STATE</td><td>State of the system at T=0 (\f$\hat{x}_0\f$).</td><td>yes</td><td>empty, meaning 1st measurement will be picked</td></tr>
 *  <tr><td>INITIAL_TRUST</td><td>Trust in the initial state (\f$P_0\f$).</td><td>yes</td><td>1, meaning an appropriate identity matrix</td></tr>
 *  <tr><td>PROCESS_COVARIANCE</td><td>Process noise covariance matrix (trust in the model, \f$Q\f$).</td><td>yes</td><td>0, meaning a matrix with all zeros</td></tr>
 *  <tr><td>ADD_OBSERVABLES</td><td>List of observables in addition to the one the filter runs on.</td><td>yes</td><td>empty, meaning we have 1 state and 1 observable</td></tr>
 *  <tr><td>FILTER_ALL_PARAMETERS</td><td>Filter only the parameter the filter runs on, or all states?</td><td>yes</td><td>FALSE, but TRUE is useful</td></tr>
 *  <tr><td>OBSERVATION_RELATION</td><td>Linear model relating the observations to the states (\f$H\f$).</td><td>yes</td><td>1, meaning identity matrix</td></tr>
 *  <tr><td>OBSERVATION_COVARIANCE</td><td>Observation noise covariance matrix (trust in the measurements, \f$R\f$).</td><td>yes</td><td>0</td></tr>
 *  <tr><td>CONTROL_SIGNAL</td><td>External control signal acting on the states (\f$\hat{u}\f$).</td><td>yes</td><td>0</td></tr>
 *  <tr><td>OUT_ESTIMATED_ERROR</td><td>Parameter name(s) to output error evolution to.</td><td>yes</td><td>empty</td></tr>
 *  <tr><td>PROCESS_COVARIANCE_PARAMS</td><td>List of parameter names to read \f$Q\f$ from.</td><td>yes</td><td>empty</td></tr>
 *  <tr><td>OBSERVATION_COVARIANCE_PARAMS</td><td>List of parameter names to read \f$R\f$ from.</td><td>yes</td><td>empty</td></tr>
 *  <tr><td>VERBOSE</td><td>Output warnings to the console.</td><td>yes</td><td>TRUE, warnings should be mitigated</td></tr>
 * </table>
 *
 * <b>Read on</b> at FilterParticle.
 *
 */

class FilterKalman : public ProcessingBlock {
	public:
		FilterKalman(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		        std::vector<MeteoData>& ovec);
	private:
		Eigen::VectorXd buildInitialStates(const std::vector<std::string>& xx_str, const std::vector<size_t>& meas_idx,
		        const std::vector<MeteoData>& ivec, const size_t& nr_observations);
		size_t buildObservationsMatrix(const unsigned int& param, const std::vector<MeteoData>& ivec, const size_t& nx,
		        Eigen::MatrixXd& zz, std::vector<size_t>& meas_idx) const;
		Eigen::MatrixXd buildControlSignal(const size_t& nx, const size_t& TT, const std::vector<MeteoData>& ivec) const;
		Eigen::MatrixXd parseMatrix(const std::string& line, const size_t& rows, const size_t& cols,
		        const std::string& block) const;
		std::vector<std::string> parseSystemMatrix(const std::string& line, const size_t& rows) const;
		Eigen::MatrixXd buildSystemMatrix(const std::vector<std::string>& AA_str, const size_t& sz, const double& dt,
		        const std::vector<MeteoData>& ivec, const size_t& kk) const;
		double substitute(const std::string& expr, const double& dt, const std::vector<MeteoData>& ivec, const size_t& kk) const;
		Eigen::MatrixXd extractInputMatrix(const std::vector<std::string>& vecParams, const std::vector<MeteoData>& mvec,
		        const size_t& kk) const;
		void assertInputCovariances(const std::vector<MeteoData>& ivec, const size_t& nx, bool& has_RR_params,
		        bool& has_QQ_params) const;
		Eigen::MatrixXd bloatMatrix(const std::string& line, const size_t& rows, const size_t& cols, const std::string& block) const;
		Eigen::VectorXd buildTimeVector(const std::vector<MeteoData>& ivec) const;
		bool checkNodata(const Eigen::VectorXd& ivec) const;
		bool findFirstDatapoint(const std::vector<MeteoData>& ivec, const size_t& param, double& retval);
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void cleanBrackets(std::string& iline);

		std::string mat_in_xx, mat_in_AA, mat_in_HH, mat_in_PP, mat_in_QQ, mat_in_RR, mat_in_BB, mat_in_uu;
		std::vector<std::string> meas_params; //parameter names of observations
		std::vector<std::string> error_params; //output PP to these parameters
		std::vector<std::string> QQ_params, RR_params; //input QQ and RR from these parameters

		bool filter_all_params; //filter only the parameter the filter runs on or all the ones used?

		bool be_verbose; //output warnings/info?
		std::string unrecognized_keys; //to warn about unknown ini keys
};

} //namespace

#endif //FILTERKALMAN_H
