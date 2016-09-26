/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <meteoio/meteoFilters/ProcIIRLowPass.h>
#include <meteoio/MathOptim.h>
#include <cmath>

using namespace std;

namespace mio {

const double ProcIIRLowPass::n = 1.; //one pass

//Butterworth
//const double ProcIIRLowPass::g = 1.;
//const double ProcIIRLowPass::p = sqrt(2.);
//const double ProcIIRLowPass::c = 1. / pow( pow(2, 1./n) - 1., 1./4. ); //3dB cutoff correction

//critically damped filter
const double ProcIIRLowPass::g = 1.;
const double ProcIIRLowPass::p = 2.;
const double ProcIIRLowPass::c = 1. / sqrt( pow(2., 1./(2.*n)) - 1. ); //3dB cutoff correction

//Bessel filter
//const double ProcIIRLowPass::g = 3.;
//const double ProcIIRLowPass::p = 3.;
//const double ProcIIRLowPass::c = 1./ sqrt( sqrt(pow(2., 1/n) - 3./4.) - 0.5 ) / sqrt(3.); //3dB cutoff correction

ProcIIRLowPass::ProcIIRLowPass(const std::vector<std::string>& vec_args, const std::string& name)
                  : ProcessingBlock(name), cutoff(0.), bidirectional(true)
{
	parse_args(vec_args);
	properties.points_before = 2;
	properties.stage = ProcessingProperties::first;
}

void ProcIIRLowPass::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ivec.size()<2 || cutoff==0.) return;

	const double days = ivec.back().date.getJulian() - ivec.front().date.getJulian();
	const size_t nr_data_pts = ivec.size();
	const double sampling_rate = static_cast<double>(nr_data_pts-1) / (days*24.*3600.); //in Hz
	double A[3], B[3];
	computeCoefficients(sampling_rate, 1./cutoff, A, B);
	std::vector<double> X(3, IOUtils::nodata), Y(3, IOUtils::nodata);

	if (!bidirectional) {
		//only forward filter
		for (size_t ii=0; ii<ovec.size(); ++ii){
			const double& raw_val = ivec[ii](param);

			//propagate in X and Y
			X[2] = X[1]; X[1] = X[0]; X[0] = raw_val;
			Y[2] = Y[1]; Y[1] = Y[0]; Y[0] = raw_val; //Y[0] will be overwritten but in case of nodata we still propagate a value
			if (X[2]==IOUtils::nodata || X[1]==IOUtils::nodata || X[0]==IOUtils::nodata) continue;
			if (Y[2]==IOUtils::nodata || Y[1]==IOUtils::nodata) continue;

			Y[0] = A[0]*X[0] + A[1]*X[1] + A[2]*X[2] + B[1]*Y[1] + B[2]*Y[2];
			ovec[ii](param) = Y[0];
		}
	} else { //bidirectional filtering
		//because bidirectional filtering introduces some overshoot, we do it twice: backward/forward, then forward/backward and
		//only accept the values that returned close values between the two. Large difference between these two indicate that
		//there is some heavy overshoot going on and therefore we keep the original data.

		std::vector<double> vecTmp(nr_data_pts), vecRes1(nr_data_pts);
		//backward filter
		for (size_t ii=ovec.size(); ii--> 0;){
			const double& raw_val = ivec[ii](param);

			//propagate in X and Y
			X[2] = X[1]; X[1] = X[0]; X[0] = raw_val;
			Y[2] = Y[1]; Y[1] = Y[0]; Y[0] = raw_val; //Y[0] will be overwritten but in case of nodata we still propagate a value
			if (X[2]==IOUtils::nodata || X[1]==IOUtils::nodata || X[0]==IOUtils::nodata) continue;
			if (Y[2]==IOUtils::nodata || Y[1]==IOUtils::nodata) continue;

			Y[0] = A[0]*X[0] + A[1]*X[1] + A[2]*X[2] + B[1]*Y[1] + B[2]*Y[2];
			vecTmp[ii] = Y[0];
		}

		//forward filter
		for (size_t ii=0; ii<ovec.size(); ++ii){
			const double& raw_val = vecTmp[ii];

			//propagate in X and Y
			X[2] = X[1]; X[1] = X[0]; X[0] = raw_val;
			Y[2] = Y[1]; Y[1] = Y[0]; Y[0] = raw_val; //Y[0] will be overwritten but in case of nodata we still propagate a value
			if (X[2]==IOUtils::nodata || X[1]==IOUtils::nodata || X[0]==IOUtils::nodata) continue;
			if (Y[2]==IOUtils::nodata || Y[1]==IOUtils::nodata) continue;

			Y[0] = A[0]*X[0] + A[1]*X[1] + A[2]*X[2] + B[1]*Y[1] + B[2]*Y[2];
			vecRes1[ii] = Y[0];
		}

		std::vector<double> vecRes2(nr_data_pts);
		//forward filter
		for (size_t ii=0; ii<ovec.size(); ++ii){
			const double& raw_val = ivec[ii](param);

			//propagate in X and Y
			X[2] = X[1]; X[1] = X[0]; X[0] = raw_val;
			Y[2] = Y[1]; Y[1] = Y[0]; Y[0] = raw_val; //Y[0] will be overwritten but in case of nodata we still propagate a value
			if (X[2]==IOUtils::nodata || X[1]==IOUtils::nodata || X[0]==IOUtils::nodata) continue;
			if (Y[2]==IOUtils::nodata || Y[1]==IOUtils::nodata) continue;

			Y[0] = A[0]*X[0] + A[1]*X[1] + A[2]*X[2] + B[1]*Y[1] + B[2]*Y[2];
			vecTmp[ii] = Y[0];
		}
		//backward filter
		for (size_t ii=ovec.size(); ii--> 0;){
			const double& raw_val = vecTmp[ii];

			//propagate in X and Y
			X[2] = X[1]; X[1] = X[0]; X[0] = raw_val;
			Y[2] = Y[1]; Y[1] = Y[0]; Y[0] = raw_val; //Y[0] will be overwritten but in case of nodata we still propagate a value
			if (X[2]==IOUtils::nodata || X[1]==IOUtils::nodata || X[0]==IOUtils::nodata) continue;
			if (Y[2]==IOUtils::nodata || Y[1]==IOUtils::nodata) continue;

			Y[0] = A[0]*X[0] + A[1]*X[1] + A[2]*X[2] + B[1]*Y[1] + B[2]*Y[2];
			vecRes2[ii] = Y[0];
		}

		for(size_t ii=0; ii<ovec.size(); ++ii) {
			const double filt1 = vecRes1[ii];
			const double filt2 = vecRes2[ii];
			const double ref = std::max(abs(filt1), abs(filt2));
			const double threshold = (ref>1e-9)? ref*1e-2 : 1e-4;
			if (abs((filt1-filt2)) < threshold)
				ovec[ii](param) = filt1;
		}

	}
}

//this computes the filter coefficients for a low pass filter.
//the filter parameters are computed based on the filter type (this gives the polynomial coefficients and the cutoff correction) and the number of passes.
void ProcIIRLowPass::computeCoefficients(const double& f_s, const double& f_0, double A[3], double B[3]) const
{
	//using the filter polynomials, the number of passes and the cutoff correction, compute the filter coefficients
	const double f_star = c * f_0 / f_s; //corrected cutoff frequency
	const double w_0 = tan(Cst::PI*f_star); //warp cutoff frequency

	const double K1 = p * w_0;
	const double K2 = g * w_0*w_0;

	A[0] = K2 / (1. + K1 + K2);
	A[1] = 2. * A[0];
	A[2] = A[0];

	B[1] = 2*A[0] * (1./K2 - 1.);
	B[2] = 1. - (A[0] + A[1] + A[2] + B[1]);
}

void ProcIIRLowPass::parse_args(std::vector<std::string> vec_args)
{
	const size_t nrArgs = vec_args.size();
	bool period_read = false;

	for (size_t ii=0; ii<nrArgs; ii++) {
		if (IOUtils::isNumeric(vec_args[ii])) {
			if (period_read==true)
				throw InvalidArgumentException("Cutoff period has been provided more than once", AT);
			if (!IOUtils::convertString(cutoff, vec_args[ii]))
				throw InvalidArgumentException("Could not parse cutoff period '"+vec_args[ii]+"'", AT);
			period_read = true;
		} else if (vec_args[ii]=="single_pass") {
			bidirectional = false;
		} else
			throw InvalidArgumentException("Invalid argument \""+vec_args[ii]+"\" for filter \""+getName()+"\"", AT);
	}
	if (!period_read)
		throw InvalidArgumentException("Please provide the cutoff period for filter " + getName(), AT);
}

}
