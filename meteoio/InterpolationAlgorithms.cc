/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/InterpolationAlgorithms.h>
#include <meteoio/meteolaws/Atmosphere.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/IOManager.h>
#include <meteoio/MathOptim.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>

using namespace std;

namespace mio {

InterpolationAlgorithm* AlgorithmFactory::getAlgorithm(const std::string& i_algoname,
                                                       Meteo2DInterpolator& i_mi,
                                                       const std::vector<std::string>& i_vecArgs, IOManager& iom)
{
	const std::string algoname( IOUtils::strToUpper(i_algoname) );

	if (algoname == "NONE"){// return a nodata grid
		return new NoneAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "CST"){// constant fill
		return new ConstAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "STD_PRESS"){// standard air pressure interpolation
		return new StandardPressureAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "CST_LAPSE"){// constant fill with an elevation lapse rate
		return new ConstLapseRateAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "IDW"){// Inverse Distance Weighting fill
		return new IDWAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "IDW_LAPSE"){// Inverse Distance Weighting with an elevation lapse rate fill
		return new IDWLapseAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "LIDW_LAPSE"){// Inverse Distance Weighting with an elevation lapse rate fill, restricted to a local scale
		return new LocalIDWLapseAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "RH"){// relative humidity interpolation
		return new RHAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "ILWR"){// long wave radiation interpolation
		return new ILWRAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "WIND_CURV"){// wind velocity interpolation (using a heuristic terrain effect)
		return new SimpleWindInterpolationAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "RYAN"){// RYAN wind direction
		return new RyanAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "WINSTRAL"){// Winstral wind exposure factor
		return new WinstralAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "ODKRIG"){// ordinary kriging
		return new OrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "ODKRIG_LAPSE"){// ordinary kriging with lapse rate
		return new LapseOrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "USER"){// read user provided grid
		return new USERInterpolation(i_mi, i_vecArgs, i_algoname, iom);
	} else if (algoname == "HNW_SNOW"){// precipitation interpolation according to (Magnusson, 2010)
		return new SnowHNWInterpolation(i_mi, i_vecArgs, i_algoname, iom);
	} else {
		throw IOException("The interpolation algorithm '"+algoname+"' is not implemented" , AT);
	}
}

size_t InterpolationAlgorithm::getData(const Date& i_date, const MeteoData::Parameters& i_param,
                                       std::vector<double>& o_vecData)
{
	iomanager.getMeteoData(i_date, vecMeteo);
	o_vecData.clear();
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		const double& val = vecMeteo[ii](i_param);
		if (val != IOUtils::nodata) {
			o_vecData.push_back(val);
		}
	}

	return o_vecData.size();
}

size_t InterpolationAlgorithm::getData(const Date& i_date, const MeteoData::Parameters& i_param,
                                       std::vector<double>& o_vecData, std::vector<StationData>& o_vecMeta)
{
	iomanager.getMeteoData(i_date, vecMeteo);
	o_vecData.clear();
	o_vecMeta.clear();
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		const double& val = vecMeteo[ii](i_param);
		if (val != IOUtils::nodata){
			o_vecData.push_back(val);
			o_vecMeta.push_back(vecMeteo[ii].meta);
		}
	}

	return o_vecData.size();
}

size_t InterpolationAlgorithm::getStationAltitudes(const std::vector<StationData>& i_vecMeta,
                                                   std::vector<double>& o_vecData)
{
	o_vecData.clear();
	for (size_t ii=0; ii<i_vecMeta.size(); ii++){
		const double& alt = i_vecMeta[ii].position.getAltitude();
		if (alt != IOUtils::nodata) {
			o_vecData.push_back(alt);
		}
	}

	return o_vecData.size();
}

/**
 * @brief Return an information string about the interpolation process
 * @return string containing some information (algorithm used, number of stations)
*/
std::string InterpolationAlgorithm::getInfo() const
{
	std::ostringstream os;
	os << algo << ", " << nrOfMeasurments;
	if(nrOfMeasurments==1)
		os << " station";
	else
		os << " stations";
	const std::string tmp( info.str() );
	if(!tmp.empty()) {
		os << ", " << tmp;
	}
	return os.str();
}

/**
 * @brief Read the interpolation arguments and compute the trend accordingly
 *
 * @param vecAltitudes altitudes sorted similarly as the data in vecDat
 * @param vecDat data for the interpolated parameter
 * @param trend object containing the fitted trend to be used for detrending/retrending
*/
void InterpolationAlgorithm::getTrend(const std::vector<double>& vecAltitudes, const std::vector<double>& vecDat, Fit1D &trend) const
{
	bool status;
	if (vecArgs.empty()) {
		trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
		status = trend.fit();
	} else if (vecArgs.size() == 1) {
		double lapse_rate;
		IOUtils::convertString(lapse_rate, vecArgs.front());
		trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
		trend.setLapseRate(lapse_rate);
		status = trend.fit();
	} else if (vecArgs.size() == 2) {
		std::string extraArg;
		IOUtils::convertString(extraArg, vecArgs[1]);
		if(extraArg=="soft") { //soft
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			status = trend.fit();
			if(!status) {
				double lapse_rate;
				IOUtils::convertString(lapse_rate, vecArgs[0]);
				trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
				trend.setLapseRate(lapse_rate);
				status = trend.fit();
			}
		} else if(extraArg=="frac") {
			double lapse_rate;
			IOUtils::convertString(lapse_rate, vecArgs[0]);
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			const double avgData = Interpol1D::arithmeticMean(vecDat);
			trend.setLapseRate(lapse_rate*avgData);
			status = trend.fit();
		} else {
			throw InvalidArgumentException("Unknown argument \""+extraArg+"\" supplied for the "+algo+" algorithm", AT);
		}
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	}

	if(!status)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param) + ": " + trend.getInfo(), AT);
}


void InterpolationAlgorithm::detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat, const double& min_alt, const double& max_alt)
{
	if(vecDat.size() != vecAltitudes.size()) {
		std::ostringstream ss;
		ss << "Number of station data (" << vecDat.size() << ") and number of elevations (" << vecAltitudes.size() << ") don't match!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for(size_t ii=0; ii<vecAltitudes.size(); ii++) {
		const double altitude = std::min( std::max(vecAltitudes[ii], min_alt), max_alt );
		double &val = vecDat[ii];
		if(val!=IOUtils::nodata)
			val -= trend( altitude );
	}
}

void InterpolationAlgorithm::retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid, const double& min_alt, const double& max_alt)
{
	const size_t nxy = grid.getNx()*grid.getNy();
	const size_t dem_nxy = dem.grid2D.getNx()*dem.grid2D.getNy();
	if(nxy != dem_nxy) {
		std::ostringstream ss;
		ss << "Dem size (" << dem.grid2D.getNx() << "," << dem.grid2D.getNy() << ") and";
		ss << "grid size (" << grid.getNx() << "," << grid.getNy() << ") don't match!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for(size_t ii=0; ii<nxy; ii++) {
		const double altitude = std::min( std::max(dem(ii), min_alt), max_alt );
		double &val = grid(ii);
		if(val!=IOUtils::nodata)
			val += trend.f( altitude );
	}
}

/**********************************************************************************/
/*                    Implementation of the various algorithms                    */
/**********************************************************************************/
double NoneAlgorithm::getQualityRating(const Date& /*i_date*/, const MeteoData::Parameters& /*in_param*/)
{
	return 1e-6;
}

void NoneAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid) {
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner, IOUtils::nodata);
}


double StandardPressureAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (param != MeteoData::P)
		return 0.0;

	if (nrOfMeasurments == 0)
		return 1.0;

	return 0.1;
}

void StandardPressureAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid) {
	Interpol2D::stdPressure(dem, grid);
}


double ConstAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0) {
		const size_t nr_args = vecArgs.size();

		if (nr_args>1)
			throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
		if (nr_args==1) {
			IOUtils::convertString(user_cst, vecArgs[0]);
			user_provided = true;
			return 0.01;
		} else
			return 0.0;
	} else if (nrOfMeasurments == 1) {
		return 0.8;
	} else if (nrOfMeasurments > 1) {
		return 0.2;
	}

	return 0.0;
}

void ConstAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid) {
	if (!user_provided)
		Interpol2D::constant(Interpol1D::arithmeticMean(vecData), dem, grid);
	else
		Interpol2D::constant(user_cst, dem, grid);
}


double ConstLapseRateAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0){
		return 0.0;
	} else if (nrOfMeasurments == 1){
		if (!vecArgs.empty())
			return 0.9; //the laspe rate is provided
		else
			return 0.0; //no lapse rate is provided and it can not be computed
	} else if (nrOfMeasurments == 2){
		//in any case, we can do at least as good as IDW_LAPSE
		return 0.71;
	} else if (nrOfMeasurments>2){
		return 0.2;
	}

	return 0.2;
}

void ConstLapseRateAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Fit1D trend;
	getTrend(vecAltitudes, vecData, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);
	Interpol2D::constant(Interpol1D::arithmeticMean(vecData), dem, grid);
	retrend(dem, trend, grid);
}


double IDWAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0){
		return 0.0;
	} else if (nrOfMeasurments == 1){
		return 0.3;
	} else if (nrOfMeasurments > 1){
		return 0.5;
	}

	return 0.2;
}

void IDWAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	Interpol2D::IDW(vecData, vecMeta, dem, grid);
}


double IDWLapseAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.7;
}

void IDWLapseAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Fit1D trend;
	getTrend(vecAltitudes, vecData, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);
	Interpol2D::IDW(vecData, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);
}


LocalIDWLapseAlgorithm::LocalIDWLapseAlgorithm(Meteo2DInterpolator& i_mi, const std::vector<std::string>& i_vecArgs,
                                               const std::string& i_algo, IOManager& iom)
                      : InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom), nrOfNeighbors(0)
{
	if (vecArgs.size() == 1) { //compute lapse rate on a reduced data set
		IOUtils::convertString(nrOfNeighbors, vecArgs[0]);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Please provide the number of nearest neighbors to use for the "+algo+" algorithm", AT);
	}
}

double LocalIDWLapseAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.7;
}

void LocalIDWLapseAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	if (nrOfMeasurments == 0)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	Interpol2D::LocalLapseIDW(vecData, vecMeta, dem, nrOfNeighbors, grid);
	info << "using nearest " << nrOfNeighbors << " neighbors";
}


double RHAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for RH
	if (in_param != MeteoData::RH)
		return 0.0;

	date = i_date;
	param = in_param;
	vecData.clear(); vecMeta.clear();
	vecDataTA.clear(); vecDataRH.clear();

	nrOfMeasurments = 0;
	iomanager.getMeteoData(date, vecMeteo);
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii](MeteoData::RH) != IOUtils::nodata) && (vecMeteo[ii](MeteoData::TA) != IOUtils::nodata)){
			vecDataTA.push_back(vecMeteo[ii](MeteoData::TA));
			vecDataRH.push_back(vecMeteo[ii](MeteoData::RH));
			vecMeta.push_back(vecMeteo[ii].meta);
			nrOfMeasurments++;
		}
	}

	if (nrOfMeasurments==0)
		return 0.0;
	if( (nrOfMeasurments<vecDataRH.size()/2) || ( nrOfMeasurments<2 ) )
		return 0.6;

	return 0.9;
}

void RHAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator

	//RH->Td, interpolations, Td->RH
	//Compute dew point temperatures at stations
	std::vector<double> vecTd(vecDataRH.size());
	for (size_t ii=0; ii<vecDataRH.size(); ii++){
		vecTd[ii] = Atmosphere::RhtoDewPoint(vecDataRH[ii], vecDataTA[ii], 1);
	}

	Fit1D trend;
	getTrend(vecAltitudes, vecTd, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecTd);
	Interpol2D::IDW(vecTd, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);

	//Recompute Rh from the interpolated td
	for (size_t jj=0; jj<grid.nrows; jj++) {
		for (size_t ii=0; ii<grid.ncols; ii++) {
			double &value = grid(ii,jj);
			if(value!=IOUtils::nodata)
				value = Atmosphere::DewPointtoRh(value, ta(ii,jj), 1);
		}
	}
}


double ILWRAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for ILWR
	if (in_param != MeteoData::ILWR)
		return 0.0;

	date = i_date;
	param = in_param;
	vecData.clear(); vecMeta.clear();
	vecDataEA.clear();

	nrOfMeasurments = 0;
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii](MeteoData::ILWR) != IOUtils::nodata) && (vecMeteo[ii](MeteoData::TA) != IOUtils::nodata)){
			vecDataEA.push_back( Atmosphere::blkBody_Emissivity( vecMeteo[ii](MeteoData::ILWR), vecMeteo[ii](MeteoData::TA)) );
			vecMeta.push_back(vecMeteo[ii].meta);
			nrOfMeasurments++;
		}
	}

	if (nrOfMeasurments==0)
		return 0.0;

	return 0.9;
}

void ILWRAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator

	Fit1D trend;
	getTrend(vecAltitudes, vecDataEA, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecDataEA);
	Interpol2D::IDW(vecDataEA, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);

	//Recompute Rh from the interpolated td
	for (size_t jj=0; jj<grid.nrows; jj++) {
		for (size_t ii=0; ii<grid.ncols; ii++) {
			double &value = grid(ii,jj);
			value = Atmosphere::blkBody_Radiation(value, ta(ii,jj));
		}
	}
}


double SimpleWindInterpolationAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for VW
	if (in_param != MeteoData::VW)
		return 0.0;

	date = i_date;
	param = in_param;
	vecData.clear(); vecMeta.clear();
	vecDataVW.clear(); vecDataDW.clear();

	nrOfMeasurments = 0;
	iomanager.getMeteoData(date, vecMeteo);
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii](MeteoData::VW) != IOUtils::nodata) && (vecMeteo[ii](MeteoData::DW) != IOUtils::nodata)){
			vecDataVW.push_back(vecMeteo[ii](MeteoData::VW));
			vecDataDW.push_back(vecMeteo[ii](MeteoData::DW));
			vecMeta.push_back(vecMeteo[ii].meta);
			nrOfMeasurments++;
		}
	}

	if (nrOfMeasurments==0)
		return 0.0;
	if( (nrOfMeasurments<vecDataVW.size()/2) || ( nrOfMeasurments<2 ) )
		return 0.6;

	return 0.9;
}

void SimpleWindInterpolationAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (Interpol2D::allZeroes(vecDataVW)) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Grid2DObject dw;
	mi.interpolate(date, dem, MeteoData::DW, dw); //get DW interpolation from call back to Meteo2DInterpolator

	Fit1D trend;
	getTrend(vecAltitudes, vecDataVW, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecDataVW);
	Interpol2D::IDW(vecDataVW, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);
	Interpol2D::SimpleDEMWindInterpolate(dem, grid, dw);
}


double RyanAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for DW (we could add VW later)
	if (in_param!=MeteoData::DW)
		return 0.0;

	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments==0)
		return 0.0;
	if (nrOfMeasurments<2)
		return 0.6;

	return 0.9;
}

void RyanAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (Interpol2D::allZeroes(vecData)) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	Interpol2D::IDW(vecData, vecMeta, dem, grid);
	Interpol2D::RyanWindDir(dem, grid);
}


const double WinstralAlgorithm::dmax = 300.;

double WinstralAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for HNW (we could add HS later)
	if (in_param!=MeteoData::HNW)
		return 0.0;

	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments==0)
		return 0.0;

	return 0.99;
}

void WinstralAlgorithm::getParameters(std::string &base_algo, double &synoptic_bearing) const
{
	//option 0: compute DW from multiple stations
	//option 1: fixed DW
	//option 2: specify the station_id to take DW from
	base_algo = "IDW_LAPSE";

	const size_t nr_args = vecArgs.size();
	if (nr_args==0) {
		synoptic_bearing = getSynopticBearing(vecMeteo);
		if (synoptic_bearing==IOUtils::nodata)
			throw IOException("Synoptic wind could not be computed by the "+algo+" algorithm", AT);
		return;
	} else if (nr_args==1) {
		if (IOUtils::isNumeric(vecArgs[0]))
			IOUtils::convertString(synoptic_bearing, vecArgs[0]);
		else
			throw InvalidArgumentException("Please provide both the base interpolation method and the station_id to use for wind direction for the "+algo+" algorithm", AT);
		return;
	} else if (nr_args==2) {
		if (IOUtils::isNumeric(vecArgs[0])) {
			IOUtils::convertString(synoptic_bearing, vecArgs[0]);
			base_algo = IOUtils::strToUpper( vecArgs[1] );
		} else if (IOUtils::isNumeric(vecArgs[1]))  {
			IOUtils::convertString(synoptic_bearing, vecArgs[1]);
			base_algo = IOUtils::strToUpper( vecArgs[0] );
		} else {
			base_algo = IOUtils::strToUpper( vecArgs[0] );
			const string ref_station( vecArgs[1] );
			synoptic_bearing = getSynopticBearing(vecMeteo, ref_station, algo);
		}
		return;
	} else
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
}

void WinstralAlgorithm::initGrid(const std::string& base_algo, const DEMObject& dem, Grid2DObject& grid)
{
	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	vector<string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, iomanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();
}

double WinstralAlgorithm::getSynopticBearing(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station, const std::string& algo)
{
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (vecMeteo[ii].meta.stationID==ref_station)
			return vecMeteo[ii](MeteoData::DW);
	}

	throw InvalidArgumentException("Station \""+ref_station+"\" should be used as reference station for the "+algo+" algorithm but has not been found", AT);
}

double WinstralAlgorithm::getSynopticBearing(const std::vector<MeteoData>& vecMeteo)
{
	// 1) locate the stations in DEM and check if they are higher than their surroundings within a given radius
	// 2) simply compute a mean or median direction
	// (2) can be used on all the stations selected in (1)

	double ve=0.0, vn=0.0;
	size_t count=0;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		const double VW = vecMeteo[ii](MeteoData::VW);
		const double DW = vecMeteo[ii](MeteoData::DW);
		if(VW!=IOUtils::nodata && DW!=IOUtils::nodata) {
			ve += VW * sin(DW*Cst::to_rad);
			vn += VW * cos(DW*Cst::to_rad);
			count++;
		}
	}

	if (count!=0) {
		ve /= static_cast<double>(count);
		vn /= static_cast<double>(count);

		//const double meanspeed = sqrt(ve*ve + vn*vn);
		const double meandirection = fmod( atan2(ve,vn) * Cst::to_deg + 360., 360.);
		return meandirection;
	}

	return IOUtils::nodata;
}

void WinstralAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (Interpol2D::allZeroes(vecData)) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	string base_algo, ref_station;
	double synoptic_bearing;
	getParameters(base_algo, synoptic_bearing);
	initGrid(base_algo, dem, grid);

	//alter the field with Winstral and the chosen wind direction
	Interpol2D::Winstral(dem, dmax, synoptic_bearing, grid);
}


std::string USERInterpolation::getGridFileName() const
{
//HACK: use read2DGrid(grid, MeteoGrid::Parameters, Date) instead?
	if (vecArgs.size() != 1){
		throw InvalidArgumentException("Please provide the path to the grids for the "+algo+" interpolation algorithm", AT);
	}
	const std::string ext(".asc");
	const std::string& grid_path = vecArgs[0];
	std::string gridname = grid_path + "/";

	if (!vecMeteo.empty()) {
		const Date& timestep = vecMeteo.at(0).date;
		gridname =  gridname + timestep.toString(Date::NUM) + "_" + MeteoData::getParameterName(param) + ext;
	} else {
		gridname = gridname + "Default" + "_" + MeteoData::getParameterName(param) + ext;
	}

	return gridname;
}

double USERInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	filename = getGridFileName();

	if (!IOUtils::validFileName(filename)) {
		cerr << "[E] Invalid grid filename for "+algo+" interpolation algorithm: " << filename << "\n";
		return 0.0;
	}
	if(IOUtils::fileExists(filename)) {
		return 1.0;
	} else {
		return 0.0;
	}
}

void USERInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	iomanager.read2DGrid(grid, filename);
	if(!grid.isSameGeolocalization(dem)) {
		throw InvalidArgumentException("[E] trying to load a grid(" + filename + ") that has not the same georeferencing as the DEM!", AT);
	}
}


double SnowHNWInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.9;
}

void SnowHNWInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//retrieve optional arguments
	std::string base_algo;
	if (vecArgs.empty()){
		base_algo=std::string("IDW_LAPSE");
	} else if (vecArgs.size() == 1){
		IOUtils::convertString(base_algo, vecArgs[0]);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	}

	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	IOUtils::toUpper(base_algo);
	vector<string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, iomanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();

	 //get TA interpolation from call back to Meteo2DInterpolator
	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta);

	//slope/curvature correction for solid precipitation
	Interpol2D::SteepSlopeRedistribution(dem, ta, grid);
	//Interpol2D::CurvatureCorrection(dem, ta, grid);
}


void OrdinaryKrigingAlgorithm::getDataForVariogram(std::vector<double> &distData, std::vector<double> &variData)
{
	distData.clear();
	variData.clear();

	for(size_t j=0; j<nrOfMeasurments; j++) {
		const Coords& st1 = vecMeta[j].position;
		const double x1 = st1.getEasting();
		const double y1 = st1.getNorthing();
		const double val1 = vecData[j];

		for(size_t i=0; i<j; i++) {
			//compute distance between stations
			const Coords& st2 = vecMeta[i].position;
			const double val2 = vecData[i];
			const double DX = x1-st2.getEasting();
			const double DY = y1-st2.getNorthing();
			const double distance = Optim::fastSqrt_Q3( Optim::pow2(DX) + Optim::pow2(DY) );

			distData.push_back(distance);
			variData.push_back( 0.5*Optim::pow2(val1-val2) );
		}
	}
}

bool OrdinaryKrigingAlgorithm::computeVariogram()
{//return variogram fit of covariance between stations i and j
	std::vector<double> distData, variData;
	getDataForVariogram(distData, variData);

	std::vector<string> vario_types( vecArgs );
	if(vario_types.empty()) vario_types.push_back("LINVARIO");

	size_t args_index=0;
	do {
		const string vario_model = IOUtils::strToUpper( vario_types[args_index] );
		const bool status = variogram.setModel(vario_model, distData, variData);
		if(status) {
			info << " - " << vario_model;
			return true;
		}

		args_index++;
	} while (args_index<vario_types.size());

	return false;
}

double OrdinaryKrigingAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if(nrOfMeasurments==0) return 0.;
	if(nrOfMeasurments>=7) return 0.9;
	return 0.1;
}

void OrdinaryKrigingAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//optimization: getRange (from variogram fit -> exclude stations that are at distances > range (-> smaller matrix)
	//or, get max range from io.ini, build variogram from this user defined max range
	if(!computeVariogram()) //only refresh once a month, or once a week, etc
		throw IOException("The variogram for parameter " + MeteoData::getParameterName(param) + " could not be computed!", AT);
	Interpol2D::ODKriging(vecData, vecMeta, dem, variogram, grid);
}


void LapseOrdinaryKrigingAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	//optimization: getRange (from variogram fit -> exclude stations that are at distances > range (-> smaller matrix)
	//or, get max range from io.ini, build variogram from this user defined max range
	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Fit1D trend(Fit1D::NOISY_LINEAR, vecAltitudes, vecData, false);
	if(!trend.fit())
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param) + ": " + trend.getInfo(), AT);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);

	if(!computeVariogram()) //only refresh once a month, or once a week, etc
		throw IOException("The variogram for parameter " + MeteoData::getParameterName(param) + " could not be computed!", AT);
	Interpol2D::ODKriging(vecData, vecMeta, dem, variogram, grid);

	retrend(dem, trend, grid);
}

} //namespace
