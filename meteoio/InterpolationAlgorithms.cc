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
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/MathOptim.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>

#include <algorithm>

using namespace std;

namespace mio {

InterpolationAlgorithm* AlgorithmFactory::getAlgorithm(const std::string& i_algoname,
                                                       Meteo2DInterpolator& i_mi,
                                                       const std::vector<std::string>& i_vecArgs, TimeSeriesManager& tsm, GridsManager& gdm)
{
	const std::string algoname( IOUtils::strToUpper(i_algoname) );

	if (algoname == "NONE"){// return a nodata grid
		return new NoneAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "STD_PRESS"){// standard air pressure interpolation
		return new StandardPressureAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "CST"){// constant fill
		return new ConstAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "AVG"){// average fill
		return new AvgAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "AVG_LAPSE"){// average fill with an elevation lapse rate
		return new AvgLapseRateAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "IDW"){// Inverse Distance Weighting fill
		return new IDWAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "IDW_LAPSE"){// Inverse Distance Weighting with an elevation lapse rate fill
		return new IDWLapseAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "LIDW_LAPSE"){// Inverse Distance Weighting with an elevation lapse rate fill, restricted to a local scale
		return new LocalIDWLapseAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "RH"){// relative humidity interpolation
		return new RHAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ILWR"){// long wave radiation interpolation
		return new ILWRAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "LISTON_WIND"){// wind velocity interpolation (using a heuristic terrain effect)
		return new ListonWindAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "RYAN"){// RYAN wind direction
		return new RyanAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "WINSTRAL"){// Winstral wind exposure factor
		return new WinstralAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "WINSTRAL++"){// Winstral/Liston wind exposure factor
		return new WinstralListonAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ODKRIG"){// ordinary kriging
		return new OrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ODKRIG_LAPSE"){// ordinary kriging with lapse rate
		return new LapseOrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "USER"){// read user provided grid
		return new USERInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "ALS_SCALING"){// scale from ALS grid
		return new ALS_Interpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "PPHASE"){// precipitation phase parametrization
		return new PPHASEInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "PSUM_SNOW"){// precipitation interpolation according to (Magnusson, 2010)
		return new SnowPSUMInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	} else if (algoname == "SWRAD"){// terrain shadding interpolation
		return new SWRadInterpolation(i_mi, i_vecArgs, i_algoname, tsm, gdm);
	}else {
		throw IOException("The interpolation algorithm '"+algoname+"' is not implemented" , AT);
	}
}

size_t InterpolationAlgorithm::getData(const Date& i_date, const MeteoData::Parameters& i_param,
                                       std::vector<double>& o_vecData)
{
	tsmanager.getMeteoData(i_date, vecMeteo);
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
	tsmanager.getMeteoData(i_date, vecMeteo);
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
	if (nrOfMeasurments==1)
		os << " station";
	else
		os << " stations";
	const std::string tmp( info.str() );
	if (!tmp.empty()) {
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
		if (extraArg=="soft") { //soft
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			status = trend.fit();
			if (!status) {
				double lapse_rate;
				IOUtils::convertString(lapse_rate, vecArgs[0]);
				trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
				trend.setLapseRate(lapse_rate);
				status = trend.fit();
			}
		} else if (extraArg=="frac") {
			double lapse_rate;
			IOUtils::convertString(lapse_rate, vecArgs[0]);
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			const double avgData = Interpol1D::arithmeticMean(vecDat);
			trend.setLapseRate(lapse_rate*avgData);
			status = trend.fit();
			if (lapse_rate*avgData==0.) trend.setInfo(trend.getInfo() + " (null average input for frac lapse rate)");
		} else {
			throw InvalidArgumentException("Unknown argument \""+extraArg+"\" supplied for the "+algo+" algorithm", AT);
		}
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	}

	if (!status)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param) + ": " + trend.getInfo(), AT);
}


void InterpolationAlgorithm::detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat, const double& min_alt, const double& max_alt)
{
	if (vecDat.size() != vecAltitudes.size()) {
		std::ostringstream ss;
		ss << "Number of station data (" << vecDat.size() << ") and number of elevations (" << vecAltitudes.size() << ") don't match!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for (size_t ii=0; ii<vecAltitudes.size(); ii++) {
		const double altitude = std::min( std::max(vecAltitudes[ii], min_alt), max_alt );
		double &val = vecDat[ii];
		if (val!=IOUtils::nodata)
			val -= trend( altitude );
	}
}

void InterpolationAlgorithm::retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid, const double& min_alt, const double& max_alt)
{
	const size_t nxy = grid.getNx()*grid.getNy();
	const size_t dem_nxy = dem.grid2D.getNx()*dem.grid2D.getNy();
	if (nxy != dem_nxy) {
		std::ostringstream ss;
		ss << "Dem size (" << dem.grid2D.getNx() << "," << dem.grid2D.getNy() << ") and";
		ss << "grid size (" << grid.getNx() << "," << grid.getNy() << ") don't match!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	for (size_t ii=0; ii<nxy; ii++) {
		const double altitude = std::min( std::max(dem(ii), min_alt), max_alt );
		double &val = grid(ii);
		if (val!=IOUtils::nodata)
			val += trend.f( altitude );
	}
}

//this interpolates VW, DW by converting to u,v and then doing IDW_LAPSE before reconverting to VW, DW
void InterpolationAlgorithm::simpleWindInterpolate(const DEMObject& dem, const std::vector<double>& vecDataVW, const std::vector<double>& vecDataDW, Grid2DObject &VW, Grid2DObject &DW)
{
	if (vecDataVW.size() != vecDataDW.size())
		throw InvalidArgumentException("VW and DW vectors should have the same size!", AT);

	//compute U,v
	std::vector<double> Ve, Vn;
	for (size_t ii=0; ii<vecDataVW.size(); ii++) {
		Ve.push_back( vecDataVW[ii]*sin(vecDataDW[ii]*Cst::to_rad) );
		Vn.push_back( vecDataVW[ii]*cos(vecDataDW[ii]*Cst::to_rad) );
	}

	//spatially interpolate U,V
	vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating wind", AT);

	if (vecDataVW.size()>=4) { //at least for points to perform detrending
		Fit1D trend;

		getTrend(vecAltitudes, Ve, trend);
		info << trend.getInfo();
		detrend(trend, vecAltitudes, Ve);
		Interpol2D::IDW(Ve, vecMeta, dem, VW);
		retrend(dem, trend, VW);

		getTrend(vecAltitudes, Vn, trend);
		info << trend.getInfo();
		detrend(trend, vecAltitudes, Vn);
		Interpol2D::IDW(Vn, vecMeta, dem, DW);
		retrend(dem, trend, DW);
	} else {
		Interpol2D::IDW(Ve, vecMeta, dem, VW);
		Interpol2D::IDW(Vn, vecMeta, dem, DW);
	}

	//recompute VW, DW in each cell
	const size_t nrCells = VW.getNx()*VW.getNy();
	for (size_t ii=0; ii<nrCells; ii++) {
		const double ve = VW(ii);
		const double vn = DW(ii);

		if (ve!=IOUtils::nodata && vn!=IOUtils::nodata) {
			VW(ii) = sqrt(ve*ve + vn*vn);
			DW(ii) = fmod( atan2(ve,vn) * Cst::to_deg + 360., 360.);
		}
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
	grid.set(dem.getNx(), dem.getNy(), dem.cellsize, dem.llcorner, IOUtils::nodata);
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

void StandardPressureAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid) 
{
	Interpol2D::stdPressure(dem, grid);
}


double ConstAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);
	
	const size_t nr_args = vecArgs.size();
	if (nr_args!=1)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	
	IOUtils::convertString(user_cst, vecArgs[0]);
	return 0.01;
}

void ConstAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid) {
	Interpol2D::constant(user_cst, dem, grid);
}


double AvgAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 1) {
		return 0.8;
	} else if (nrOfMeasurments > 1) {
		return 0.2;
	}

	return 0.0;
}

void AvgAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid) {
	Interpol2D::constant(Interpol1D::arithmeticMean(vecData), dem, grid);
}


double AvgLapseRateAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
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
		// as good as IDW_LAPSE
		return 0.7;
	} else if (nrOfMeasurments>2){
		return 0.2;
	}

	return 0.2;
}

void AvgLapseRateAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
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

	if (nrOfMeasurments == 0) return 0.0;

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
                                               const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                      : InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), nrOfNeighbors(0)
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
	tsmanager.getMeteoData(date, vecMeteo);
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
	if ( (nrOfMeasurments<vecDataRH.size()/2) || ( nrOfMeasurments<2 ) )
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
		const double rh = vecDataRH[ii];
		if (rh<0. || rh>1.) {
			ostringstream ss;
			ss << "Invalid relative humidity: " << rh << " on " << date.toString(Date::ISO) << "\n";
			throw InvalidArgumentException(ss.str(), AT);
		}
		vecTd[ii] = Atmosphere::RhtoDewPoint(rh, vecDataTA[ii], 1);
	}

	Fit1D trend;
	getTrend(vecAltitudes, vecTd, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecTd);
	Interpol2D::IDW(vecTd, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);

	//Recompute Rh from the interpolated td
	for (size_t jj=0; jj<grid.getNy(); jj++) {
		for (size_t ii=0; ii<grid.getNx(); ii++) {
			double &value = grid(ii,jj);
			if (value!=IOUtils::nodata)
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
	for (size_t jj=0; jj<grid.getNy(); jj++) {
		for (size_t ii=0; ii<grid.getNx(); ii++) {
			double &value = grid(ii,jj);
			value = Atmosphere::blkBody_Radiation(value, ta(ii,jj));
		}
	}
}


double ListonWindAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for VW or DW
	if (in_param != MeteoData::VW && in_param != MeteoData::DW)
		return 0.0;

	date = i_date;
	param = in_param;
	vecData.clear(); vecMeta.clear();
	vecDataVW.clear(); vecDataDW.clear();

	nrOfMeasurments = 0;
	tsmanager.getMeteoData(date, vecMeteo);
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
	
	if ( (param==MeteoData::VW && Interpol2D::allZeroes(vecDataVW)) ||
	     (param==MeteoData::DW && Interpol2D::allZeroes(vecDataDW)) ) {
		inputIsAllZeroes = true;
		return 0.9;
	}

	if ( (nrOfMeasurments<vecDataVW.size()/2) || ( nrOfMeasurments<2 ) )
		return 0.6;

	return 0.9;
}

void ListonWindAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	if (param==MeteoData::VW) {
		Grid2DObject DW;
		simpleWindInterpolate(dem, vecDataVW, vecDataDW, grid, DW);
		Interpol2D::ListonWind(dem, grid, DW);
	}
	if (param==MeteoData::DW) {
		Grid2DObject VW;
		simpleWindInterpolate(dem, vecDataVW, vecDataDW, VW, grid);
		Interpol2D::ListonWind(dem, VW, grid);
	}
}


double RyanAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for VW or DW
	if (in_param != MeteoData::VW && in_param != MeteoData::DW)
		return 0.0;

	date = i_date;
	param = in_param;
	vecData.clear(); vecMeta.clear();
	vecDataVW.clear(); vecDataDW.clear();

	nrOfMeasurments = 0;
	tsmanager.getMeteoData(date, vecMeteo);
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
	
	if ( (param==MeteoData::VW && Interpol2D::allZeroes(vecDataVW)) ||
	     (param==MeteoData::DW && Interpol2D::allZeroes(vecDataDW)) ) {
		inputIsAllZeroes = true;
		return 0.9;
	}
	
	if (nrOfMeasurments<2)
		return 0.6;

	return 0.9;
}

void RyanAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	if (param==MeteoData::VW) {
		Grid2DObject DW;
		simpleWindInterpolate(dem, vecDataVW, vecDataDW, grid, DW);
		Interpol2D::RyanWind(dem, grid, DW);
	}
	if (param==MeteoData::DW) {
		Grid2DObject VW;
		simpleWindInterpolate(dem, vecDataVW, vecDataDW, VW, grid);
		Interpol2D::RyanWind(dem, VW, grid);
	}
}


const double WinstralAlgorithm::dmax = 300.;

WinstralAlgorithm::WinstralAlgorithm(Meteo2DInterpolator& i_mi, const std::vector<std::string>& i_vecArgs,
                                     const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                  : InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), base_algo("IDW_LAPSE"), ref_station(),
                    user_synoptic_bearing(IOUtils::nodata), inputIsAllZeroes(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args==1) {
		if (IOUtils::isNumeric(vecArgs[0]))
			IOUtils::convertString(user_synoptic_bearing, vecArgs[0]);
		else
			throw InvalidArgumentException("Please provide both the base interpolation method and the station_id to use for wind direction for the "+algo+" algorithm", AT);
		return;
	} else if (nr_args==2) {
		if (IOUtils::isNumeric(vecArgs[0])) {
			IOUtils::convertString(user_synoptic_bearing, vecArgs[0]);
			base_algo = IOUtils::strToUpper( vecArgs[1] );
		} else if (IOUtils::isNumeric(vecArgs[1]))  {
			IOUtils::convertString(user_synoptic_bearing, vecArgs[1]);
			base_algo = IOUtils::strToUpper( vecArgs[0] );
		} else {
			base_algo = IOUtils::strToUpper( vecArgs[0] );
			ref_station = vecArgs[1];
		}
		return;
	} else if (nr_args>2)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
}

double WinstralAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for PSUM (we could add HS later)
	if (in_param!=MeteoData::PSUM)
		return 0.0;

	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);
	inputIsAllZeroes = Interpol2D::allZeroes(vecData);

	if (inputIsAllZeroes)
		return 0.99;

	if (nrOfMeasurments==0)
		return 0.0;

	if (nrOfMeasurments==1 && ref_station.empty()) { //ie: still using default base_algo
		base_algo = "AVG";
	}

	//check that the necessary wind data is available
	if (user_synoptic_bearing==IOUtils::nodata) {
		if (!windIsAvailable(vecMeteo, ref_station))
			return 0.0;
	}

	return 0.99;
}

void WinstralAlgorithm::initGrid(const DEMObject& dem, Grid2DObject& grid)
{
	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	vector<string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, tsmanager, gridsmanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();
}

bool WinstralAlgorithm::windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station)
{
	if (ref_station.empty()) {
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			const double VW = vecMeteo[ii](MeteoData::VW);
			const double DW = vecMeteo[ii](MeteoData::DW);
			if (VW!=IOUtils::nodata && DW!=IOUtils::nodata)
				return true; //at least one station is enough
		}
	} else {
		if (getSynopticBearing(vecMeteo, ref_station) != IOUtils::nodata)
			return true;
	}

	return false;
}

double WinstralAlgorithm::getSynopticBearing(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station)
{
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (vecMeteo[ii].meta.stationID==ref_station)
			return vecMeteo[ii](MeteoData::DW);
	}

	return IOUtils::nodata;
}

//this method scans a square centered on the station for summmits
//that would shelter the station for wind.
//the sheltering criteria is: if ( height_difference*shade_factor > distance ) sheltered=true
bool WinstralAlgorithm::isExposed(const DEMObject& dem, Coords location)
{
	if (!dem.gridify(location)) {
		return false;
	}

	const int i_ref = location.getGridI();
	const int j_ref = location.getGridJ();
	const double alt_ref = location.getAltitude()+4.; //4 m mast added so flat terrain behaves properly
	const size_t ii_ref = static_cast<size_t>(i_ref);
	const size_t jj_ref = static_cast<size_t>(j_ref);

	const double shade_factor = 5.;
	const double cellsize = dem.cellsize;
	const double min_dh = cellsize/shade_factor; //min_dist=cellsize -> min_dh
	const double search_dist = (dem.grid2D.getMax() - alt_ref) * shade_factor;
	if (search_dist<=cellsize) return true;
	const int search_idx = static_cast<int>( Optim::ceil( search_dist/cellsize ) );

	const size_t i_min = static_cast<size_t>(std::max( 0, i_ref-search_idx ));
	const size_t i_max = static_cast<size_t>(std::min( static_cast<int>(dem.getNx()-1), i_ref+search_idx ));
	const size_t j_min = static_cast<size_t>(std::max( 0, j_ref-search_idx ));
	const size_t j_max = static_cast<size_t>(std::min( static_cast<int>(dem.getNy()-1), j_ref+search_idx ));

	for (size_t jj=j_min; jj<=j_max; jj++) {
		for (size_t ii=i_min; ii<=i_max; ii++) {
			if (ii==ii_ref && jj==jj_ref) continue; //skip the cell containing the station!

			const double dh = dem.grid2D(ii,jj) - alt_ref;
			if (dh<=min_dh) continue; //for negative or too small dh
			const double distance = Optim::fastSqrt_Q3( Optim::pow2(static_cast<int>(ii)-i_ref) + Optim::pow2(static_cast<int>(jj)-j_ref) ) * cellsize;
			if (distance<dh*shade_factor) {
				return false;
			}
		}
	}

	return true;
}

double WinstralAlgorithm::getSynopticBearing(const std::vector<MeteoData>& vecMeteo)
{
	double ve=0.0, vn=0.0;
	size_t count=0;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		const double VW = vecMeteo[ii](MeteoData::VW);
		const double DW = vecMeteo[ii](MeteoData::DW);
		if (VW!=IOUtils::nodata && DW!=IOUtils::nodata) {
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
		return static_cast<double>(Optim::round(meandirection/10.))*10.; //round to nearest 10 deg
	}

	return IOUtils::nodata;
}

double WinstralAlgorithm::getSynopticBearing(const DEMObject& dem, const std::vector<MeteoData>& vecMeteo)
{
	// 1) locate the stations in DEM and check if they are higher than their surroundings within a given radius
	// 2) simply compute a mean or median direction
	// (2) can be used on all the stations selected in (1)

	std::vector<MeteoData> stationsSubset;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (isExposed(dem, vecMeteo[ii].meta.position))
			stationsSubset.push_back( vecMeteo[ii] );
	}

	if (!stationsSubset.empty()) {
		return getSynopticBearing(stationsSubset);
	} else {
		//std::cerr << "[W] Synoptic wind direction computed from wind-sheltered stations only\n";
		return getSynopticBearing(vecMeteo);
	}
}

void WinstralAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	double synoptic_bearing = user_synoptic_bearing;
	if (synoptic_bearing==IOUtils::nodata) {
		if (!ref_station.empty())
			synoptic_bearing = getSynopticBearing(vecMeteo, ref_station);
		else
			synoptic_bearing = getSynopticBearing(dem, vecMeteo);
	}
	info << "DW=" << synoptic_bearing << " - ";
	initGrid(dem, grid);

	//get TA interpolation from call back to Meteo2DInterpolator
	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta);

	//alter the field with Winstral and the chosen wind direction
	Interpol2D::Winstral(dem, ta,  dmax, synoptic_bearing, grid);
}


const double WinstralListonAlgorithm::dmax = 300.;

WinstralListonAlgorithm::WinstralListonAlgorithm(Meteo2DInterpolator& i_mi, const std::vector<std::string>& i_vecArgs,
                                     const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                  : InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), base_algo("IDW_LAPSE"), ref_station(),
                    inputIsAllZeroes(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args==2) {
		base_algo = IOUtils::strToUpper( vecArgs[0] );
		ref_station = vecArgs[1];
		return;
	} else if (nr_args!=2)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
}

double WinstralListonAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for PSUM (we could add HS later)
	if (in_param!=MeteoData::PSUM)
		return 0.0;

	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);
	inputIsAllZeroes = Interpol2D::allZeroes(vecData);

	if (inputIsAllZeroes) return 0.99;

	if (nrOfMeasurments==0) return 0.0;

	if (nrOfMeasurments==1 && ref_station.empty()) { //ie: still using default base_algo
		base_algo = "AVG";
	}

	//check that the necessary wind data is available
	if (!ref_station.empty()) {
		if (!windIsAvailable(vecMeteo, ref_station))
			return 0.0;
	}

	return 0.99;
}

void WinstralListonAlgorithm::initGrid(const DEMObject& dem, Grid2DObject& grid)
{
	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	vector<string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, tsmanager, gridsmanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();
}

bool WinstralListonAlgorithm::windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station)
{
	if (ref_station.empty()) {
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			const double VW = vecMeteo[ii](MeteoData::VW);
			const double DW = vecMeteo[ii](MeteoData::DW);
			if (VW!=IOUtils::nodata && DW!=IOUtils::nodata)
				return true; //at least one station is enough
		}
	} else {
		double VW, DW;
		getSynopticWind(vecMeteo, ref_station, VW, DW);
		if (VW!=IOUtils::nodata && DW!=IOUtils::nodata)
			return true;
	}

	return false;
}

void WinstralListonAlgorithm::getSynopticWind(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station, double& VW, double& DW)
{
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (vecMeteo[ii].meta.stationID==ref_station) {
			VW = vecMeteo[ii](MeteoData::VW);
			DW = vecMeteo[ii](MeteoData::DW);
			return;
		}
	}
}

void WinstralListonAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	//get initial PSUM grid
	initGrid(dem, grid);

	//get meteo fields interpolation from call back to Meteo2DInterpolator
	Grid2DObject ta, dw, vw;
	mi.interpolate(date, dem, MeteoData::TA, ta);
	mi.interpolate(date, dem, MeteoData::DW, dw);
	mi.interpolate(date, dem, MeteoData::VW, vw);

	//alter the field with Winstral and the chosen wind direction
	Interpol2D::Winstral(dem, ta, dw, vw, dmax, grid);
}


std::string USERInterpolation::getGridFileName() const
{
	const size_t nrArgs = vecArgs.size();
	if (nrArgs > 2){
		throw InvalidArgumentException("Too many arguments for the "+algo+" interpolation algorithm", AT);
	}
	const std::string prefix = (nrArgs==1)? vecArgs[0] + "/" : "";
	const std::string ext = (nrArgs==2)? vecArgs[1] : ".asc";
	const std::string gridname =  prefix + date.toString(Date::NUM) + "_" + MeteoData::getParameterName(param) + ext;

	return gridname;
}

double USERInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	filename = getGridFileName();
	
	if (grid2d_path.empty()) {
		const Config cfg( gridsmanager.getConfig() );
		cfg.getValue("GRID2DPATH", "Input", grid2d_path);
	}

	if (!FileUtils::validFileAndPath(grid2d_path+"/"+filename)) {
		cerr << "[E] Invalid grid filename for "+algo+" interpolation algorithm: " << grid2d_path+"/"+filename << "\n";
		return 0.0;
	}
	
	return (FileUtils::fileExists(grid2d_path+"/"+filename))? 1. : 0.;
}

void USERInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	gridsmanager.read2DGrid(grid, filename);
	if (!grid.isSameGeolocalization(dem)) {
		throw InvalidArgumentException("[E] trying to load a grid(" + filename + ") that does not have the same georeferencing as the DEM!", AT);
	} else {
		info << FileUtils::getFilename(filename);
	}
}

ALS_Interpolation::ALS_Interpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
		  : InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), ALS_scan(), filename(), 
		    grid2d_path(), base_algo(), base_algo_user(), ta_thresh(IOUtils::nodata), als_mean(IOUtils::nodata), inputIsAllZeroes(false) 
{
	const Config cfg( gridsmanager.getConfig() );
	cfg.getValue("GRID2DPATH", "Input", grid2d_path);
	
	const size_t nr_args = vecArgs.size();
	if (nr_args<2 || nr_args>3)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);

	for (size_t ii=0; ii<nr_args; ii++) {
		if (IOUtils::isNumeric(vecArgs[ii])) {
			if (ta_thresh!=IOUtils::nodata)
				throw InvalidArgumentException("Only one threshold temperature can be provided to the "+algo+" algorithm", AT);
			IOUtils::convertString(ta_thresh, vecArgs[ii]);
		} else {
			if (base_algo_user.empty())
				base_algo_user = IOUtils::strToUpper( vecArgs[ii] );
			else if (filename.empty())
				filename = vecArgs[ii];
			else 
				throw InvalidArgumentException("Wrong arguments for the "+algo+" algorithm: ta_thresh, base_algo and filename are required", AT);
		}
	}

	if (!FileUtils::validFileAndPath(grid2d_path+"/"+filename)) {
		throw InvalidNameException("[E] Invalid grid filename for "+algo+" interpolation algorithm: "+grid2d_path+"/"+filename, AT);
	}
}


void ALS_Interpolation::initGrid(const DEMObject& dem, Grid2DObject& grid)
{
	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	vector<string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, tsmanager, gridsmanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();
}

double ALS_Interpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);
	inputIsAllZeroes = Interpol2D::allZeroes(vecData);

	if (nrOfMeasurments==0) return 0.;
	if (inputIsAllZeroes) return 1.;
	
	base_algo = (nrOfMeasurments>1)? base_algo_user : "AVG";
	
	return (FileUtils::fileExists(grid2d_path+"/"+filename))? 1. : 0.;
}

void ALS_Interpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	
	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}
	
	if (ALS_scan.empty()) { //read the ALS scan if necessary
		gridsmanager.read2DGrid(ALS_scan, filename);
		if (ta_thresh==IOUtils::nodata) {
			als_mean = ALS_scan.grid2D.getMean();
			if (als_mean==0. || als_mean==IOUtils::nodata)
				throw InvalidArgumentException("[E] the scaling grid(" + filename + ") can not have a nul or nodata mean for the '"+algo+"' method!", AT);
		}
	}
	
	//check that the ALS scan matches the provided DEM
	if (!ALS_scan.isSameGeolocalization(dem))
		throw InvalidArgumentException("[E] trying to load a grid(" + filename + ") that does not have the same georeferencing as the DEM!", AT);
	else
		info << FileUtils::getFilename(filename) << " - ";
	
	initGrid(dem, grid);
	const size_t nxy = grid.getNx()*grid.getNy();
	
	// create map of TA to differ between solid and liquid precipitation
	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator
	
	//Compute scaling factors: psum_mean and als_mean according to the pixels selection criteria
	double psum_mean = 0.;
	if (ta_thresh==IOUtils::nodata) { //simple case: no TA_THRESH
		double psum_sum = 0.;
		size_t count = 0;
		for (size_t jj=0; jj<nxy; jj++) {
			const double val = grid(jj);
			const bool has_Scan = (ALS_scan(jj)!=IOUtils::nodata);
			if (val!=IOUtils::nodata && has_Scan ) {
				psum_sum += val;
				count++;
			}
		}
		if (count==0) return; //no overlap between ALS, TA and initial grid
		psum_mean = psum_sum / static_cast<double>( count );
	} else { //use local air temperature
		double psum_sum = 0., als_sum = 0.;
		size_t count = 0;
		for (size_t jj=0; jj<nxy; jj++) {
			const double val = grid(jj);
			const bool has_Scan = (ALS_scan(jj)!=IOUtils::nodata);
			const bool has_TA = (ta(jj)!=IOUtils::nodata);
			if (val!=IOUtils::nodata && has_Scan && has_TA && ta(jj) < ta_thresh) {
				psum_sum += val;
				als_sum += ALS_scan(jj);
				count++;
			}
		}
		if (count==0) return; //no overlap between ALS, TA and initial grid
		psum_mean = psum_sum / static_cast<double>( count );
		als_mean = als_sum / static_cast<double>( count );
		if (als_mean==0.)
			throw InvalidArgumentException("[E] the scaling grid(" + filename + ") can not have a nul mean for the '"+algo+"' method!", AT);
	}
	
	//pixels that are nodata are kept such as computed by "base_algo", otherwise we take the newly computed values
	if (ta_thresh==IOUtils::nodata) { //simple case: no TA_THRESH
		for (size_t jj=0; jj<nxy; jj++) {
			double &val = grid(jj);
			const bool has_Scan = (ALS_scan(jj)!=IOUtils::nodata);
			if (val!=IOUtils::nodata && has_Scan) 
				val = ALS_scan(jj) * ( psum_mean / als_mean );
		}
	} else { //use local air temperature
		for (size_t jj=0; jj<nxy; jj++) {
			double &val = grid(jj);
			const bool has_Scan = (ALS_scan(jj)!=IOUtils::nodata);
			const bool has_TA = (ta(jj)!=IOUtils::nodata);
			if (val!=IOUtils::nodata && has_Scan && has_TA && ta(jj) < ta_thresh)
				val = ALS_scan(jj) * ( psum_mean / als_mean );
		}
	}
}


double PPHASEInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	
	const size_t nArgs = vecArgs.size();
	
	if (nArgs<1 || IOUtils::isNumeric(vecArgs[0]))
		throw InvalidArgumentException("Wrong arguments supplied to the "+algo+" interpolation. Please provide the method to use and its arguments!", AT);
	
	const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);
	if (user_algo=="THRESH") {
		if (nArgs!=2)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" interpolation for the "+user_algo+" method", AT);
		IOUtils::convertString(fixed_thresh, vecArgs[1]);
		model = THRESH;
	} else if (user_algo=="RANGE") {
		if (nArgs!=3)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" interpolation for the "+user_algo+" method", AT);
		double range_thresh1, range_thresh2;
		IOUtils::convertString(range_thresh1, vecArgs[1]);
		IOUtils::convertString(range_thresh2, vecArgs[2]);
		if (range_thresh1==range_thresh2)
			throw InvalidArgumentException(algo+" interpolation, "+user_algo+" method: the two provided threshold must be different", AT);
		if (range_thresh1>range_thresh2) 
			std::swap(range_thresh1, range_thresh2);
		range_start = range_thresh1;
		range_norm = 1. / (range_thresh2-range_thresh1);
		model = RANGE;
	} else
		throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied to the "+algo+" interpolation", AT);
	
	return 0.1;
}

void PPHASEInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	
	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator
	
	grid.set(dem, IOUtils::nodata);
	const size_t nrCells = dem.getNx() * dem.getNy();
	
	if (model==THRESH) {
		for (size_t ii=0; ii<nrCells; ii++) {
			const double TA=ta(ii);
			if (TA==IOUtils::nodata) continue;
			grid(ii) = (TA>=fixed_thresh)? 1. : 0.;
		}
	} else if (model==RANGE) {
		for (size_t ii=0; ii<nrCells; ii++) {
			const double TA=ta(ii);
			if (TA==IOUtils::nodata) continue;
			const double tmp_rainfraction = range_norm * (TA - range_start);
			grid(ii) = (tmp_rainfraction>1)? 1. : (tmp_rainfraction<0.)? 0. : tmp_rainfraction;
		}
	}
}


double SnowPSUMInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.9;
}

void SnowPSUMInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//retrieve optional arguments
	std::string base_algo("IDW_LAPSE");
	const size_t nrArgs = vecArgs.size();
	if (nrArgs == 1){
		IOUtils::convertString(base_algo, vecArgs[0]);
	} else if (nrArgs>1){ //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	}

	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	IOUtils::toUpper(base_algo);
	vector<string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, tsmanager, gridsmanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();

	//get TA interpolation from call back to Meteo2DInterpolator
	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta);

	//slope/curvature correction for solid precipitation
	const double orig_mean = grid.grid2D.getMean();
	Interpol2D::PrecipSnow(dem, ta, grid);
	//HACK: correction for precipitation sum over the whole domain
	//this is a cheap/crappy way of compensating for the spatial redistribution of snow on the slopes
	const double new_mean = grid.grid2D.getMean();
	if (new_mean!=0.) grid *= orig_mean/new_mean;
	
	//Interpol2D::SteepSlopeRedistribution(dem, ta, grid);
	//Interpol2D::CurvatureCorrection(dem, ta, grid);
}


const double SWRadInterpolation::soil_albedo = .23; //grass
const double SWRadInterpolation::snow_albedo = .85; //snow
const double SWRadInterpolation::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo
double SWRadInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	
	const size_t nrArgs = vecArgs.size();
	if ( nrArgs>1 ) {
		throw InvalidArgumentException("Wrong arguments supplied to the "+algo+" interpolation.", AT);
	} else if (nrArgs==1) {
		if (IOUtils::strToUpper(vecArgs[0])=="NO_SHADING") {
			shading = false;
		} else
			throw InvalidArgumentException("Unknown argument \""+vecArgs[0]+"\" supplied to the "+algo+" interpolation", AT);
	}
	
	tsmanager.getMeteoData(i_date, vecMeteo);
	
	//fill vecIdx with the indices of the stations that can be used and set the Sun coordinates to the middle of the stations
	double avg_lat = 0., avg_lon = 0., avg_alt = 0.;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		const Coords &location( vecMeteo[ii].meta.position );
		const bool has_meta = (location.getLat()!=IOUtils::nodata) && (location.getLon()!=IOUtils::nodata) && (location.getAltitude()!=IOUtils::nodata);
		const bool has_meteo = (vecMeteo[ii](MeteoData::ISWR)!=IOUtils::nodata) 
		                                  && (vecMeteo[ii](MeteoData::TA)!=IOUtils::nodata) 
		                                  && (vecMeteo[ii](MeteoData::RH)!=IOUtils::nodata);
		if (has_meta && has_meteo) {
			vecMeta.push_back( vecMeteo[ii].meta );
			avg_lat += location.getLat();
			avg_lon += location.getLon();
			avg_alt += location.getAltitude();
			vecIdx.push_back( ii );
		}
	}
	
	nrOfMeasurments = vecIdx.size();
	if (nrOfMeasurments == 0) return 0.0;

	avg_lat /= static_cast<double>(nrOfMeasurments);
	avg_lon /= static_cast<double>(nrOfMeasurments);
	avg_alt /= static_cast<double>(nrOfMeasurments);
	Sun.setLatLon(avg_lat, avg_lon, avg_alt);
	Sun.setDate(vecMeteo[ vecIdx[0] ].date.getJulian(), vecMeteo[ vecIdx[0] ].date.getTimeZone()); //we have at least one station
	
	return 0.9;
}

void SWRadInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	double solarAzimuth, solarElevation;
	Sun.position.getHorizontalCoordinates(solarAzimuth, solarElevation);
	const double tan_sun_elev = tan(solarElevation*Cst::to_rad);
	
	//compute the correction factors at every station
	std::vector<double> vecMd(nrOfMeasurments, IOUtils::nodata);
	std::vector<double> vecCorr(nrOfMeasurments, IOUtils::nodata);
	bool glob_day=true, glob_night=true;
	for (size_t ii=0; ii<nrOfMeasurments; ii++) {
		const size_t idx = vecIdx[ii];
		const double HS = vecMeteo[idx](MeteoData::HS);
		double albedo = 0.5;
		if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
			albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
		Sun.calculateRadiation(vecMeteo[idx](mio::MeteoData::TA), vecMeteo[idx](mio::MeteoData::RH), vecMeteo[idx](mio::MeteoData::P), albedo);

		bool day, night;
		vecCorr[ii] = Sun.getCorrectionFactor(vecMeteo[idx](mio::MeteoData::ISWR), vecMd[ii], day, night);
		glob_day = glob_day && day;
		glob_night = glob_night && night;
	}
	
	if (glob_night) {
		Interpol2D::constant(0., dem, grid); //no iswr at night
		return;
	}
	
	//compute the distributed  splitting and correction coefficient fields
	Grid2DObject Md;
	Interpol2D::IDW(vecMd, vecMeta, dem, Md);
	Grid2DObject Corr;
	Interpol2D::IDW(vecCorr, vecMeta, dem, Corr);
	
	//get TA, RH and P interpolation from call back to Meteo2DInterpolator
	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta);
	Grid2DObject rh;
	mi.interpolate(date, dem, MeteoData::RH, rh);
	Grid2DObject p;
	mi.interpolate(date, dem, MeteoData::P, p);
	
	//fill the final results with the proper radiation (with shading)
	grid.set(dem, IOUtils::nodata);
	for (size_t jj=0; jj<dem.getNy(); jj++) {
		for (size_t ii=0; ii<dem.getNx(); ii++) {
			if (dem(ii,jj)==IOUtils::nodata) continue;
			
			Sun.resetAltitude( dem(ii,jj) );
			Sun.calculateRadiation(ta(ii,jj), rh(ii,jj), p(ii,jj), .5); //we don't have any albedo, so use .5
			double cell_toa, cell_direct, cell_diffuse;
			Sun.getHorizontalRadiation(cell_toa, cell_direct, cell_diffuse);
			
			if (glob_day && shading) { //at dawn/dusk, we consider it to be all diffuse, so no shading
				const double tan_horizon = dem.getHorizon(ii, jj, solarAzimuth);
				
				//redo the splitting using the distributed splitting coefficient
				const double global = cell_direct + cell_diffuse;
				cell_direct = global * Md(ii,jj);
				cell_diffuse = global * Md(ii,jj);

				if ( tan_sun_elev<tan_horizon ) cell_direct = 0.;//cell is shaded
			}
			grid(ii,jj) = Corr(ii,jj) * (cell_direct+cell_diffuse);
		}
	}
}


//since this uses vecData, the optional detrending has already been done
void OrdinaryKrigingAlgorithm::getDataForEmpiricalVariogram(std::vector<double> &distData, std::vector<double> &variData) const
{
	distData.clear();
	variData.clear();

	for (size_t j=0; j<nrOfMeasurments; j++) {
		const Coords& st1 = vecMeta[j].position;
		const double x1 = st1.getEasting();
		const double y1 = st1.getNorthing();
		const double val1 = vecData[j];

		for (size_t i=0; i<j; i++) {
			//compute distance between stations
			const Coords& st2 = vecMeta[i].position;
			const double val2 = vecData[i];
			const double DX = x1-st2.getEasting();
			const double DY = y1-st2.getNorthing();
			//const double distance = Optim::fastSqrt_Q3( Optim::pow2(DX) + Optim::pow2(DY) );
			const double inv_distance = Optim::invSqrt( Optim::pow2(DX) + Optim::pow2(DY) );

			distData.push_back( 1./inv_distance );
			variData.push_back( inv_distance*Optim::pow2(val1-val2) );
		}
	}

	if (distData.size()>40)
		Interpol1D::equalCountBin(10, distData, variData);
}

size_t OrdinaryKrigingAlgorithm::getTimeSeries(const bool& detrend_data, std::vector< std::vector<double> > &vecVecData) const
{
	//get all the data between "date" and "date-daysBefore"
	const double daysBefore = 1.5;
	const double Tstep = 1./24.;
	Date d1 = date - daysBefore;

	vecVecData.clear();
	std::vector<MeteoData> Meteo;
	tsmanager.getMeteoData(d1, Meteo);
	const size_t nrStations = Meteo.size();
	vecVecData.insert(vecVecData.begin(), nrStations, std::vector<double>()); //allocation for the vectors

	//get the stations altitudes
	vector<double> vecAltitudes;
	for (size_t ii=0; ii<nrStations; ii++){
		const double& alt = Meteo[ii].meta.position.getAltitude();
		if (alt != IOUtils::nodata) {
			vecAltitudes.push_back(alt);
		}
	}

	//fill time series
	for (; d1<=date; d1+=Tstep) {
		tsmanager.getMeteoData(d1, Meteo);
		if (Meteo.size()!=nrStations) {
			std::ostringstream ss;
			ss << "Number of stations varying between " << nrStations << " and " << Meteo.size();
			ss << ". This is currently not supported for " << algo << "!";
			throw InvalidArgumentException(ss.str(), AT);
		}

		if (detrend_data) { //detrend data
			//find trend
			std::vector<double> vecDat;
			for (size_t ii=0; ii<nrStations; ii++)
				vecDat.push_back( Meteo[ii](param) );

			Fit1D trend;
			trend.setModel(Fit1D::NOISY_LINEAR, vecAltitudes, vecDat, false);
			const bool status = trend.fit();
			if (!status)
				throw InvalidArgumentException("Could not fit variogram model to the data", AT);

			//detrend the data
			for (size_t ii=0; ii<nrStations; ii++) {
				double val = Meteo[ii](param);
				if (val!=IOUtils::nodata)
					val -= trend( vecAltitudes[ii] );

				vecVecData.at(ii).push_back( val );
			}
		} else { //do not detrend data
			for (size_t ii=0; ii<nrStations; ii++)
				vecVecData.at(ii).push_back( Meteo[ii](param) );
		}
	}

	return nrStations;
}

//this gets the full data over a preceeding period
//we can not rely on the data in vecData/vecMeta since they filter nodata points
//and we need to do our own nodata handling here (to keep stations' indices constant)
void OrdinaryKrigingAlgorithm::getDataForVariogram(std::vector<double> &distData, std::vector<double> &variData, const bool& detrend_data) const
{
	distData.clear();
	variData.clear();

	std::vector< std::vector<double> > vecTimeSeries;
	const size_t nrStations = getTimeSeries(detrend_data, vecTimeSeries);

	//for each station, compute distance to other stations and
	// variance of ( Y(current) - Y(other station) )
	for (size_t j=0; j<nrStations; j++) {
		const Coords& st1 = vecMeta[j].position;
		const double x1 = st1.getEasting();
		const double y1 = st1.getNorthing();

		for (size_t i=0; i<j; i++) { //compare with the other stations
			//compute distance between stations
			const Coords& st2 = vecMeta[i].position;
			const double DX = x1-st2.getEasting();
			const double DY = y1-st2.getNorthing();
			const double distance = Optim::fastSqrt_Q3( Optim::pow2(DX) + Optim::pow2(DY) );

			std::vector<double> Y;
			for (size_t dt=0; dt<vecTimeSeries[j].size(); dt++) {
				if (vecTimeSeries[j][dt]!=IOUtils::nodata && vecTimeSeries[i][dt]!=IOUtils::nodata)
					Y.push_back( vecTimeSeries[j][dt] - vecTimeSeries[i][dt] );
			}

			distData.push_back( distance );
			variData.push_back( Interpol1D::variance(Y) );
		}
	}

	if (distData.size()>40)
		Interpol1D::equalCountBin(10, distData, variData);
}

bool OrdinaryKrigingAlgorithm::computeVariogram(const bool& /*detrend_data*/)
{//return variogram fit of covariance between stations i and j
	std::vector<double> distData, variData;
	getDataForEmpiricalVariogram(distData, variData);
	//getDataForVariogram(distData, variData, detrend_data);

	std::vector<string> vario_types( vecArgs );
	if (vario_types.empty()) vario_types.push_back("LINVARIO");

	size_t args_index=0;
	do {
		const string vario_model = IOUtils::strToUpper( vario_types[args_index] );
		const bool status = variogram.setModel(vario_model, distData, variData);
		if (status) {
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

	if (nrOfMeasurments==0) return 0.;
	if (nrOfMeasurments>=7) return 0.9;
	return 0.1;
}

void OrdinaryKrigingAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//optimization: getRange (from variogram fit -> exclude stations that are at distances > range (-> smaller matrix)
	//or, get max range from io.ini, build variogram from this user defined max range
	if (!computeVariogram(false)) //only refresh once a month, or once a week, etc
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
	if (!trend.fit())
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param) + ": " + trend.getInfo(), AT);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);

	if (!computeVariogram(true)) //only refresh once a month, or once a week, etc
		throw IOException("The variogram for parameter " + MeteoData::getParameterName(param) + " could not be computed!", AT);
	Interpol2D::ODKriging(vecData, vecMeta, dem, variogram, grid);

	retrend(dem, trend, grid);
}

} //namespace
