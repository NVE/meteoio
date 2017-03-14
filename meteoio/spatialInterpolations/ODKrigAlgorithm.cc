/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include <meteoio/spatialInterpolations/ODKrigAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/MathOptim.h>

namespace mio {

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
	std::vector<double> vecAltitudes;
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

	std::vector<std::string> vario_types( vecArgs );
	if (vario_types.empty()) vario_types.push_back("LINVARIO");

	size_t args_index=0;
	do {
		const std::string vario_model( IOUtils::strToUpper( vario_types[args_index] ) );
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

} //namespace