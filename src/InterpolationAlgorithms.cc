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
#include "InterpolationAlgorithms.h"

using namespace std;

namespace mio {

std::set<std::string> AlgorithmFactory::setAlgorithms;
const bool AlgorithmFactory::__init = AlgorithmFactory::initStaticData();

bool AlgorithmFactory::initStaticData()
{
	/* 
	 * Keywords for selecting the spatial interpolation algorithm among the 
	 * available methods for single source and multiple sources interpolations. 
	 * More details about some of these algorithms can be found in "A Meteorological 
	 * Distribution System for High-Resolution Terrestrial Modeling (MicroMet)", Liston and Elder, 2006.
	 * Please don't forget to document InterpolationAlgorithms.h for the available algorithms!
	 */

	setAlgorithms.insert("CST");       // constant fill
	setAlgorithms.insert("STD_PRESS"); // standard air pressure interpolation
	setAlgorithms.insert("CST_LAPSE"); // constant fill with an elevation lapse rate
	setAlgorithms.insert("IDW");       // Inverse Distance Weighting fill
	setAlgorithms.insert("IDW_LAPSE"); // Inverse Distance Weighting with an elevation lapse rate fill
	setAlgorithms.insert("RH");        // relative humidity interpolation
	setAlgorithms.insert("WIND_CURV"); // wind velocity interpolation (using a heuristic terrain effect)
	setAlgorithms.insert("USER"); // read user provided grid

	return true;
}

InterpolationAlgorithm* AlgorithmFactory::getAlgorithm(const std::string& _algoname, 
                                                       const Meteo2DInterpolator& _mi,
                                                       const DEMObject& _dem,
                                                       const std::vector<MeteoData>& _vecMeteo,
                                                       const std::vector<StationData>& _vecStation,
                                                       const std::vector<std::string>& _vecArgs)
{
	std::string algoname(_algoname);
	IOUtils::toUpper(algoname);

	//Check whether algorithm theoretically exists
	if (setAlgorithms.find(algoname) == setAlgorithms.end())
		throw UnknownValueException("The interpolation algorithm '"+algoname+"' does not exist" , AT);

	if (algoname == "CST"){
		return new ConstAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "STD_PRESS"){
		return new StandardPressureAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "CST_LAPSE"){
		return new ConstLapseRateAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "IDW"){
		return new IDWAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "IDW_LAPSE"){
		return new IDWLapseAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "RH"){
		return new RHAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "WIND_CURV"){
		return new SimpleWindInterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else if (algoname == "USER"){
		return new USERinterpolation(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algoname);
	} else {
		throw IOException("The interpolation algorithm '"+algoname+"' is not implemented" , AT);
	}

	return NULL;
}

unsigned int InterpolationAlgorithm::getData(const MeteoData::Parameters& param, 
                                             std::vector<double>& vecData) const
{
	vector<StationData> vecMeta;
	return getData(param, vecData, vecMeta);
	
}

unsigned int InterpolationAlgorithm::getData(const MeteoData::Parameters& param, 
                                             std::vector<double>& vecData, std::vector<StationData>& vecMeta) const
{
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		const double& val = vecMeteo[ii].param(param);
		if (val != IOUtils::nodata){
			vecData.push_back(val);
			vecMeta.push_back(vecStation[ii]);
		}
	}

	return vecData.size();
}

unsigned int InterpolationAlgorithm::getStationAltitudes(const std::vector<StationData>& vecMeta,
                                                         std::vector<double>& vecData) const
{
	for (unsigned int ii=0; ii<vecMeta.size(); ii++){
		const double& val = vecMeta[ii].position.getAltitude();
		if (val != IOUtils::nodata)
			vecData.push_back(val);
	}

	return vecData.size();
}

void InterpolationAlgorithm::printInfo(const MeteoData::Parameters& param, const unsigned int& nrOfMeasurments) const
{
	std::cout << "[i] Interpolating " << MeteoData::getParameterName(param) << " (" << algo << ") ";
	std::cout << "on " << nrOfMeasurments << " stations" << std::endl;
}

/**********************************************************************************/

double ConstAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	const unsigned int nrOfMeasurments = getData(param, vecData);

	if (nrOfMeasurments == 0){
		return 0.0;
	} else if (nrOfMeasurments == 1){
		return 0.8;
	} else if (nrOfMeasurments > 1){
		return 0.2;
	}

	return 0.0;
}


void ConstAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{
	//get a vector of all the relevant data
	vector<double> vecData;
	const unsigned int nrOfMeasurments = getData(param, vecData);
	printInfo(param, nrOfMeasurments);

	//check how many data points there are (not including nodata)
	if (nrOfMeasurments == 0)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	//run algorithm
	Interpol2D::constantGrid2DFill(Interpol1D::arithmeticMean(vecData), dem, grid);
}


double StandardPressureAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	if (param != MeteoData::P)
		return 0.0;

	//check incoming data for number of data points
	vector<double> vecData;
	const unsigned int nrOfMeasurments = getData(param, vecData);

	if (nrOfMeasurments == 0)
		return 1.0;

	return 0.1;
}

void StandardPressureAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{
	printInfo(param, 0);
	//run algorithm
	Interpol2D::stdPressureGrid2DFill(dem, grid);
}


double ConstLapseRateAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	const unsigned int nrOfMeasurments = getData(param, vecData);

	if (nrOfMeasurments == 0){
		return 0.0;
	} else if (nrOfMeasurments == 1){
		return 0.9;
	} else if (nrOfMeasurments>1){
		return 0.2;
	}

	return 0.2;
}

void ConstLapseRateAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	vector<double> vecAltitudes, vecData;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = getData(param, vecData, vecMeta);
	printInfo(param, nrOfMeasurments);

	getStationAltitudes(vecMeta, vecAltitudes);
	LapseRateProjectPtr funcptr = &Interpol2D::LinProject;	

	if ((vecAltitudes.size() == 0) || (nrOfMeasurments == 0))
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	double avgAltitudes = Interpol1D::arithmeticMean(vecAltitudes);
	double avgData = Interpol1D::arithmeticMean(vecData);

	if (vecData.size() == 1)
		funcptr = &Interpol2D::LinProject;

	//Set regression coefficients
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);

	//Get the argument for the algorithm: the default lapse rate, otherwise throw an exception
	if (vecArgs.size() == 1){
		IOUtils::convertString(vecCoefficients[1], vecArgs[0]);
	} else {
		throw InvalidArgumentException("Wrong number of arguments supplied for the CST_LAPSE algorithm", AT);
	}

	//run algorithm
	Interpol2D::constantLapseGrid2DFill(avgData, avgAltitudes, dem, vecCoefficients, funcptr, grid);
}


double IDWAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	const unsigned int nrOfMeasurments = getData(param, vecData);

	if (nrOfMeasurments == 0){
		return 0.0;
	} else if (nrOfMeasurments == 1){
		return 0.3;
	} else if (nrOfMeasurments > 1){
		return 0.5;
	}

	return 0.2;
}

void IDWAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	vector<double> vecData;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = getData(param, vecData, vecMeta);
	printInfo(param, nrOfMeasurments);

	if (nrOfMeasurments == 0)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	//run algorithm
	Interpol2D::IDW(vecData, vecMeta, dem, grid);
}

double IDWLapseAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	const unsigned int nrOfMeasurments = getData(param, vecData);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.7;
}

void IDWLapseAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	vector<double> vecAltitudes, vecData;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = getData(param, vecData, vecMeta);
	printInfo(param, nrOfMeasurments);

	if (nrOfMeasurments == 0)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	//Set regression coefficients
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);

	getStationAltitudes(vecMeta, vecAltitudes);

	//run algorithm
	Interpol2D::LinRegression(vecAltitudes, vecData, vecCoefficients);
	Interpol2D::LapseIDW(vecData, vecMeta, dem, vecCoefficients, &Interpol2D::LinProject, grid);
}

double RHAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//This algorithm is only valid for RH
	if (param != MeteoData::RH)
		return 0.0;

	//check incoming data for number of data points

	vector<double> vecDataTA, vecDataRH;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = 0;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].rh != IOUtils::nodata) && (vecMeteo[ii].ta != IOUtils::nodata)){
			vecDataTA.push_back(vecMeteo[ii].ta);
			vecDataRH.push_back(vecMeteo[ii].rh);
			vecMeta.push_back(vecStation[ii]);
			nrOfMeasurments++;
		}
	}

	if (vecDataTA.size() == 0)
		return 0.0;
	if( ( nrOfMeasurments<(unsigned int)(0.5*vecDataRH.size()) ) || ( nrOfMeasurments<2 ) )
		return 0.6;

	return 0.9;
}

void RHAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	//This algorithm is only valid for RH
	if (param != MeteoData::RH)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	vector<double> vecAltitudes, vecDataTA, vecDataRH;
	vector<StationData> vecMeta;

	unsigned int nrOfMeasurments = 0;
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].rh != IOUtils::nodata) && (vecMeteo[ii].ta != IOUtils::nodata)){
			vecDataTA.push_back(vecMeteo[ii].ta);
			vecDataRH.push_back(vecMeteo[ii].rh);
			vecMeta.push_back(vecStation[ii]);
			nrOfMeasurments++;
		}
	}
	printInfo(param, nrOfMeasurments);
	getStationAltitudes(vecMeta, vecAltitudes);

	if (vecDataTA.size() == 0) //No matching data
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	Grid2DObject ta;
	mi.interpolate(MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator

	//here, RH->Td, interpolations, Td->RH
	std::vector<double> vecTd(vecDataRH.size(), 0.0); // init to 0.0

	//Compute dew point temperatures at stations
	for (unsigned int ii=0; ii<vecDataRH.size(); ii++){
		vecTd[ii] = Interpol2D::RhtoDewPoint(vecDataRH[ii], vecDataTA[ii], 1);
	}
			
	//Krieging on Td
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);

	//run algorithm
	Interpol2D::LinRegression(vecAltitudes, vecTd, vecCoefficients);
	Interpol2D::LapseIDW(vecTd, vecMeta, dem, vecCoefficients, &Interpol2D::LinProject, grid);

	//Recompute Rh from the interpolated td
	for (unsigned int ii=0; ii<grid.ncols; ii++) {
		for (unsigned int jj=0; jj<grid.nrows; jj++) {
			grid.grid2D(ii,jj) = Interpol2D::DewPointtoRh(grid.grid2D(ii,jj), ta.grid2D(ii,jj), 1);
		}
	}
}

double SimpleWindInterpolationAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//This algorithm is only valid for VW
	if (param != MeteoData::VW)
		return 0.0;

	vector<double> vecAltitudes, vecDataVW, vecDataDW;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = 0;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].vw != IOUtils::nodata) && (vecMeteo[ii].dw != IOUtils::nodata)){
			vecDataVW.push_back(vecMeteo[ii].vw);
			vecDataDW.push_back(vecMeteo[ii].dw);
			vecMeta.push_back(vecStation[ii]);
			nrOfMeasurments++;
		}
	}
	
	if (vecDataVW.size() == 0)
		return 0.0;
	if( ( nrOfMeasurments<(unsigned int)(0.5*vecDataVW.size()) ) || ( nrOfMeasurments<2 ) )
		return 0.6;

	return 0.9;
}

void SimpleWindInterpolationAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	//This algorithm is only valid for VW
	if (param != MeteoData::VW)
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	vector<double> vecAltitudes, vecDataVW, vecDataDW;
	vector<StationData> vecMeta;

	unsigned int nrOfMeasurments = 0;
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].vw != IOUtils::nodata) && (vecMeteo[ii].dw != IOUtils::nodata)){
			vecDataVW.push_back(vecMeteo[ii].vw);
			vecDataDW.push_back(vecMeteo[ii].dw);
			vecMeta.push_back(vecStation[ii]);
			nrOfMeasurments++;
		}
	}
	printInfo(param, nrOfMeasurments);
	getStationAltitudes(vecMeta, vecAltitudes);

	if( vecDataDW.size() == 0) 
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param), AT);

	Grid2DObject dw;
	mi.interpolate(MeteoData::DW, dw); //get DW interpolation from call back to Meteo2DInterpolator
	
	//Krieging
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);
	
	Interpol2D::LinRegression(vecAltitudes, vecDataVW, vecCoefficients);
	Interpol2D::LapseIDW(vecDataVW, vecMeta, dem, vecCoefficients, &Interpol2D::LinProject, grid);
	Interpol2D::SimpleDEMWindInterpolate(dem, grid, dw);
}

std::string USERinterpolation::getGridFileName(const MeteoData::Parameters& param)
{
	const std::string ext=std::string(".asc");
	if (vecArgs.size() != 1){
		throw InvalidArgumentException("Please provide the path to the grids for the USER interpolation algorithm", AT);
	}
	const std::string& grid_path = vecArgs[0];
	std::string gridname = grid_path + std::string("/");

	if(vecMeteo.size()>0) {
		const Date& timestep = vecMeteo.at(0).date;
		gridname = gridname + timestep.toString(Date::NUM) + std::string("_") + MeteoData::getParameterName(param) + ext;
	} else {
		gridname = gridname + std::string("Default") + std::string("_") + MeteoData::getParameterName(param) + ext;
	}
	
	return gridname;
}

double USERinterpolation::getQualityRating(const MeteoData::Parameters& param)
{
	const std::string filename = getGridFileName(param);

	if (!IOUtils::validFileName(filename)) {
		std::cout << "[E] Invalid grid filename for USER interpolation algorithm: " << filename << "\n";
		return 0.0;
	}
	if(IOUtils::fileExists(filename)) {
		return 1.0;
	} else {
		return 0.0;
	}
}


void USERinterpolation::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{
	const std::string filename = getGridFileName(param);
	//read2DGrid(grid, filename);
	throw IOException("USER interpolation algorithm not yet implemented...", AT);
}

} //namespace

