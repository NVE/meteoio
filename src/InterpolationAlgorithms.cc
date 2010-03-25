#include "InterpolationAlgorithms.h"

using namespace std;

std::set<std::string> AlgorithmFactory::setAlgorithms;
const bool AlgorithmFactory::__init = AlgorithmFactory::initStaticData();

bool AlgorithmFactory::initStaticData()
{
	setAlgorithms.insert("CST");
	setAlgorithms.insert("STD_PRESS");
	setAlgorithms.insert("CST_LAPSE");
	setAlgorithms.insert("IDW");
	setAlgorithms.insert("IDW_LAPSE");
	setAlgorithms.insert("RH");
	setAlgorithms.insert("WIND");
	
	return true;
}

InterpolationAlgorithm* AlgorithmFactory::getAlgorithm(const string& _algoname, 
											const Meteo2DInterpolator& _mi, 
											const DEMObject& _dem, 
											const std::vector<MeteoData>& _vecMeteo, 
											const std::vector<StationData>& _vecStation,
											const std::vector<string>& _vecArgs)
{
	string algoname(_algoname);
	IOUtils::toUpper(algoname);

	//Check whether algorithm theoretically exists
	if (setAlgorithms.find(algoname) == setAlgorithms.end())
		throw UnknownValueException("The interpolation algorithm '"+algoname+"' does not exist" , AT);

	if (algoname == "CST"){
		return new ConstAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
	} else if (algoname == "STD_PRESS"){
		return new StandardPressureAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
	} else if (algoname == "CST_LAPSE"){
		return new ConstLapseRateAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
	} else if (algoname == "IDW"){
		return new IDWAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
	} else if (algoname == "IDW_LAPSE"){
		return new IDWLapseAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
	} else if (algoname == "RH"){
		return new RHAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
	} else if (algoname == "WIND"){
		return new SimpleWindInterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs);
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

unsigned int InterpolationAlgorithm::getStationAltitudes(const vector<StationData>& vecMeta,
											  vector<double>& vecData) const
{
	for (unsigned int ii=0; ii<vecMeta.size(); ii++){
		const double& val = vecMeta[ii].position.getAltitude();
		if (val != IOUtils::nodata)
			vecData.push_back(val);
	}

	return vecData.size();
}


/**********************************************************************************/


double ConstAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	unsigned int counter = getData(param, vecData);

	if (counter == 0){
		return 0.0;
	} else if (counter == 1){
		return 0.8;
	} else if (counter > 1){
		return 0.2;
	}

	return 0.0;
}


bool ConstAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{
	//get a vector of all the relevant data
	vector<double> vecData;
	unsigned int counter = getData(param, vecData);

	//check how many data points there are (not including nodata)
	if (counter == 0)
		return false;

	//run algorithm
	Interpol2D::constantGrid2DFill(Interpol1D::arithmeticMean(vecData), dem, grid);

	return true;
}


double StandardPressureAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	if (param != MeteoData::P)
		return 0.0;

	//check incoming data for number of data points
	vector<double> vecData;
	unsigned int counter = getData(param, vecData);

	if (counter == 0)
		return 1.0;

	return 0.3;
}

bool StandardPressureAlgorithm::calculate(const MeteoData::Parameters&, Grid2DObject& grid)
{
	//run algorithm
	Interpol2D::stdPressureGrid2DFill(dem, grid);

	return true;
}


double ConstLapseRateAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	unsigned int counter = getData(param, vecData);

	if (counter == 0){
		return 0.0;
	} else if (counter == 1){
		return 0.7;
	} else if (counter>1){
		return 0.2;
	}

	return 0.2;
}

bool ConstLapseRateAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	vector<double> vecAltitudes, vecData;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = getData(param, vecData, vecMeta);
	getStationAltitudes(vecMeta, vecAltitudes);

	LapseRateProjectPtr funcptr = &Interpol2D::LinProject;	

	if ((vecAltitudes.size() == 0) || (nrOfMeasurments == 0))
		return false;

	double avgAltitudes = Interpol1D::arithmeticMean(vecAltitudes);
	double avgData = Interpol1D::arithmeticMean(vecData);

	if (vecData.size() == 1)
		funcptr = &Interpol2D::ConstProject;

	//Set regression coefficients
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);

	//Get the argument for the filter: the default temperature lapse rate, otherwise use predefined
	if (vecArgs.size() < 1){
		vecCoefficients[1] = Interpol2D::dflt_temperature_lapse_rate;
	} else if (vecArgs.size() == 1){
		IOUtils::convertString(vecCoefficients[1], vecArgs[0]);
	} else {
		throw InvalidArgumentException("Wrong number of arguments supplied for the CST_LAPSE algorithm", AT);
	}

	//run algorithm
	Interpol2D::constantLapseGrid2DFill(avgData, avgAltitudes, dem, vecCoefficients, funcptr, grid);

	return true;
}


double IDWAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	unsigned int counter = getData(param, vecData);

	if (counter == 0){
		return 0.0;
	} else if (counter == 1){
		return 0.3;
	} else if (counter > 1){
		return 0.5;
	}

	return 0.2;
}

bool IDWAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	vector<double> vecData;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = getData(param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return false;

	//run algorithm
	Interpol2D::IDWKrieging(dem, vecData, vecMeta, grid);

	return true;
}

double IDWLapseAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//check incoming data for number of data points
	vector<double> vecData;
	unsigned int counter = getData(param, vecData);

	if (counter == 0)
		return 0.0;

	return 0.7;
}

bool IDWLapseAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	vector<double> vecAltitudes, vecData;
	vector<StationData> vecMeta;
	unsigned int nrOfMeasurments = getData(param, vecData, vecMeta);
	getStationAltitudes(vecMeta, vecAltitudes);

	if (nrOfMeasurments == 0)
		return false;

	//Set regression coefficients
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);

	//run algorithm
	Interpol2D::LinRegression(vecAltitudes, vecData, vecCoefficients);
	Interpol2D::LapseIDWKrieging(dem, &Interpol2D::LinProject, vecCoefficients, vecData, vecMeta, grid);

	return true;
}

double RHAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//This algorithm is only valid for RH
	if (param != MeteoData::RH)
		return 0.0;

	//check incoming data for number of data points

	vector<double> vecDataTA, vecDataRH;
	vector<StationData> vecMeta;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].rh != IOUtils::nodata) && (vecMeteo[ii].ta != IOUtils::nodata)){
			vecDataTA.push_back(vecMeteo[ii].ta);
			vecDataRH.push_back(vecMeteo[ii].rh);
			vecMeta.push_back(vecStation[ii]);
		}
	}

	if (vecDataTA.size() == 0)
		return 0.0;

	return 0.8;
}

bool RHAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	//This algorithm is only valid for RH
	if (param != MeteoData::RH)
		return 0.0;

	vector<double> vecAltitudes, vecDataTA, vecDataRH;
	vector<StationData> vecMeta;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].rh != IOUtils::nodata) && (vecMeteo[ii].ta != IOUtils::nodata)){
			vecDataTA.push_back(vecMeteo[ii].ta);
			vecDataRH.push_back(vecMeteo[ii].rh);
			vecMeta.push_back(vecStation[ii]);
		}
	}

	getStationAltitudes(vecMeta, vecAltitudes);

	if (vecDataTA.size() == 0) //No matching data
		return false;

	Grid2DObject ta;
	mi.interpolate(MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator

	//here, RH->Td, interpolations, Td->RH
	std::vector<double> vecTdStations(vecDataTA.size(), 0.0); // init to 0.0

	//Compute dew point temperatures at stations
	for (unsigned int ii=0; ii<vecDataRH.size(); ii++){
		vecTdStations[ii] = Interpol2D::RhtoDewPoint(vecDataRH[ii], vecDataTA[ii], 1);
	}
			
	//Krieging on Td
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);

	//run algorithm
	Interpol2D::LinRegression(vecAltitudes, vecTdStations, vecCoefficients);
	Interpol2D::LapseIDWKrieging(dem, &Interpol2D::LinProject, vecCoefficients, vecTdStations, vecStation, grid);

	//Recompute Rh from the interpolated td
	for (unsigned int ii=0; ii<grid.ncols; ii++) {
		for (unsigned int jj=0; jj<grid.nrows; jj++) {
			grid.grid2D(ii,jj) = Interpol2D::DewPointtoRh(grid.grid2D(ii,jj), ta.grid2D(ii,jj), 1);
		}
	}

	return true;
}

double SimpleWindInterpolationAlgorithm::getQualityRating(const MeteoData::Parameters& param)
{
	//This algorithm is only valid for VW
	if (param != MeteoData::VW)
		return 0.0;

	vector<double> vecAltitudes, vecDataVW, vecDataDW;
	vector<StationData> vecMeta;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].vw != IOUtils::nodata) && (vecMeteo[ii].dw != IOUtils::nodata)){
			vecDataVW.push_back(vecMeteo[ii].vw);
			vecDataDW.push_back(vecMeteo[ii].dw);
			vecMeta.push_back(vecStation[ii]);
		}
	}
	
	if (vecDataVW.size() == 0)
		return 0.0;

	return 1.0;
}

bool SimpleWindInterpolationAlgorithm::calculate(const MeteoData::Parameters& param, Grid2DObject& grid)
{		
	//This algorithm is only valid for VW
	if (param != MeteoData::VW)
		return 0.0;

	vector<double> vecAltitudes, vecDataVW, vecDataDW;
	vector<StationData> vecMeta;

	for (unsigned int ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii].vw != IOUtils::nodata) && (vecMeteo[ii].dw != IOUtils::nodata)){
			vecDataVW.push_back(vecMeteo[ii].vw);
			vecDataDW.push_back(vecMeteo[ii].dw);
			vecMeta.push_back(vecStation[ii]);
		}
	}

	getStationAltitudes(vecMeta, vecAltitudes);

	if( vecDataDW.size() == 0) return false;

	Grid2DObject dw;
	mi.interpolate(MeteoData::DW, dw); //get DW interpolation from call back to Meteo2DInterpolator
	
	//Krieging
	std::vector<double> vecCoefficients;
	vecCoefficients.resize(4, 0.0);
	
	Interpol2D::LinRegression(vecAltitudes, vecDataVW, vecCoefficients);
	Interpol2D::LapseIDWKrieging(dem, &Interpol2D::LinProject, vecCoefficients, vecDataVW, vecMeta, grid);
	Interpol2D::SimpleDEMWindInterpolate(dem, grid, dw);

	return true;
}

