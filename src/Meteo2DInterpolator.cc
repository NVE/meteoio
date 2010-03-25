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

#include "Meteo2DInterpolator.h"

using namespace std;

Meteo2DInterpolator::Meteo2DInterpolator(const ConfigReader& _cfg, const DEMObject& _dem, 
								 const std::vector<MeteoData>& _vecMeteo, 
								 const std::vector<StationData>& _vecStation) 
	: cfg(_cfg), dem(_dem), vecMeteo(_vecMeteo), vecStation(_vecStation)
{
	/*
	 * By reading the ConfigReader object build up a list of user configured algorithms
	 * for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, hnw, ...)
	 * Concept of this constructor: loop over all MeteoData::Parameters and then look
	 * for configuration of interpolation algorithms within the ConfigReader object.
	 */
	cout << "In Meteo2DInterpolator constructor" << endl;
	for (unsigned int ii=0; ii < MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		std::vector<std::string> tmpAlgorithms;
		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name
		unsigned int nrOfAlgorithms = getAlgorithmsForParameter(parname, tmpAlgorithms);

		cout << "Looking for interpoltions algorithms valid for: " << parname << endl;

		for (unsigned int jj=0; jj<nrOfAlgorithms; jj++){
			cout << parname<<" "<<jj<<": " << tmpAlgorithms[jj] << endl;
		}

		if (nrOfAlgorithms > 0)
			mapAlgorithms[parname] = tmpAlgorithms;
	}

	//check whether the size of the two vectors is equal
	if (vecMeteo.size() != vecStation.size())
		throw IOException("Size of vector<MeteoData> and vector<StationData> are not equal", AT);

	//check that the stations are using the same projection as the dem
	/*for (unsigned int i=0; i<(unsigned int)vecStation.size(); i++) {
		if(!vecStation[i].position.isSameProj(dem.llcorner)) {
			throw IOException("Some stations are not using the same geographic projection as the DEM", AT);
		}
		}*/
}

void Meteo2DInterpolator::interpolate(const MeteoData::Parameters& meteoparam, Grid2DObject& result) const
{

	//Show algorithms to be used for this parameter
	map<string, vector<string> >::const_iterator it = mapAlgorithms.find(MeteoData::getParameterName(meteoparam));
	
	if (it != mapAlgorithms.end()){
		cout << "Algorithms to be used for parameter " << MeteoData::getParameterName(meteoparam) << ":" << endl;

		double maxQualityRating = 0.0;
		auto_ptr<InterpolationAlgorithm> bestalgorithm(NULL);
		vector<string> vecArgs;

		for (unsigned int ii=0; ii < it->second.size(); ii++){
			const string& algoname = it->second.at(ii);
			getArgumentsForAlgorithm(meteoparam, algoname, vecArgs);
			
			//Get the configured algorithm
			auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(algoname, *this, dem, 
																	    vecMeteo, vecStation, vecArgs)); 
			//Get the quality rating and compare to already quality ratings
			double rating = algorithm->getQualityRating(meteoparam);
			if ((rating != 0.0) && (rating > maxQualityRating)){
				bestalgorithm = algorithm; //remember this algorithm: ownership belongs to bestalgorithm
			}

			cout << "\t" << it->second.at(ii) << "  rating:" << rating << endl;
		}

		//finally execute the algorithm with the best quality rating or throw an exception
		if (bestalgorithm.get() == NULL) {
			//result.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);
			throw IOException("No interpolation algorithm with quality rating >0 found", AT);
		} 
		
		if (!bestalgorithm->calculate(meteoparam, result)){
			throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(meteoparam), AT);
		}
	} else {
		//Some default message, that interpolation for this parameter needs configuration
		throw IOException("You need to configure the interpolation algorithms for parameter " + 
					   MeteoData::getParameterName(meteoparam), AT);
		
	}

	//check that the output grid is using the same projection as the dem
	if(!result.llcorner.isSameProj(dem.llcorner)) {
		throw IOException("The output grid is not using the same geographic projection as the DEM", AT);
	}

	/*
	if (meteoparam == MeteoData::P){
		interpolateP(result);
	} else if (meteoparam == MeteoData::HNW){
		interpolateHNW(result);
	} else if (meteoparam == MeteoData::TA){
		interpolateTA(result);
	} else if (meteoparam == MeteoData::RH){
		Grid2DObject ta(result.ncols, result.nrows, result.cellsize, result.llcorner);
		interpolateTA(ta);
		interpolateRH(result, ta);
	} else if (meteoparam == MeteoData::VW){
		interpolateVW(result);
	} else if (meteoparam == MeteoData::ISWR){
		interpolateISWR(result);
	} else {
		throw IOException("No interpolation algorithm for parameter " + MeteoData::getParameterName(meteoparam), AT);
		}*/
}

unsigned int Meteo2DInterpolator::getAlgorithmsForParameter(const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	/* 
	 * This function retrieves the user defined interpolation algorithms for 
	 * parameter 'parname' by querying the ConfigReader object
	 */
	std::vector<std::string> vecKeys;

	vecAlgorithms.clear();
	cfg.findKeys(vecKeys, parname+"::algorithms", "Interpolations2D");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::algorithms in config file", AT);;

	if (vecKeys.size() == 0)
		return 0;


	cfg.getValue(vecKeys.at(0), "Interpolations2D", vecAlgorithms, ConfigReader::nothrow);

	return vecAlgorithms.size();
}

unsigned int Meteo2DInterpolator::getArgumentsForAlgorithm(const MeteoData::Parameters& param, 
											    const string& algorithm, vector<string>& vecArgs) const
{
	vecArgs.clear();
	string keyname = MeteoData::getParameterName(param) +"::"+ algorithm;		
	cfg.getValue(keyname, "Interpolations2D", vecArgs, ConfigReader::nothrow);

	return vecArgs.size();
}

/***********************************************/
/*LEGACY                                       */
/***********************************************/
/*
void Meteo2DInterpolator::interpolateHNW(Grid2DObject& hnw) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecMeteo.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].hnw != IOUtils::nodata) {
			//cout << vecMeteo[ii].hnw << endl;
			vecSelectedStations.push_back(vecStation[ii]);
			vecInput.push_back(vecMeteo[ii].hnw);
		}
	}

	printf("[i] interpolating HNW using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D HNW(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	HNW.calculate(hnw);
}

void Meteo2DInterpolator::interpolateRH(Grid2DObject& rh, Grid2DObject& ta) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecExtraInput;
	std::vector<double> vecInput;
	const unsigned int datacount = vecMeteo.size();
	unsigned int rh_count=0;

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].rh != IOUtils::nodata) {
			//cout << vecMeteo[ii].rh << endl;
			rh_count++;
			if(vecMeteo[ii].ta != IOUtils::nodata) {
				vecSelectedStations.push_back(vecStation[ii]);
				vecInput.push_back(vecMeteo[ii].rh);
				vecExtraInput.push_back(vecMeteo[ii].ta);
			}
		}
	}

	if( ((int)vecSelectedStations.size() > (int)(0.5*rh_count)) && ((int)vecSelectedStations.size() >= 2) ) {
		printf("[i] interpolating RH using %d stations\n", (int)vecSelectedStations.size());
		Interpol2D RH(Interpol2D::I_CST, Interpol2D::I_RH, vecInput, vecSelectedStations, dem);
		RH.calculate(rh, vecExtraInput, ta);
	} else { //we are loosing too many stations when trying to get both rh and ta, trying a different strategy
		printf("[W] not enough stations with both TA and RH for smart RH interpolation (only %d from %d), using simpler IDWK\n",	
			(int)vecSelectedStations.size(),rh_count);
		vecSelectedStations.clear();
		vecInput.clear();
		vecExtraInput.clear();
		for (unsigned int ii=0; ii<datacount; ii++) {
			if(vecMeteo[ii].rh != IOUtils::nodata) {
				//cout << vecMeteo[ii].rh << endl;
				vecSelectedStations.push_back(vecStation[ii]);
				vecInput.push_back(vecMeteo[ii].rh);
			}
		}
		printf("[i] interpolating RH using %d stations\n", (int)vecSelectedStations.size());
		Interpol2D RH(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
		RH.calculate(rh);
	}

}

void Meteo2DInterpolator::interpolateTA(Grid2DObject& ta) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecMeteo.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].ta != IOUtils::nodata) {
			//cout << vecMeteo[ii].ta << endl;
			vecSelectedStations.push_back(vecStation[ii]);
			vecInput.push_back(vecMeteo[ii].ta);
		}
	}

	printf("[i] interpolating TA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, vecInput, vecSelectedStations, dem);
	TA.calculate(ta);
}

void Meteo2DInterpolator::interpolateDW(Grid2DObject& dw) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecMeteo.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].dw != IOUtils::nodata) {
			vecSelectedStations.push_back(vecStation[ii]);
			vecInput.push_back(vecMeteo[ii].dw);
		}
	}

	printf("[i] interpolating DW using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D DW(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	DW.calculate(dw);
}

void Meteo2DInterpolator::interpolateVW(Grid2DObject& vw) const
{	//HACK this is a quick and dirty fix for the wind interpolation...
	//HACK we *really* need a better design for the interpolations...
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecMeteo.size();
	unsigned int countDataDir = 0;
	std::vector<double> vecEmpty;

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].vw != IOUtils::nodata) {
			vecSelectedStations.push_back(vecStation[ii]);
			vecInput.push_back(vecMeteo[ii].vw);
		}
		if(vecMeteo[ii].dw != IOUtils::nodata) {
			countDataDir++;
		}
	}

	countDataDir=0; //HACK, to prevent using the enhanced method...
	printf("[i] interpolating VW using %d stations\n", (int)vecSelectedStations.size());
	// If direction doesn't exist, use the kriging
	if( countDataDir > 0.) {
		Grid2DObject dw(vw, 0, 0, vw.ncols, vw.nrows);
		interpolateDW(dw);
		Interpol2D VW(Interpol2D::I_CST, Interpol2D::I_VW, vecInput, vecSelectedStations, dem);
		VW.calculate(vw, vecEmpty, dw);
	} else {
		Interpol2D VW(Interpol2D::I_CST, Interpol2D::I_LAPSE_IDWK, vecInput, vecSelectedStations, dem);
		VW.calculate(vw);
	}

}

void Meteo2DInterpolator::interpolateP(Grid2DObject& p) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;

	printf("[i] interpolating P using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D P(Interpol2D::I_PRESS, Interpol2D::I_PRESS, vecInput, vecSelectedStations, dem);
	P.calculate(p);
}

void Meteo2DInterpolator::interpolateISWR(Grid2DObject& iswr) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecMeteo.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].iswr != IOUtils::nodata) {
			vecSelectedStations.push_back(vecStation[ii]);
			vecInput.push_back(vecMeteo[ii].iswr);
		}
	}

	printf("[i] interpolating ISWR using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D ISWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	ISWR.calculate(iswr);
}

void Meteo2DInterpolator::interpolateLWR(Grid2DObject& lwr) const
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecMeteo.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecMeteo[ii].lwr != IOUtils::nodata) {
			vecSelectedStations.push_back(vecStation[ii]);
			vecInput.push_back(vecMeteo[ii].lwr);
		}
	}

	printf("[i] interpolating EA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D LWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	LWR.calculate(lwr);
}
*/

#ifdef _POPC_
#include "marshal_meteoio.h"
void Meteo2DInterpolator::Serialize(POPBuffer &buf, bool pack)
{
//TODO: check this serialization!! It seems dubious that it would work at all...
	if (pack)
	{
		//buf.Pack(&cfg,1);
		/*buf.Pack(&dem,1);
		marshal_METEO_DATASET(buf, vecMeteo, 0, FLAG_MARSHAL, NULL);
		marshal_STATION_DATASET(buf, vecStation, 0, FLAG_MARSHAL, NULL);*/
		marshal_map_str_vecstr(buf, mapAlgorithms, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		//buf.UnPack(&cfg,1);
		/*buf.UnPack(&dem,1);
		marshal_METEO_DATASET(buf, vecMeteo, 0, !FLAG_MARSHAL, NULL);
		marshal_STATION_DATASET(buf, vecStation, 0, !FLAG_MARSHAL, NULL);*/
		marshal_map_str_vecstr(buf, mapAlgorithms, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif
