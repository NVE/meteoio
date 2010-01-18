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

Meteo2DInterpolator::Meteo2DInterpolator(const ConfigReader& _cfg, const DEMObject& _dem, 
								 const std::vector<MeteoData>& _vecData, 
								 const std::vector<StationData>& _vecMeta) 
	: cfg(_cfg), dem(_dem), vecData(_vecData), vecMeta(_vecMeta)
{
	/*
	 * By reading the ConfigReader object build up a list of user configured algorithms
	 * for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, hnw, ...)
	 * Concept of this constructor: loop over all MeteoData::Parameters and then look
	 * for configuration of interpolation algorithms within the ConfigReader object.
	 */
	for (unsigned int ii=0; ii < MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		std::vector<std::string> tmpAlgorithms;
		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name
		unsigned int nrOfAlgorithms = getAlgorithmsForParameter(parname, tmpAlgorithms);

		for (unsigned int jj=0; jj<nrOfAlgorithms; jj++){
			cout << parname<<" "<<jj<<": " << tmpAlgorithms[jj] << endl;
		}

		if (nrOfAlgorithms > 0)
			mapAlgorithms[parname] = tmpAlgorithms;
	}

	//check whether the size of the two vectors is equal
	if (vecData.size() != vecMeta.size())
		throw IOException("Size of vector<MeteoData> and vector<StationData> are no equal", AT);
}


void Meteo2DInterpolator::interpolate(const MeteoData::Parameters& meteoparam, Grid2DObject& result)
{
	if (meteoparam == MeteoData::P){
		interpolateP(result);
	} else if (meteoparam == MeteoData::HNW){
		interpolateHNW(result);
	} else if (meteoparam == MeteoData::TA){
		interpolateTA(result);
	} else if (meteoparam == MeteoData::RH){
		Grid2DObject ta;
		interpolateTA(ta);
		interpolateRH(result, ta);
	} else if (meteoparam == MeteoData::VW){
		interpolateVW(result);
	} else if (meteoparam == MeteoData::ISWR){
		interpolateISWR(result);
	} else {
		throw IOException("No interpolation algorithm for parameter " + MeteoData::getParameterName(meteoparam), AT);
	}
}

unsigned int Meteo2DInterpolator::getAlgorithmsForParameter(const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	/* 
	 * This function retrieves the user defined interpolation algorithms for 
	 * parameter 'parname' by querying the ConfigReader object
	 */
	std::vector<std::string> vecKeys;
	std::string tmp;

	vecAlgorithms.clear();
	cfg.findKeys(vecKeys, parname+"::algorithms", "Interpolations");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::algorithms in config file", AT);;

	if (vecKeys.size() == 0)
		return 0;

	cfg.getValue(vecKeys[0], "Interpolations", tmp, ConfigReader::nothrow);
	vecAlgorithms.push_back(tmp);

	return vecAlgorithms.size();
}

unsigned int Meteo2DInterpolator::getArgumentsForAlgorithm(const std::string& keyname, std::vector<std::string>& vecArguments)
{
	/*
	 * Retrieve the values for a given 'keyname' and store them in a vector calles 'vecArguments'
	 */
	cfg.getValue(keyname, "Interpolations", vecArguments, ConfigReader::nothrow);
	return vecArguments.size();
}


/***********************************************/
/*LEGACY                                       */
/***********************************************/
void Meteo2DInterpolator::interpolateHNW(Grid2DObject& hnw)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].hnw != nodata) {
			//cout << vecData[ii].hnw << endl;
			vecSelectedStations.push_back(vecMeta[ii]);
			vecInput.push_back(vecData[ii].hnw);
		}
	}

	printf("[i] interpolating HNW using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D HNW(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	HNW.calculate(hnw);
}

void Meteo2DInterpolator::interpolateRH(Grid2DObject& rh, Grid2DObject& ta)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecExtraInput;
	std::vector<double> vecInput;
	const unsigned int datacount = vecData.size();
	unsigned int rh_count=0;

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].rh != nodata) {
			//cout << vecData[ii].rh << endl;
			rh_count++;
			if(vecData[ii].ta != nodata) {
				vecSelectedStations.push_back(vecMeta[ii]);
				vecInput.push_back(vecData[ii].rh);
				vecExtraInput.push_back(vecData[ii].ta);
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
			if(vecData[ii].rh != nodata) {
				//cout << vecData[ii].rh << endl;
				vecSelectedStations.push_back(vecMeta[ii]);
				vecInput.push_back(vecData[ii].rh);
			}
		}
		printf("[i] interpolating RH using %d stations\n", (int)vecSelectedStations.size());
		Interpol2D RH(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
		RH.calculate(rh);
	}

}

void Meteo2DInterpolator::interpolateTA(Grid2DObject& ta)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].ta != nodata) {
			//cout << vecData[ii].ta << endl;
			vecSelectedStations.push_back(vecMeta[ii]);
			vecInput.push_back(vecData[ii].ta);
		}
	}

	printf("[i] interpolating TA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, vecInput, vecSelectedStations, dem);
	TA.calculate(ta);
}

void Meteo2DInterpolator::interpolateDW(Grid2DObject& dw)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].dw != nodata) {
			vecSelectedStations.push_back(vecMeta[ii]);
			vecInput.push_back(vecData[ii].dw);
		}
	}

	printf("[i] interpolating DW using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D DW(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	DW.calculate(dw);
}

void Meteo2DInterpolator::interpolateVW(Grid2DObject& vw)
{	//HACK this is a quick and dirty fix for the wind interpolation...
	//HACK we *really* need a better design for the interpolations...
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecData.size();
	unsigned int countDataDir = 0;
	std::vector<double> vecEmpty;

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].vw != nodata) {
			vecSelectedStations.push_back(vecMeta[ii]);
			vecInput.push_back(vecData[ii].vw);
		}
		if(vecData[ii].dw != nodata) {
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

void Meteo2DInterpolator::interpolateP(Grid2DObject& p)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;

	printf("[i] interpolating P using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D P(Interpol2D::I_PRESS, Interpol2D::I_PRESS, vecInput, vecSelectedStations, dem);
	P.calculate(p);
}

void Meteo2DInterpolator::interpolateISWR(Grid2DObject& iswr)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].iswr != nodata) {
			vecSelectedStations.push_back(vecMeta[ii]);
			vecInput.push_back(vecData[ii].iswr);
		}
	}

	printf("[i] interpolating ISWR using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D ISWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	ISWR.calculate(iswr);
}

void Meteo2DInterpolator::interpolateLWR(Grid2DObject& lwr)
{
	std::vector<StationData> vecSelectedStations;
	std::vector<double> vecInput;
	unsigned int datacount = vecData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(vecData[ii].lwr != nodata) {
			vecSelectedStations.push_back(vecMeta[ii]);
			vecInput.push_back(vecData[ii].lwr);
		}
	}

	printf("[i] interpolating EA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D LWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	LWR.calculate(lwr);
}

