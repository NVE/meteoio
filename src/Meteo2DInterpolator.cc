#include "Meteo2DInterpolator.h"
#include "DEMObject.h"

using namespace std;
using namespace IOUtils;

Meteo2DInterpolator::Meteo2DInterpolator(const DEMObject& _dem, const vector<MeteoData>& vecData,
					 const vector<StationData>& vecMeta) : dem(_dem), SourcesData(vecData), SourcesMeta(vecMeta)
{
	//check whether the size of the two vectors is equal
	if (vecData.size() != vecMeta.size()) {
		throw IOException("Size of vector<MeteoData> and vector<StationData> are no equal", AT);
	}

}	

//This function calls the interpolation class for each individual meteo parameter.
//It also builds a list of valid data sources for the given parameter.
void Meteo2DInterpolator::interpolate(Grid2DObject& hnw, Grid2DObject& rh, Grid2DObject& ta, 
				      Grid2DObject& vw, Grid2DObject& p)
{

	interpolateP(p);
	interpolateHNW(hnw);
	interpolateTA(ta);
	interpolateRH(rh,ta);
	interpolateVW(vw);
}

void Meteo2DInterpolator::interpolate(Grid2DObject& hnw, Grid2DObject& rh, Grid2DObject& ta,
					Grid2DObject& vw, Grid2DObject& p, Grid2DObject& iswr/*, Grid2DObject& ea*/)
{

	interpolateP(p);
	interpolateHNW(hnw);
	interpolateTA(ta);
	interpolateRH(rh,ta);
	interpolateVW(vw);
	interpolateISWR(iswr);
	//interpolateEA(ea);
}

void Meteo2DInterpolator::interpolateHNW(Grid2DObject& hnw)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].hnw != nodata) {
			//cout << SourcesData[ii].hnw << endl;
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].hnw);
		}
	}

	printf("[i] interpolating HNW using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D HNW(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	HNW.calculate(hnw);
}

void Meteo2DInterpolator::interpolateRH(Grid2DObject& rh, Grid2DObject& ta)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecExtraInput;
	vector<double> vecInput;
	const unsigned int datacount = SourcesData.size();
	unsigned int rh_count=0;

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].rh != nodata) {
			//cout << SourcesData[ii].rh << endl;
			rh_count++;
			if(SourcesData[ii].ta != nodata) {
				vecSelectedStations.push_back(SourcesMeta[ii]);
				vecInput.push_back(SourcesData[ii].rh);
				vecExtraInput.push_back(SourcesData[ii].ta);
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
			if(SourcesData[ii].rh != nodata) {
				//cout << SourcesData[ii].rh << endl;
				vecSelectedStations.push_back(SourcesMeta[ii]);
				vecInput.push_back(SourcesData[ii].rh);
			}
		}
		printf("[i] interpolating RH using %d stations\n", (int)vecSelectedStations.size());
		Interpol2D RH(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
		RH.calculate(rh);
	}

}

void Meteo2DInterpolator::interpolateTA(Grid2DObject& ta)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].ta != nodata) {
			//cout << SourcesData[ii].ta << endl;
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].ta);
		}
	}

	printf("[i] interpolating TA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, vecInput, vecSelectedStations, dem);
	TA.calculate(ta);
}

void Meteo2DInterpolator::interpolateDW(Grid2DObject& dw)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].dw != nodata) {
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].dw);
		}
	}

	printf("[i] interpolating DW using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D DW(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	DW.calculate(dw);
}

void Meteo2DInterpolator::interpolateVW(Grid2DObject& vw)
{	//HACK this is a quick and dirty fix for the wind interpolation...
	//HACK we *really* need a better design for the interpolations...
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();
	unsigned int countDataDir = 0;
	vector<double> vecEmpty;

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].vw != nodata) {
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].vw);
		}
		if(SourcesData[ii].dw != nodata) {
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
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;

	printf("[i] interpolating P using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D P(Interpol2D::I_PRESS, Interpol2D::I_PRESS, vecInput, vecSelectedStations, dem);
	P.calculate(p);
}

void Meteo2DInterpolator::interpolateISWR(Grid2DObject& iswr)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].iswr != nodata) {
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].iswr);
		}
	}

	printf("[i] interpolating ISWR using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D ISWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	ISWR.calculate(iswr);
}

void Meteo2DInterpolator::interpolateLWR(Grid2DObject& lwr)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].lwr != nodata) {
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].lwr);
		}
	}

	printf("[i] interpolating EA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D LWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	LWR.calculate(lwr);
}

