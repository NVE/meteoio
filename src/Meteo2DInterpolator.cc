#include "Meteo2DInterpolator.h"
#include "DEMObject.h"

using namespace std;
using namespace IOUtils;

Meteo2DInterpolator::Meteo2DInterpolator(const DEMObject& dem_in, const vector<MeteoData>& vecData,
					 const vector<StationData>& vecMeta) : dem(dem_in), SourcesData(vecData), SourcesMeta(vecMeta)
{
	//check whether the size of the two vectors is equal
	if (vecData.size() != vecMeta.size()) {
		throw IOException("Size of vector<MeteoData> and vector<StationData> are no equal", AT);
	}

}	

//This function calls the interpolation class for each individual meteo parameter.
//It also builds a list of valid data sources for the given parameter.
void Meteo2DInterpolator::interpolate(CArray2D<double>& nswc, CArray2D<double>& rh, CArray2D<double>& ta, 
				      CArray2D<double>& vw, CArray2D<double>& p)
{

	interpolateP(p);
	interpolateNSWC(nswc);
	interpolateTA(ta);
	interpolateRH(rh,ta);
	interpolateVW(vw);
}
void Meteo2DInterpolator::interpolate(CArray2D<double>& nswc, CArray2D<double>& rh, CArray2D<double>& ta,
					CArray2D<double>& vw, CArray2D<double>& p, CArray2D<double>& iswr/*, CArray2D<double>& ea*/)
{

	interpolateP(p);
	interpolateNSWC(nswc);
	interpolateTA(ta);
	interpolateRH(rh,ta);
	interpolateVW(vw);
	interpolateISWR(iswr);
	//interpolateEA(ea);
}

void Meteo2DInterpolator::interpolateNSWC(CArray2D<double>& nswc)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].nswc != nodata) {
			//cout << SourcesData[ii].nswc << endl;
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].nswc);
		}
	}

	printf("[i] interpolating NSWC using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D NSWC(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	NSWC.calculate(nswc);
}

void Meteo2DInterpolator::interpolateRH(CArray2D<double>& rh, CArray2D<double>& ta)
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

void Meteo2DInterpolator::interpolateTA(CArray2D<double>& ta)
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

	/*
	int nx, ny;
	ta.GetSize(nx,ny);
	std::cerr << "SourcesSelect.size()==" << vecSelectedStations.size() << "\nvecInput.size()==" << vecInput.size()
	<< "\nta size: " << nx << " x " << ny << std::endl;
	*/
	printf("[i] interpolating TA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, vecInput, vecSelectedStations, dem);
	TA.calculate(ta);
}

void Meteo2DInterpolator::interpolateDW(CArray2D<double>& dw)
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

void Meteo2DInterpolator::interpolateVW(CArray2D<double>& vw)
{	//HACK this is a quick and dirty fix for the wind interpolation...
	//HACK we *really* need a better design for the interpolations...
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();
	unsigned int countDataDir = 0.;
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

	//countDataDir=0.; //HACK, to prevent using the enhanced method...
	printf("[i] interpolating VW using %d stations\n", (int)vecSelectedStations.size());
	// If direction doesn't exist, use the kriging
	if( countDataDir > 0.) {
		unsigned int nx, ny;
		vw.GetSize(nx, ny);
		CArray2D<double> dw;
		dw.Create(nx,ny);
		interpolateDW(dw);
		Interpol2D VW(Interpol2D::I_CST, Interpol2D::I_VW, vecInput, vecSelectedStations, dem);
		VW.calculate(vw, vecEmpty, dw);
	} else {
		Interpol2D VW(Interpol2D::I_CST, Interpol2D::I_LAPSE_IDWK, vecInput, vecSelectedStations, dem);
		VW.calculate(vw);
	}

}

void Meteo2DInterpolator::interpolateP(CArray2D<double>& p)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;

	printf("[i] interpolating P using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D P(Interpol2D::I_PRESS, Interpol2D::I_PRESS, vecInput, vecSelectedStations, dem);
	P.calculate(p);
}

void Meteo2DInterpolator::interpolateISWR(CArray2D<double>& iswr)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].iswr != nodata) {
			//cout << SourcesData[ii].iswr << endl;
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].iswr);
		}
	}

	printf("[i] interpolating ISWR using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D ISWR(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	ISWR.calculate(iswr);
}

/*void Meteo2DInterpolator::interpolateEA(CArray2D<double>& ea)
{
	vector<StationData> vecSelectedStations;
	vector<double> vecInput;
	unsigned int datacount = SourcesData.size();

	for (unsigned int ii=0; ii<datacount; ii++) {
		if(SourcesData[ii].ea != nodata) {
			//cout << SourcesData[ii].iswr << endl;
			vecSelectedStations.push_back(SourcesMeta[ii]);
			vecInput.push_back(SourcesData[ii].ea);
		}
	}

	printf("[i] interpolating EA using %d stations\n", (int)vecSelectedStations.size());
	Interpol2D EA(Interpol2D::I_CST, Interpol2D::I_IDWK, vecInput, vecSelectedStations, dem);
	EA.calculate(ea);
}*/

