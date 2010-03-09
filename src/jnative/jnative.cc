/*
 * jnative.cc
 *
 *  Created on: 08.01.2010
 *      Author: perot
 */


#if  defined(_METEOIO_JNI) ||  defined(_METEOIO_JNA)

#include "jnative.h"
#include "DEMLoader.h"
#include "DEMObject.h"
#include "libinterpol2D.h"

#include <time.h>


void loadMeteoAndStationData(double* cMetadata, double* cData,
		const int nbStation, const int nbDataPerStation,
		const std::string algorithm, const std::string metaCoordinateSystem,
		std::vector<StationData>* vecStation, std::vector<double>* vecData, std::vector<double>* vecExtraData){

	Date_IO date_in;
	Coords position(metaCoordinateSystem, "");

	for (int i = 0; i < nbStation; i++) {
		const double latitude = cMetadata[3*i];
		const double longitude = cMetadata[3*i+1];
		const double altitude = cMetadata[3*i+2];

		if (altitude<-5000 )
			continue;

		position.setXY(latitude, longitude, altitude);
		StationData station(position, "");
		vecStation->push_back(station);

		if(algorithm =="P" )
			vecData->push_back(cData[nbDataPerStation*i]);
		else if(algorithm == "HNW" )
			vecData->push_back(cData[nbDataPerStation*i]);
		else if(algorithm =="TA" )
			vecData->push_back(cData[nbDataPerStation*i]);
		else if(algorithm =="RH" ){
			vecData->push_back(cData[nbDataPerStation*i]);
			if (nbDataPerStation>1)
				vecExtraData->push_back(cData[nbDataPerStation*i+1]);
		}
		else if(algorithm =="VW" ){
			vecData->push_back(cData[nbDataPerStation*i]);
			if (nbDataPerStation>1)
				vecExtraData->push_back(cData[nbDataPerStation*i+1]);
		}
		else if(algorithm =="DW" )
			vecData->push_back(cData[nbDataPerStation*i]);
		else if(algorithm =="ISWR" )
			vecData->push_back(cData[nbDataPerStation*i]);
		else if(algorithm =="LWR" )
			vecData->push_back(cData[nbDataPerStation*i]);
		else//TA
			vecData->push_back(cData[nbDataPerStation*i]);
	}
}

void processInterpolation(const  std::string algorithm,
		Grid2DObject&  p, const  DEMObject& dem,
		std::vector<StationData>* vecStation,
		std::vector<double>* vecData,
		std::vector<double>* vecExtraData){

	std::cout << "processInterpolation "  << algorithm<< "*"  << vecStation->size() <<  std::endl;

	if(algorithm =="P"){
		Interpol2D P(Interpol2D::I_PRESS, Interpol2D::I_PRESS, *vecData, *vecStation, dem);
		P.calculate(p);
	}
	else if(algorithm =="HNW" ){
		Interpol2D HNW(Interpol2D::I_CST, Interpol2D::I_IDWK, *vecData, *vecStation, dem);
		HNW.calculate(p);
	}
	else if(algorithm =="TA" ){
		Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, *vecData, *vecStation, dem);
		TA.calculate(p);
	}
	else if(algorithm =="RH" ) {
		Grid2DObject ta(p, 0, 0, p.ncols, p.nrows);
		Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, *vecExtraData, *vecStation, dem);
		TA.calculate(ta);
		Interpol2D RH(Interpol2D::I_CST, Interpol2D::I_RH, *vecData, *vecStation, dem);
		RH.calculate(p, *vecExtraData, ta);
	}
	else if(algorithm =="VW" ){
		if( vecExtraData->size() > 0) { //extraData is DW
			std::vector<double> vecEmpty;
			Grid2DObject dw(p, 0, 0, p.ncols, p.nrows);
			Interpol2D DW(Interpol2D::I_CST, Interpol2D::I_IDWK, *vecExtraData, *vecStation, dem);
			DW.calculate(dw);
			Interpol2D VW(Interpol2D::I_CST, Interpol2D::I_VW, *vecData, *vecStation, dem);
			VW.calculate(p, vecEmpty, dw);
		} else {
			Interpol2D VW(Interpol2D::I_CST, Interpol2D::I_LAPSE_IDWK, *vecData, *vecStation, dem);
			VW.calculate(p);
		}
	}
	else if(algorithm =="DW" ){
		Interpol2D DW(Interpol2D::I_CST, Interpol2D::I_IDWK, *vecData, *vecStation, dem);
		DW.calculate(p);
	}
	else if(algorithm =="ISWR" ){
		Interpol2D ISWR(Interpol2D::I_CST, Interpol2D::I_IDWK, *vecData, *vecStation, dem);
		ISWR.calculate(p);
	}
	else if(algorithm =="LWR"){
		Interpol2D LWR(Interpol2D::I_CST, Interpol2D::I_IDWK, *vecData, *vecStation, dem);
		LWR.calculate(p);
	}
	else {
		Interpol2D TA(Interpol2D::I_LAPSE_CST, Interpol2D::I_LAPSE_IDWK, *vecData, *vecStation, dem);
		TA.calculate(p);
	}
}



void fulfillDoubleArray(const Grid2DObject&  p,
		const std::string& cellOrder,
		double* dest){

	dest[0] = 1.;//code for success
	dest[1] = p.ncols; //width
	dest[2] = p.nrows; //height
	dest[3] = 1.; //reserved ...
	dest[4] = 1.; //reserved ...
	dest[5] = 1.; //reserved ...

	if(cellOrder == "llur"){
		for (unsigned int kk = 0; kk < p.nrows; kk++)
			for (unsigned int ll=0; ll < p.ncols; ll++)
				dest[6+kk*p.ncols + ll] = p.grid2D(ll, kk);
			}
	else if(cellOrder == "urll" ){
		for (int kk = (signed)p.nrows-1; kk >=0; kk--)
			for (int ll=(signed)p.ncols -1; ll >=0; ll--)
				dest[6+(unsigned)kk*p.ncols + (unsigned)ll] = p.grid2D(ll, (unsigned)kk);
			}
	else if(cellOrder == "lrul" ){
		for (int kk = 0; kk < (signed)p.nrows; kk++)
			for (int ll=(signed)p.ncols -1; ll >=0; ll--)
				dest[6+(unsigned)kk*p.ncols + (unsigned)ll] = p.grid2D((unsigned)ll, (unsigned)kk);
			}
	else if(cellOrder == "ullr"){
		for (int kk = (signed)p.nrows-1; kk >=0; kk--)
			for (unsigned int ll=0; ll < p.ncols; ll++)
				dest[6+(unsigned)kk*p.ncols + ll] = p.grid2D(ll, (unsigned)kk);
			}
	else{
		for (unsigned int kk = 0; kk < p.nrows; kk++)
			for (unsigned int ll=0; ll < p.ncols; ll++)
				dest[6+kk*p.ncols + ll] = p.grid2D(ll, kk);
	}
}

#endif


#ifdef _METEOIO_JNA

double* makeError ( double errorCode){
	//PS. the array allocated here should be automatically deleted (and memory freed) when the
	//Java mapped object (com.sun.jna.ptr.DoubleByReference here) is "Garbage Collected"
	double* dest = (double*) malloc( sizeof(double));
	dest[0] = errorCode;
	return dest;
}

double* executeInterpolationSubDem(char* algorithm, char* iointerface,
		  char* demFile,char* demCoordSystem,
		  double demXll, double demYll, double demXrt, double demYrt,
		  double* metadata, int nbStation,
		  double* data, int nbDataPerStation, char* metaCoordSystem, char* cellOrder){

	clock_t tmpStart;
	clock_t tmpEnd;

	//get dem
	tmpStart = clock(); //start
	const DEMObject& dem = (demXll > -1 && demYrt> -1)?
			DEMLoader::loadSubDEM(demFile, demCoordSystem, iointerface, demXll, demYll, demXrt, demYrt) :
				DEMLoader::loadFullDEM(demFile, demCoordSystem, iointerface);
	if (dem.ncols<1 || dem.ncols<2  ){
		std::cout << "Problem with DEM creation : "  << std::endl;
		//error
		return makeError(-2.);
	}
	tmpEnd = clock(); //end
	double msDemLoading = (tmpEnd - tmpStart)/1000.0;


	//Create MeteoData and StationData vectors
	tmpStart = clock(); //start
	std::vector<double> vecData;
	std::vector<double> vecExtraData;
	std::vector<StationData> vecStation;
	//initialize MeteoData and StationData vectors
	loadMeteoAndStationData(metadata, data, nbStation, nbDataPerStation, algorithm,
			metaCoordSystem, &vecStation, &vecData, &vecExtraData);
	tmpEnd = clock(); //end
	double msDataLoading = (tmpEnd - tmpStart)/1000.0;


	//Interpolation
	tmpStart = clock(); //start
	Grid2DObject  p(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);
	bool success = true;
	try{
		processInterpolation(algorithm, p, dem, &vecStation, &vecData, &vecExtraData);
	}
	catch(IOException e){
		std::cout << "Interpolation failed : " << e.exception() << std::endl;
		success = false;
	}
	catch(...){
		std::cout << "Interpolation failed for some reason ?!? " <<  std::endl;
		success = false;
	}

	//PS. the array allocated here should be automatically deleted (and memory freed) when the
	//Java mapped object (com.sun.jna.ptr.DoubleByReference here) is "Garbage Collected"
	double* out = (success)?
			(double*) malloc( (6 + p.nrows*p.ncols)* sizeof(double)):
			(double*) malloc( 6 * sizeof(double));
	//copy the interpolation result into a double array
	if (success)
		fulfillDoubleArray( p, cellOrder, out);
	else{
		out[0] = -1.;//code for failing
		out[1] = 0;
		out[2] = 0;
	}
	tmpEnd = clock(); //end
	double msInterpolation = (tmpEnd - tmpStart)/1000.0;

	std::cout << " - time to load DEM : "  << msDemLoading << std::endl;
	std::cout << " - time to load Data : "  << msDataLoading << std::endl;
	std::cout << " - time to interpolate: "  << msInterpolation << std::endl;
	//put the different process in the result
	out[3] = msDemLoading;
	out[4] = msDataLoading;
	out[5] = msInterpolation;

	vecData.clear();
	vecExtraData.clear();
	vecStation.clear();
	return out;
}

double* executeInterpolation(char* algorithm, char* iointerface,
		  char* demFile,char* demCoordSystem,
		  double* metadata, int nbStation,
		  double* data, int nbDataPerStation,
		  char* metaCoordSystem, char* cellOrder){

	return executeInterpolationSubDem(algorithm,  iointerface, demFile,demCoordSystem,
			  -1.,-1.,-1.,-1.,
			 metadata, nbStation,
			 data, nbDataPerStation,
			 metaCoordSystem, cellOrder);
}

#endif











