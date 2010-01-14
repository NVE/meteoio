/*
 * jnative.cc
 *
 *  Created on: 08.01.2010
 *      Author: perot
 */


#if  defined(_METEOIO_JNI) ||  defined(_METEOIO_JNA)

#include "plugins/ARCIO.h"
//#include "plugins/BoschungIO.h"
#include "plugins/GeotopIO.h"
#include "plugins/GrassIO.h"
#include "IOInterface.h"
#include "ConfigReader.h"
#include "Meteo2DInterpolator.h"
#include "DEMObject.h"


IOInterface* getIOInterface(const std::string cDemFile,
		const std::string cDemCoordSystem, const std::string interfaceType){

	IOInterface *io = NULL;
	try {
		ConfigReader cfg;
		cfg.addKey("DEMFILE", cDemFile);
		cfg.addKey("COORDIN", cDemCoordSystem);
		cfg.addKey("COORDPARAM","");
		if(interfaceType == "ARCIO")
			io = new ARCIO(cfg);
		//else if(interfaceType ==  "BoschungIO" ): io = new BoschungIO(cfg);
		else if(interfaceType == "GeotopIO" )
			io = new GeotopIO(cfg);
		else if(interfaceType == "GrassIO" )
			io = new GrassIO(cfg);
		else
			io = new ARCIO(cfg); //default IOinterface
	}catch (IOException& e){
		std::cout << "Problem with ARCIO creation, cause: " << e.what() << std::endl;
		//error
		return NULL ;
	}
	return io;
}

DEMObject loadSubDEM(IOInterface* io, const std::string cDemCoordSystem,
		const double demXll, const double demYll, const double demXur, const double demYur){
	//reading initial dem
	DEMObject dem;
	io->readDEM(dem);
	if (dem.ncols<1 || dem.ncols<2  )
		return dem;

	//compute WGS coordinates (considered as the true reference)
	double latll, longll, latur, longur;
	MapProj mymapproj( cDemCoordSystem, "");
	mymapproj.convert_to_WGS84(demXll, demYll, latll, longll);
	mymapproj.convert_to_WGS84(demXur, demYur, latur, longur);
	//retrieving grid coordinates of a real world point
	unsigned int i0,j0,i1,j1;
	dem.WGS84_to_grid(latll, longll, i0,j0);
	dem.WGS84_to_grid(latur, longur, i1,j1);
	//extracting a sub-dem
	DEMObject sub_dem(dem, i0, j0, i1-i0+1, j0-j1+1);
	return sub_dem;
}

DEMObject loadFullDEM(IOInterface* io){
	//reading initial dem
	DEMObject dem;
	io->readDEM(dem);
	return dem;
}

void loadMeteoAndStationData(double* cMetadata, double* cData,
		const int nbStation, const int nbDataPerStation,
		const std::string algorithm, const std::string metaCoordinateSystem,
		std::vector<StationData>* vecStation, std::vector<MeteoData>* vecMeteo){

	Date_IO date_in;
	for (int i = 0; i < nbStation; i++){

		double xllcorner = cMetadata[3*i];
		double yllcorner = cMetadata[3*i+1];
		double altitude = cMetadata[3*i+2];
		if (altitude<-5000 )
			continue;
		double latitude, longitude;
		if (metaCoordinateSystem != "WGS84"){
			MapProj mymapproj(metaCoordinateSystem, "");
			mymapproj.convert_to_WGS84(xllcorner, yllcorner, latitude, longitude);
		}
		else{
			latitude = yllcorner;
			longitude = xllcorner;
		}
		StationData station(xllcorner, yllcorner,
				altitude, "",latitude, longitude);
		vecStation->push_back(station);

		if(algorithm =="P" ) {
			MeteoData meteoData(date_in,nodata,nodata,nodata,nodata,nodata,nodata,nodata,nodata,nodata,nodata,nodata,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm == "HNW" ){
			MeteoData meteoData(date_in,nodata,nodata,nodata,nodata,nodata,nodata,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm =="TA" ) {
			MeteoData meteoData(date_in,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm =="RH" ){ //Air temperature is defined before RH, so if both are delivered, TA comes first
			MeteoData meteoData(date_in,(nbDataPerStation>1)?cData[nbDataPerStation*i]:nodata,nodata,nodata,nodata,(nbDataPerStation>1)?cData[nbDataPerStation*i+1]:cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm =="VW" ){//DW optional ?
			MeteoData meteoData(date_in,nodata,nodata,cData[nbDataPerStation*i],(nbDataPerStation>1)?cData[nbDataPerStation*i+1]:nodata);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm =="DW" ){
			MeteoData meteoData(date_in,nodata,nodata,nodata,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm =="ISWR" ){
			MeteoData meteoData(date_in,nodata,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else if(algorithm =="LWR" ){
			MeteoData meteoData(date_in,nodata,nodata,nodata,nodata,nodata,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
		else{//TA
			MeteoData meteoData(date_in,cData[nbDataPerStation*i]);
			vecMeteo->push_back(meteoData);
		}
	}
}

void processInterpolation(const  std::string algorithm,Grid2DObject&  p, const  DEMObject& dem,
		std::vector<StationData>* vecStation, std::vector<MeteoData>* vecMeteo){

	std::cout << "processInterpolation "  << algorithm<< "*"  << vecStation->size() <<  std::endl;
	Meteo2DInterpolator mi(dem, *vecMeteo, *vecStation);
	if(algorithm =="P")
		mi.interpolateP(p);
	else if(algorithm =="HNW" )
		mi.interpolateHNW(p);
	else if(algorithm =="TA" )
		mi.interpolateTA(p);
	else if(algorithm =="RH" ) {
		Grid2DObject  ta(dem.ncols, dem.nrows,
				dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
		mi.interpolateRH(p,ta);
	}
	else if(algorithm =="VW" )
		mi.interpolateVW(p);
	else if(algorithm =="DW" )
		mi.interpolateDW(p);
	else if(algorithm =="ISWR" )
		mi.interpolateISWR(p);
	else if(algorithm =="LWR")
		mi.interpolateLWR(p);
	else 	mi.interpolateTA(p);
}



void fulfillDoubleArray(const Grid2DObject&  p, const std::string& cellOrder,
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
		for (unsigned int kk = p.nrows-1; kk >=0; kk--)
			for (unsigned int ll=p.ncols -1; ll >=0; ll--)
				dest[6+kk*p.ncols + ll] = p.grid2D(ll, kk);
			}
	else if(cellOrder == "lrul" ){
		for (unsigned int kk = 0; kk < p.nrows; kk++)
			for (unsigned int ll=p.ncols -1; ll >=0; ll--)
				dest[6+kk*p.ncols + ll] = p.grid2D(ll, kk);
			}
	else if(cellOrder == "ullr"){
		for (unsigned int kk = p.nrows-1; kk >=0; kk--)
			for (unsigned int ll=0; ll < p.ncols; ll++)
				dest[6+kk*p.ncols + ll] = p.grid2D(ll, kk);
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


	IOInterface* io = getIOInterface(demFile, demCoordSystem, iointerface);
	if(io==NULL)
		return makeError(-1.);

	//reading initial dem
	DEMObject dem = (demXll > -1 && demYrt> -1)?
			loadSubDEM(io,demCoordSystem,  demXll, demYll, demXrt, demYrt) :
			loadFullDEM(io);
	if (dem.ncols<1 || dem.ncols<2  ){
		std::cout << "Problem with DEM creation : "  << std::endl;
		//error
		return makeError(-2.);
	}


	//Create MeteoData and StationData vectors
	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;
	//initialize MeteoData and StationData vectors
	loadMeteoAndStationData(metadata, data, nbStation, nbDataPerStation, algorithm,
			metaCoordSystem, &vecStation, &vecMeteo);


	Grid2DObject  p(dem.ncols, dem.nrows,
			dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	processInterpolation(algorithm, p, dem, &vecStation, &vecMeteo);

	//PS. the array allocated here should be automatically deleted (and memory freed) when the
	//Java mapped object (com.sun.jna.ptr.DoubleByReference here) is "Garbage Collected"
	double* out = (double*) malloc( (6 + p.nrows*p.ncols)* sizeof(double));
	//copy the interpolation result into a double array
	fulfillDoubleArray( p, cellOrder, out);

	delete io;
	vecMeteo.clear();
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











