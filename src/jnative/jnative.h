/*
 * jnative.h
 *
 *  Created on: 08.01.2010
 *      Author: perot
 */


#ifndef _Included_JNative
#define _Included_JNative

#if  defined(_METEOIO_JNI) ||  defined(_METEOIO_JNA)

#include "plugins/ARCIO.h"
//#include "plugins/BoschungIO.h"
#include "plugins/GeotopIO.h"
//#include "plugins/GSNIO.h"
#include "plugins/GrassIO.h"
#include "IOInterface.h"
#include "ConfigReader.h"
#include "DEMObject.h"


void loadMeteoAndStationData(double* cMetadata, double* cData,
		const int nbStation,const int nbDataPerStation,
		const std::string algorithm,const std::string metaCoordinateSystem,
		std::vector<StationData>* vecStation,
		 std::vector<double>* vecData, std::vector<double>* vecExtraData);

void processInterpolation(const std::string algorithm, Grid2DObject&  p, const  DEMObject& dem,
		std::vector<StationData>* vecStation,
		std::vector<double>* vecData, std::vector<double>* vecExtraData);

void fulfillDoubleArray(const Grid2DObject&  p, const std::string& cellOrder,
		double* dest);


#endif // defined(_METEOIO_JNI) ||  defined(_METEOIO_JNA)


#ifdef _METEOIO_JNA

/**
 *
 * Originally, these methods are dedicated to be called from JAVA with JNA framework.
 *
 *
 */


double* executeInterpolationSubDem
  (char*, char*, char*,char*, double, double, double, double,double*, int, double*, int, char*, char*);


double* executeInterpolation
(char*, char*, char*, char*, double*, int, double*, int, char*, char*);


#endif //_METEOIO_JNA

#endif//_Included_JNative
