/*
 * ch_slf_gin_jnative_MeteoIOJNIInterface.cc
 *
 *  Created on: 08.01.2010
 *      Author: perot
 */

#ifdef _METEOIO_JNI

#include <jni.h>
#include "jnative.h"
#include "IOInterface.h"
#include "Grid2DObject.h"
#include "DEMObject.h"
#include "ch_slf_gin_jnative_MeteoIOJNIInterface.h"


jdoubleArray jMakeError (JNIEnv *env, float errorCode){
	jdoubleArray out = env->NewDoubleArray(1);
	jboolean isCopyOut;
	jdouble* dest = env->GetDoubleArrayElements(out, &isCopyOut);
	dest[0] = errorCode;
	if (isCopyOut == JNI_TRUE)
		env->ReleaseDoubleArrayElements(out, dest, 0);
	return out;
}

jdoubleArray convert_JNIArray(JNIEnv *env, const Grid2DObject&  p, const std::string& cellOrder){

    jboolean isCopyOut;
	jdoubleArray out = env->NewDoubleArray(p.nrows*p.ncols + 6);//size = dem.nrow*dem.ncols+6
	jdouble* dest = env->GetDoubleArrayElements(out, &isCopyOut);

	fulfillDoubleArray( p, cellOrder, dest);

	if (isCopyOut == JNI_TRUE)
		env->ReleaseDoubleArrayElements(out, dest, 0);

	return out;
}

JNIEXPORT jdoubleArray JNICALL Java_ch_slf_gin_jnative_MeteoIOJNIInterface_executeInterpolationSubDem
			  (JNIEnv *env, jclass theClass,
			  jstring jAlgorithm, jstring jIOinterface,
			  jstring jDemFile,jstring jDemCoordSystem,
			  jdouble demXll, jdouble demYll,
			  jdouble demXrt, jdouble demYrt,
			  jdoubleArray jMetadata, jdoubleArray jData,
			  jstring jMetaCoordSystem, jstring jcellOrder){

	const char * cDemFile = env->GetStringUTFChars(jDemFile,0);
	const char * cDemCoordSystem = env->GetStringUTFChars(jDemCoordSystem,0);
	const char * cIOInterface = env->GetStringUTFChars(jIOinterface,0);
	const char * cAlgorithm = env->GetStringUTFChars(jAlgorithm,0);
	const char * cMetaCoordSystem = env->GetStringUTFChars(jMetaCoordSystem,0);
	const char * cCellOrder = env->GetStringUTFChars(jcellOrder,0);

	IOInterface* io = getIOInterface(cDemFile, cDemCoordSystem, cIOInterface);
	if(io==NULL)
		return jMakeError(env, -1.f);

	//reading initial dem
	DEMObject dem = (demXll > -1 && demYrt> -1)?
			loadSubDEM(io,cDemCoordSystem,  demXll, demYll, demXrt, demYrt) :
			loadFullDEM(io);
	if (dem.nrows<1 || dem.ncols<2  ){
		std::cout << "Problem with DEM creation : "  << std::endl;
		//error
		return jMakeError(env, -2.f);
	}


	//Create MeteoData and StationData vectors
	int nbStation = (int)env->GetArrayLength(jMetadata)/3;
	int nbDataPerStation = (int)env->GetArrayLength(jData)/nbStation;
    jboolean isCopyMetadata;
    jboolean isCopyData;
	double *cMetadata = env->GetDoubleArrayElements(jMetadata,&isCopyMetadata);
	double *cData = env->GetDoubleArrayElements(jData,&isCopyData);
	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;
	//initialize MeteoData and StationData vectors
	loadMeteoAndStationData(cMetadata, cData, nbStation, nbDataPerStation, cAlgorithm,
			cMetaCoordSystem, &vecStation, &vecMeteo);


	Grid2DObject  p(dem.ncols, dem.nrows,
			dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	processInterpolation(cAlgorithm, p, dem, &vecStation, &vecMeteo);
	//copy the interpolation result into a jdoubleArray
	jdoubleArray out = convert_JNIArray(env, p, cCellOrder);

	//release cMetadata
	if (isCopyMetadata == JNI_TRUE)
		env->ReleaseDoubleArrayElements(jMetadata, cMetadata, JNI_ABORT);
	//release cData
	if (isCopyData == JNI_TRUE)
		env->ReleaseDoubleArrayElements(jData, cData, JNI_ABORT);
	delete io;

    env->ReleaseStringUTFChars(jDemFile, cDemFile);
    env->ReleaseStringUTFChars(jDemCoordSystem, cDemCoordSystem);
    env->ReleaseStringUTFChars(jIOinterface, cIOInterface);
    env->ReleaseStringUTFChars(jAlgorithm, cAlgorithm);
    env->ReleaseStringUTFChars(jMetaCoordSystem, cMetaCoordSystem);
    env->ReleaseStringUTFChars(jcellOrder, cCellOrder);
	vecMeteo.clear();
	vecStation.clear();
	return out;
}

JNIEXPORT jdoubleArray JNICALL Java_ch_slf_gin_jnative_MeteoIOJNIInterface_executeInterpolation
			(JNIEnv *env, jclass theClass,
			jstring jAlgorithm, jstring jIOinterface, jstring jDemFile, jstring jDemCoordSystem,
			jdoubleArray jMetadata, jdoubleArray jData,
			jstring jMetaCoordSystem, jstring jcellOrder){


	return Java_ch_slf_gin_jnative_MeteoIOJNIInterface_executeInterpolationSubDem(
			env, theClass,jAlgorithm, jIOinterface, jDemFile, jDemCoordSystem,
			-1, -1,-1, -1,jMetadata, jData, jMetaCoordSystem, jcellOrder);
}






#endif
