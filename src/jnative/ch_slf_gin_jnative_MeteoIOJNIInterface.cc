/*
 * ch_slf_gin_jnative_MeteoIOJNIInterface.cc
 *
 *  Created on: 08.01.2010
 *      Author: perot
 */

#ifdef _METEOIO_JNI

#include <jni.h>
#include "jnative.h"
#include "Grid2DObject.h"
#include "DEMObject.h"
#include "DEMLoader.h"
#include "ch_slf_gin_jnative_MeteoIOJNIInterface.h"

#include <stdlib.h>
#include <time.h>


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

	clock_t tmpStart;
	clock_t tmpEnd;

	tmpStart = clock(); //start
	//get dem
	const DEMObject& dem = (demXll > -1 && demYrt> -1)?
			DEMLoader::loadSubDEM(cDemFile, cDemCoordSystem, cIOInterface, demXll, demYll, demXrt, demYrt) :
				DEMLoader::loadFullDEM(cDemFile, cDemCoordSystem, cIOInterface);
	if (dem.nrows<1 || dem.ncols<2  ){
		std::cout << "Problem with DEM creation : "  << std::endl;
		//error
		return jMakeError(env, -2.f);
	}
	tmpEnd = clock(); //end
	double msDemLoading = (tmpEnd - tmpStart)/1000.0;


	//Create MeteoData and StationData vectors
	tmpStart = clock(); //start
	int nbStation = (int)env->GetArrayLength(jMetadata)/3;
	int nbDataPerStation = (int)env->GetArrayLength(jData)/nbStation;
    jboolean isCopyMetadata;
    jboolean isCopyData;
	double *cMetadata = env->GetDoubleArrayElements(jMetadata,&isCopyMetadata);
	double *cData = env->GetDoubleArrayElements(jData,&isCopyData);
	std::vector<double> vecData;
	std::vector<double> vecExtraData;
	std::vector<StationData> vecStation;
	//initialize MeteoData and StationData vectors
	loadMeteoAndStationData(cMetadata, cData, nbStation, nbDataPerStation, cAlgorithm,
			cMetaCoordSystem, &vecStation, &vecData, &vecExtraData);
	tmpEnd = clock(); //end
	double msDataLoading = (tmpEnd - tmpStart)/1000.0;

	//Interpolation
	tmpStart = clock(); //start
	Grid2DObject  p(dem.ncols, dem.nrows,
			dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	processInterpolation(cAlgorithm, p, dem, &vecStation, &vecData, &vecExtraData);
	//copy the interpolation result into a jdoubleArray
	jdoubleArray out = convert_JNIArray(env, p, cCellOrder);
	tmpEnd = clock(); //end
	double msInterpolation = (tmpEnd - tmpStart)/1000.0;

	//put the different process in the result
	double* times = (double*) malloc( 3* sizeof(double));
	times[0] = msDemLoading;
	times[1] = msDataLoading;
	times[2] = msInterpolation;
	env->SetDoubleArrayRegion(out, 3, 3, times);
	free(times);
	std::cout << " - time to load DEM : "  << msDemLoading << std::endl;
	std::cout << " - time to load Data : "  << msDataLoading << std::endl;
	std::cout << " - time to interpolate: "  << msInterpolation << std::endl;

	//release cMetadata
	if (isCopyMetadata == JNI_TRUE)
		env->ReleaseDoubleArrayElements(jMetadata, cMetadata, JNI_ABORT);
	//release cData
	if (isCopyData == JNI_TRUE)
		env->ReleaseDoubleArrayElements(jData, cData, JNI_ABORT);
    env->ReleaseStringUTFChars(jDemFile, cDemFile);
    env->ReleaseStringUTFChars(jDemCoordSystem, cDemCoordSystem);
    env->ReleaseStringUTFChars(jIOinterface, cIOInterface);
    env->ReleaseStringUTFChars(jAlgorithm, cAlgorithm);
    env->ReleaseStringUTFChars(jMetaCoordSystem, cMetaCoordSystem);
    env->ReleaseStringUTFChars(jcellOrder, cCellOrder);
    vecData.clear();
	vecExtraData.clear();
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
