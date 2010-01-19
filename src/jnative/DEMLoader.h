/*
 * DEMLoader.h
 *
 *  Created on: 19.01.2010
 *      Author: perot
 */
#ifndef _Included_DEMLoader
#define _Included_DEMLoader


#include "IOInterface.h"
#include "Grid2DObject.h"
#include "DEMObject.h"
#include "Grid2DObject.h"
#include "DEMObject.h"


typedef std::map<std::string, DEMObject>  demMapType;
typedef std::pair<std::string, DEMObject>  demPairType;


/**
 * Class with an unique instance loading and storing DEMobjects.
 */
class DEMLoader
{
private:
  DEMLoader () { }
  ~DEMLoader () {
	  demMap.clear();
  }


  static DEMLoader *getInstance (){
    if (NULL == _singleton){
        _singleton =  new DEMLoader;
     }
    return _singleton;
  }

  const DEMObject& internal_loadFullDEM(const std::string cDemFile,
 			const std::string  cDemCoordSystem, const std::string cInterfaceType);

  const DEMObject& internal_loadSubDEM(const std::string cDemFile,
  			const std::string  cDemCoordSystem, const std::string cInterfaceType ,
  			const double xll , const double yll, const double xur, const double yur);

  /**
   * Returns an implementation of IOInterface with the appropriate
   * configuration. Each call generates a new Instance -> delete it after use !
   * TODO : Store the instances and reuse them, even in a concurent acceses context.
   *        Pb of the ConfigReader member which cannot be mutualized throught different interplation !
   *
   */
  IOInterface* generateIOInterface(const std::string cDemFile,
		  const std::string  cDemCoordSystem,
		  const std::string cInterfaceType);

public:

	static const DEMObject& loadFullDEM(const std::string cDemFile,
			const std::string  cDemCoordSystem, const std::string cInterfaceType){
		return DEMLoader::getInstance()->internal_loadFullDEM(cDemFile, cDemCoordSystem, cInterfaceType);
	}

	static const DEMObject& loadSubDEM(const std::string cDemFile,
			const std::string  cDemCoordSystem, const std::string cInterfaceType ,
			const double xll , const double yll, const double xur, const double yur){

		return DEMLoader::getInstance()->internal_loadSubDEM(cDemFile,
				cDemCoordSystem,cInterfaceType , xll , yll, xur,yur);
	}

	/**
	* Public static access to destroy (interest ?!?)
	*/
	static void kill (){
		if (NULL != _singleton) {
			delete _singleton; //call destructor
			_singleton = NULL;
		  }
	  }

private:

	/**
	 * Structure containing the DEMObjects. The key is the dem file name.
	 */
	demMapType  demMap;

	/**
	 * Static pointer to the unique instance, if null the instance had not been initialized yet
	 */
	static DEMLoader *_singleton;
};

#endif
