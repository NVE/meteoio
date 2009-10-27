#ifndef __MAPPROJ_H__
#define __MAPPROJ_H__

#include "StationData.h"
#include "MeteoData.h"
#include "ConfigReader.h"
#include "IOExceptions.h"

#include <string>
#include <map>

class MapProj; //forward declaration

typedef void(MapProj::*convfunc)(const double&, const double&, double&, double&);

class MapProj {
 public:
	MapProj(const std::string& coordinatesystem, const std::string& parameters="");

	void convert_to_WGS84(const double& easting, const double& northing, double& latitude, double& longitude);
	void convert_from_WGS84(const double& latitude, const double& longitude, double& easting, double& northing);

	
	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to Swiss grid
	* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
	* @param lat_in Decimal Latitude 
	* @param long_in Decimal Longitude
	* @param east_out easting coordinate (Swiss system)
	* @param north_out northing coordinate (Swiss system)
	*/
	void WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out);

	/**
	* @brief Coordinate conversion: from Swiss grid to WGS84 Lat/Long
	* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
	* @param east_in easting coordinate (Swiss system)
	* @param north_in northing coordinate (Swiss system)
	* @param lat_out Decimal Latitude 
	* @param long_out Decimal Longitude
	*/
	void CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out);

	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to proj4 parameters
	* @param lat_in Decimal Latitude 
	* @param long_in Decimal Longitude
	* @param east_out easting coordinate (target system)
	* @param north_out northing coordinate (target system)
	*/
	void WGS84_to_PROJ4(const double& lat_in, const double& long_in, double& east_out, double& north_out);

	/**
	* @brief Coordinate conversion: from proj4 parameters to WGS84 Lat/Long
	* @param east_in easting coordinate (Swiss system)
	* @param north_in northing coordinate (Swiss system)
	* @param lat_out Decimal Latitude 
	* @param long_out Decimal Longitude
	*/
	void PROJ4_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out);


 private:
	void initializeMaps();

	std::map<std::string, convfunc> to_wgs84;
	std::map<std::string, convfunc> from_wgs84;
	std::string coordsystem;
	std::string coordparam;
	convfunc convToWGS84, convFromWGS84;
};

#endif
