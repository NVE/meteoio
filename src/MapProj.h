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
#ifndef __MAPPROJ_H__
#define __MAPPROJ_H__

#include "StationData.h"
#include "MeteoData.h"
#include "ConfigReader.h"
#include "IOExceptions.h"

#include <string>
#include <map>

class MapProj; //forward declaration

typedef void(MapProj::*convfunc)(const double&, const double&, double&, double&) const;

#ifdef _POPC_
class MapProj : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class MapProj {
#endif
 public:
	/**
	* @enum GEO_DISTANCES
	* keywords for selecting the algorithm for computing geodesic distances
	*/
	enum GEO_DISTANCES {
		GEO_COSINE, ///< Spehrical law of cosine
		GEO_VINCENTY ///< Vincenty ellispoid formula
	};
	
	/**
	* @brief Default constructor
	* This constructor builds a dummy object that performs no conversions but can be used for comparison
	* purpose. This is more or less the euqivalent of NULL for a pointer...
	*/
	MapProj();

	/**
	* @brief Regular constructor: usually, this is the constructor to use
	* @param[in] coordinatesystem string identifying the coordinate system to use
	* @param[in] parameters string giving some additional parameters for the projection (optional)
	*/
	MapProj(const std::string& coordinatesystem, const std::string& parameters="");

	/**
	* @brief Local projection onstructor: this constructor is only suitable for building a local projection.
	* Such a projection defines easting and northing as the distance (in meters) to a reference point
	* which coordinates have to be provided here.
	* @param[in] _lat_ref latitude of the reference point
	* @param[in] _long_ref longitude of the reference point
	*/
	MapProj(const double& _lat_ref, const double& _long_ref);

	/**
	* @brief Method converting towards WGS84
	* @param[in] easting easting of the coordinate to convert
	* @param[in] northing northing of the coordinate to convert
	* @param[out] latitude converted latitude
	* @param[out] longitude converted longitude
	*/
	void convert_to_WGS84(const double& easting, const double& northing, double& latitude, double& longitude) const;

	/**
	* @brief Method converting towards WGS84
	* @param[in] latitude latitude of the coordinate to convert
	* @param[in] longitude longitude of the coordinate to convert
	* @param[out] easting converted easting
	* @param[out] northing converted northing
	*/
	void convert_from_WGS84(const double& latitude, const double& longitude, double& easting, double& northing) const;

	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to Swiss grid
	* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
	* @param[in] lat_in Decimal Latitude 
	* @param[in] long_in Decimal Longitude
	* @param[out] east_out easting coordinate (Swiss system)
	* @param[out] north_out northing coordinate (Swiss system)
	*/
	void WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out) const;

	/**
	* @brief Coordinate conversion: from Swiss grid to WGS84 Lat/Long
	* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
	* @param[in] east_in easting coordinate (Swiss system)
	* @param[in] north_in northing coordinate (Swiss system)
	* @param[out] lat_out Decimal Latitude 
	* @param[out] long_out Decimal Longitude
	*/
	void CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out) const;

	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to proj4 parameters
	* @param lat_in Decimal Latitude 
	* @param long_in Decimal Longitude
	* @param east_out easting coordinate (target system)
	* @param north_out northing coordinate (target system)
	*/
	void WGS84_to_PROJ4(const double& lat_in, const double& long_in, double& east_out, double& north_out) const;

	/**
	* @brief Coordinate conversion: from proj4 parameters to WGS84 Lat/Long
	* @param east_in easting coordinate (Swiss system)
	* @param north_in northing coordinate (Swiss system)
	* @param lat_out Decimal Latitude 
	* @param long_out Decimal Longitude
	*/
	void PROJ4_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out) const;

	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to local grid as given in coordparam
	* @param lat_in Decimal Latitude 
	* @param long_in Decimal Longitude
	* @param east_out easting coordinate (target system)
	* @param north_out northing coordinate (target system)
	*/
	void WGS84_to_local(const double& lat_in, const double& long_in, double& east_out, double& north_out) const;

	/**
	* @brief Coordinate conversion: from ocal grid as given in coordparam to WGS84 Lat/Long
	* @param east_in easting coordinate (Swiss system)
	* @param north_in northing coordinate (Swiss system)
	* @param lat_out Decimal Latitude 
	* @param long_out Decimal Longitude
	*/
	void local_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out) const;

	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to local metric grid
	* @param lat_ref Decimal Latitude of origin (const double&)
	* @param lon_ref Decimal Longitude of origin (const double&)
	* @param lat Decimal Latitude (const double&)
	* @param lon Decimal Longitude(const double&)
	* @param easting easting coordinate in local grid (double&)
	* @param northing northing coordinate in local grid (double&)
	* @param algo select the algorith to be used for distances calculations (optional, default=GEO_VINCENTY)
	*/
	static void WGS84_to_local(const double& lat_ref, const double& lon_ref, const double& lat, const double& lon, double& easting, double& northing, const enum GEO_DISTANCES algo=GEO_VINCENTY);

	/**
	* @brief Coordinate conversion: from local metric grid to WGS84 Lat/Long
	* @param lat_ref Decimal Latitude of origin (const double&)
	* @param lon_ref Decimal Longitude of origin (const double&)
	* @param easting easting coordinate in local grid (double&)
	* @param northing northing coordinate in local grid (double&)
	* @param lat Decimal Latitude (double&)
	* @param lon Decimal Longitude(double&)
	* @param algo select the algorith to be used for distances calculations (optional, default=GEO_VINCENTY)
	*/
	static void local_to_WGS84(const double& lat_ref, const double& lon_ref, const double& easting, const double& northing, double& lat, double& lon, const enum GEO_DISTANCES algo=GEO_VINCENTY);

	/**
	* @brief Spherical law of cosine Distance calculation between points in WGS84 (decimal Lat/Long)
	* See http://www.movable-type.co.uk/scripts/latlong.html for more
	* @param[in] lat1 Decimal Latitude (const double&)
	* @param[in] lon1 Decimal Longitude (const double&)
	* @param[in] lat2 Decimal Latitude (const double&)
	* @param[in] lon2 Decimal Longitude (const double&)
	* @param[out] alpha average bearing
	* @return distance (double)
	*/
	static double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2);

	/**
	* @brief Spherical law of cosine Distance calculation between points in WGS84 (decimal Lat/Long)
	* See http://www.movable-type.co.uk/scripts/latlong.html for more
	* @param[in] lat_ref Decimal Latitude (const double&)
	* @param[in] lon_ref Decimal Longitude (const double&)
	* @param[in] distance Distance in meters (const double&)
	* @param[in] bearing bearing in degrees, 0 being north (const double&)
	* @param[out] lat Decimal latitude of target point (double&)
	* @param[out] lon Decimal longitude of target point (double&)
	*/
	static void cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

	/**
	* @brief Vincenty Distance calculation between points in WGS84 (decimal Lat/Long)
	* See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems", 
	* Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599, 
	* see http://www.springerlink.com/content/y7108u6862473583 for more
	* @param[in] lat1 Decimal Latitude (const double&)
	* @param[in] lon1 Decimal Longitude (const double&)
	* @param[in] lat2 Decimal Latitude (const double&)
	* @param[in] lon2 Decimal Longitude (const double&)
	* @param[in] alpha average bearing (double&)
	* @return distance (double)
	*/
	static double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2);

	/**
	* @brief Vincenty Inverse calculation giving WGS84 (decimal Lat/Long) position
	* given a start location (lat,lon) a distance and a bearing
	* See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems", 
	* Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599, 
	* see http://www.springerlink.com/content/y7108u6862473583 for more
	* @param[in] lat_ref Decimal Latitude (const double&)
	* @param[in] lon_ref Decimal Longitude (const double&)
	* @param[in] distance Distance in meters (const double&)
	* @param[in] bearing bearing in degrees, 0 being north (const double&)
	* @param[out] lat Decimal latitude of target point (double&)
	* @param[out] lon Decimal longitude of target point (double&)
	*/
	static void VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

	bool operator==(const MapProj&) const; ///<Operator that tests for equality
	bool operator!=(const MapProj&) const; ///<Operator that tests for inequality

 private:
	void initializeMaps();
	void setFunctionPointers();
	void parseLocalParameters(double& lat_ref, double& lon_ref) const;
	static double normalizeBearing(double angle);

	std::map<std::string, convfunc> to_wgs84;
	std::map<std::string, convfunc> from_wgs84;
	std::string coordsystem;
	std::string coordparam;
	convfunc convToWGS84, convFromWGS84;
	double ref_latitude, ref_longitude;
	
	enum ELLIPSOIDS_NAMES {
		E_WGS84,
		E_GRS80,
		E_AIRY,
		E_INTL1924,
		E_CLARKE1880,
		E_GRS67
	};
	struct ELLIPSOID {
		double a;
		double b;
	};
	static const struct ELLIPSOID ellipsoids[];

};

#endif
