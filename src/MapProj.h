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

typedef void(MapProj::*convfunc)(double, double, double&, double&) const;

#ifdef _POPC_
class MapProj : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class MapProj {
#endif
 public:
	///Keywords for selecting the algorithm for computing geodesic distances
	enum GEO_DISTANCES {
		GEO_COSINE, ///< Spherical law of cosine
		GEO_VINCENTY ///< Vincenty ellispoid formula
	};

	//Constructors
	MapProj();
	MapProj(const std::string& coordinatesystem, const std::string& parameters="");
	MapProj(const double& _lat_ref, const double& _long_ref);

	//Operators
	bool operator==(const MapProj&) const; ///<Operator that tests for equality
	bool operator!=(const MapProj&) const; ///<Operator that tests for inequality

	//Coordinates conversions
	void convert_to_WGS84(double easting, double northing, double& latitude, double& longitude) const;
	void convert_from_WGS84(double latitude, double longitude, double& easting, double& northing) const;

	void WGS84_to_CH1903(double lat_in, double long_in, double& east_out, double& north_out) const;
	void CH1903_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_UTM(double lat_in, double long_in, double& east_out, double& north_out) const;
	void UTM_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;
	void WGS84_to_PROJ4(double lat_in, double long_in, double& east_out, double& north_out) const;
	void PROJ4_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;

	void WGS84_to_local(double lat_in, double long_in, double& east_out, double& north_out) const;
	void local_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;

	//Static methods
	static void WGS84_to_local(double lat_ref, double lon_ref, const double& lat, const double& lon, double& easting, double& northing, const enum GEO_DISTANCES algo=GEO_VINCENTY);
	static void local_to_WGS84(double lat_ref, double lon_ref, const double& easting, const double& northing, double& lat, double& lon, const enum GEO_DISTANCES algo=GEO_VINCENTY);

	static double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2);
	static void cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

	static double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	static double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2);
	static void VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

	static double dms_to_decimal(const std::string& dms);
	static void parseLatLon(const std::string& coordinates, double& lat, double& lon);
	static std::string decimal_to_dms(const double& decimal);

 private:
	void initializeMaps();
	void setFunctionPointers();
	int getUTMZone(const double latitude, const double longitude, std::string& zone_out) const;

	std::map<std::string, convfunc> to_wgs84;
	std::map<std::string, convfunc> from_wgs84;
	std::string coordsystem;
	std::string coordparam;
	convfunc convToWGS84, convFromWGS84;
	double ref_latitude, ref_longitude;
	
	///Keywords for selecting an ellipsoid to use
	enum ELLIPSOIDS_NAMES {
		E_WGS84, ///<Globally useable WGS84 ellipsoid
		E_GRS80, ///<GRS80 ellispoid, equivalent to WGS84 but used by America and Australia
		E_AIRY, ///<Airy ellispoid, good fit for the UK
		E_INTL1924, ///<International 1924 ellispoid, good for most of Europe
		E_CLARKE1880, ///<Clarke 1880, good for Africa
		E_GRS67 ///<GRS67 ellispoid, good for South America
	};
	struct ELLIPSOID {
		double a;
		double b;
	};
	static const struct ELLIPSOID ellipsoids[];

};

#endif
