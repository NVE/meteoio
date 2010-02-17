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
#ifndef __COORDS_H__
#define __COORDS_H__

#include "IOExceptions.h"
#include "IOUtils.h"
#include "ConfigReader.h"
#include <string>
#include <map>

class Coords; //forward declaration

typedef void(Coords::*convfunc)(double, double, double&, double&) const;

#ifdef _POPC_
class Coords : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Coords {
#endif
 public:
	///Keywords for selecting the algorithm for computing geodesic distances
	typedef enum GEO_DISTANCES {
		GEO_COSINE, ///< Spherical law of cosine
		GEO_VINCENTY ///< Vincenty ellispoid formula
	} geo_distances;

	//Constructors
	Coords();
	Coords(const std::string& coordinatesystem, const std::string& parameters="");
	Coords(const ConfigReader& cfg);
	Coords(const double& _lat_ref, const double& _long_ref);

	//Operators
	Coords& operator=(const Coords&); ///<Assignement operator
	bool operator==(const Coords&) const; ///<Operator that tests for equality
	bool operator!=(const Coords&) const; ///<Operator that tests for inequality

	//Getter methods
	double getEasting() const;
	double getNorthing() const;
	double getLat() const;
	double getLon() const;
	double getAltitude() const;
	std::string printLatLon() const;

	//Setter methods
	void setLatLon(const double _latitude, const double _longitude, const double _altitude=IOUtils::nodata, const bool _update=true);
	void setLatLon(const std::string& _coordinates, const double _altitude=IOUtils::nodata, const bool _update=true);
	void setXY(const double _easting, const double _northing, const double _altitude=IOUtils::nodata, const bool _update=true);
	void setProj(const std::string& _coordinatesystem, const std::string& _parameters="");
	void setLocalRef(const double _ref_latitude, const double _ref_longitude);
	void setLocalRef(const std::string _coordparam);
	void setDistances(const geo_distances _algo);

	void check();
	double distance(const Coords& destination) const;
	bool isSameProj(const Coords& target) const;
	void copyProj(const Coords& target, const bool _update=true);

	//Static helper methods
	static double dms_to_decimal(const std::string& dms);
	static std::string decimal_to_dms(const double& decimal);

 private:
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
	void WGS84_to_NULL(double lat_in, double long_in, double& east_out, double& north_out) const;
	void NULL_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const;

	//Distances calculations
	void distance(const Coords& destination, double& distance, double& bearing) const;
	double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha) const;
	void cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon) const;
	double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha) const;
	void VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon) const;

 private:
	void setDefaultValues();
	void initializeMaps();
	void setFunctionPointers();
	int getUTMZone(const double latitude, const double longitude, std::string& zone_out) const;
	void parseLatLon(const std::string& coordinates, double& lat, double& lon) const;

 private:
	double altitude;
	double latitude, longitude;
	double easting, northing;

	std::map<std::string, convfunc> to_wgs84;
	std::map<std::string, convfunc> from_wgs84;
	std::string coordsystem;
	std::string coordparam;
	convfunc convToWGS84, convFromWGS84;
	double ref_latitude, ref_longitude;
	geo_distances distance_algo;
	
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
