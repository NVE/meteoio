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
#include <string>
#include <iostream>
#include <map>
#include <cmath>

namespace mio {

class Coords; //forward declaration
typedef void(Coords::*convfunc)(double, double, double&, double&) const;

/**
 * @class Coords
 * @brief A class to handle geographic coordinate systems.
 * This class offers an easy way to transparently convert between various coordinate systems. For any
 * given object, as soon as a latitude/longitude can be computed/read, it will be used as a reference.
 * This means that every subsequent change of projection system or read will be the conversion of this
 * reference lat/lon position (until a new "set" is called). See Coords::setProj for the supported coordinate systems.
 * @author Mathias Bavay
 * @date   2009-01-20
*/

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
		GEO_COSINE, ///< Spherical law of cosine (See http://www.movable-type.co.uk/scripts/latlong.html)
		GEO_VINCENTY ///< Vincenty ellispoid formula (See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems", Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599, or http://www.springerlink.com/content/y7108u6862473583 for more)
	} geo_distances;

	//Constructors
	Coords();
	Coords(const std::string& coordinatesystem, const std::string& parameters="");
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
	int getGridI() const;
	int getGridJ() const;
	int getGridK() const;
	void getProj(std::string& proj_type, std::string& proj_args) const;
	std::string printLatLon() const;

	friend std::ostream& operator<<(std::ostream& os, const Coords& coord);

	//Setter methods
	void setLatLon(const double _latitude, const double _longitude, const double _altitude, const bool _update=true);
	void setLatLon(const std::string& _coordinates, const double _altitude, const bool _update=true);
	void setXY(const double _easting, const double _northing, const double _altitude, const bool _update=true);
	void setGridIndex(const int _grid_i, const int _grid_j, const int _grid_k, const bool _invalidate=true);
	void setProj(const std::string& _coordinatesystem, const std::string& _parameters="");
	void setLocalRef(const double _ref_latitude, const double _ref_longitude);
	void setLocalRef(const std::string _coordparam);
	void setDistances(const geo_distances _algo);

	void check();
	double distance(const Coords& destination) const;
	bool isSameProj(const Coords& target) const;
	void copyProj(const Coords& source, const bool _update=true);

	//Static helper methods
	static double dms_to_decimal(const std::string& dms);
	static std::string decimal_to_dms(const double& decimal);
	static void parseLatLon(const std::string& coordinates, double& lat, double& lon);

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
	void clearCoordinates();
	void setDefaultValues();
	void initializeMaps();
	void setFunctionPointers();
	int getUTMZone(const double latitude, const double longitude, std::string& zone_out) const;

 private:
	double altitude; ///<altitude of the point (the altitude is currently NOT dependant on the projection)
	double latitude; ///<latitude of the point
	double longitude; ///<longitude of the point
	double easting; ///<east coordinate of the point in a cartesian grid
	double northing; ///<north coordinate of the point in a cartesian grid
	int grid_i; ///<grid index i (please notice that this index is NOT automatically regenerated NOR checked)
	int grid_j; ///<grid index j (please notice that this index is NOT automatically regenerated NOR checked)
	int grid_k; ///<grid index k (please notice that this index is NOT automatically regenerated NOR checked)
	double ref_latitude, ref_longitude;

	std::map<std::string, convfunc> to_wgs84;
	std::map<std::string, convfunc> from_wgs84;
	std::string coordsystem;
	std::string coordparam;
	convfunc convToWGS84, convFromWGS84;
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
} //end namespace

#endif
