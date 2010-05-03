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
#include "Coords.h"

#ifdef PROJ4
	#include <proj_api.h>
#endif

using namespace std;
using namespace mio;

 /**
 * @page coords Available coordinate systems
 * Geographic coordinates will be transparently and automatically converted to lat/lon and any other coordinate system that
 * the client program uses. However, in order to do so, the input coordinate system must be specified. In order to output
 * geolocalized data, the desired coordinate system must also be specified for the outputs (in the output section).
 * This is done through the use of the COORDIN and COORDPARAM keys (see the documentation for each plugin).
 *
 * There are two ways of supporting a given coordinate system: through the use of an adhoc implementation
 * (that becomes part of MeteoIO) or through the use of an external library, Proj4 [ref: http://trac.osgeo.org/proj/].
 * The current internal implementations are the following (given by their keyword):
 * - CH1903 for coordinates in the Swiss Grid [ref: http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf]
 * - UTM for UTM coordinates (the zone must be specified in the parameters, for example 31T) [ref: http://www.oc.nps.edu/oc2902w/maps/utmups.pdf]
 * - LOCAL for local coordinate system (using the horizontal and vertical distance from a reference point, see Coords::geo_distances for the available choice of distance algorithms)
 *
 * Such an example of use is the following:
 * @code
 * COORDSYS	= UTM
 * COORDPARAM	= 31T
 * @endcode
 *
 * On the other hand, when using the Proj4 library for handling the coordinate conversion, the proj4 conversion
 * string must be specified in the parameters. For example (for UTM zone 10 coordinates):
 * @code
 * COORDSYS	= PROJ4
 * COORDPARAM	= +proj=utm +ellps=WGS84 +zone=10
 * @endcode
 * It is also possible to use EPSG codes (such codes can be found at http://spatialreference.org/ref/epsg/?page=1)
 * as illustrated below (21781 is the EPSG code for the CH1903 coordinate system):
 * @code
 * COORDSYS	= PROJ4
 * COORDPARAM	= +init=epsg:21781
 * @endcode
 *
 */

const struct Coords::ELLIPSOID Coords::ellipsoids[] = {
	{ 6378137.,	6356752.3142 }, ///< E_WGS84
	{ 6378137.,	6356752.3141 }, ///< E_GRS80
	{ 6377563.396,	6356256.909 }, ///< E_AIRY
	{ 6378388.,	6356911.946 }, ///< E_INTL1924
	{ 6378249.145,	6356514.86955 }, ///< E_CLARKE1880
	{ 6378160.,	6356774.719 } ///< E_GRS67
};

/**
* @brief Equality operator that checks that lat/lon match
* @param[in] in Coord object to compare to
* @return true or false
*/
bool Coords::operator==(const Coords& in) const {
//check that two Coords objects represent the same location

	if(latitude==IOUtils::nodata || longitude==IOUtils::nodata || easting==IOUtils::nodata || northing==IOUtils::nodata) {
		//only available information is grid indices
		const bool comparison = ( grid_i==in.grid_i && grid_j==in.grid_j && grid_k==in.grid_k );
		return comparison;
	} else {
		const bool comparison = ( IOUtils::checkEpsilonEquality(getLat(), in.getLat(), IOUtils::lat_epsilon) &&
				IOUtils::checkEpsilonEquality(getLon(), in.getLon(), IOUtils::lon_epsilon) );
		return comparison;
	}
}

/**
* @brief Inequality operator that checks that lat/lon don't match
* @param[in] in Coord object to compare to
* @return true or false
*/
bool Coords::operator!=(const Coords& in) const {
	return !(*this==in);
}

Coords& Coords::operator=(const Coords& source) {
	if(this != &source) {
		altitude = source.altitude;
		latitude = source.latitude;
		longitude = source.longitude;
		easting = source.easting;
		northing = source.northing;
		ref_latitude = source.ref_latitude;
		ref_longitude = source.ref_longitude;
		distance_algo = source.distance_algo;
		coordsystem = source.coordsystem;
		coordparam = source.coordparam;
		grid_i = source.grid_i;
		grid_j = source.grid_j;
		grid_k = source.grid_k;
		initializeMaps();
		setFunctionPointers();
	}
	return *this;
}

namespace mio {
/**
* @brief Print the content of the Coords object (usefull for debugging)
* The Coords is bound by "<Coords>" and "</Coords>" on separate lines
*/
	std::ostream& operator<<(std::ostream &os, const Coords& coord)
	{
		os << "<Coords>\n";
		os << "Lat/Long\t" << coord.printLatLon() << "\n";
		os << "X/Y_coords\t" << "(" << coord.getEasting() << " , " << coord.getNorthing() << ")" << "\n";
		os << "I/J_indices\t" << "(" << coord.getGridI() << " , " << coord.getGridJ() << ")" << "\n";
		os << "Projection\t" << coord.coordsystem << " " << coord.coordparam << "\n";
		os << "</Coords>\n";
		return os;
	}
}

/**
* @brief Default constructor
* This constructor builds a dummy object that performs no conversions but can be used for comparison
* purpose. This is more or less the equivalent of NULL for a pointer...
*/
Coords::Coords() {
	initializeMaps();
	setDefaultValues();
	setProj("NULL", "NULL");
}

/**
* @brief Regular constructor: usually, this is the constructor to use
* @param[in] _coordinatesystem string identifying the coordinate system to use
* @param[in] _parameters string giving some additional parameters for the projection (optional)
*
* See setProj() for a full description of these strings
*/
Coords::Coords(const std::string& _coordinatesystem, const std::string& _parameters) {
	setDefaultValues();
	initializeMaps();
	setProj(_coordinatesystem, _parameters);
}

/**
* @brief Local projection onstructor: this constructor is only suitable for building a local projection.
* Such a projection defines easting and northing as the distance (in meters) to a reference point
* which coordinates have to be provided here.
* @param[in] _lat_ref latitude of the reference point
* @param[in] _long_ref longitude of the reference point
*/
Coords::Coords(const double& _lat_ref, const double& _long_ref)
{
	setDefaultValues();
	setLocalRef(_lat_ref, _long_ref);
	initializeMaps();
	setProj("LOCAL", "");
}

/**
* @brief Returns the East coordinate in the configured projection system
* @return easting
*/
double Coords::getEasting() const {
	return easting;
}

/**
* @brief Returns the North coordinate in the configured projection system
* @return northing
*/
double Coords::getNorthing() const {
	return northing;
}

/**
* @brief Returns the Latitude in the configured projection system
* @return latitude
*/
double Coords::getLat() const {
	return latitude;
}

/**
* @brief Returns the Latitude in the configured projection system
* @return longitude
*/
double Coords::getLon() const {
	return longitude;
}

/**
* @brief Returns the Altitude. This is currently independent of the configured projection system
* @return altitude
*/
double Coords::getAltitude() const {
	return altitude;
}

/**
* @brief Returns the grid index along the X axis
* @return grid index i
*/
int Coords::getGridI() const {
	return grid_i;
}

/**
* @brief Returns the grid index along the Y axis
* @return grid index j
*/
int Coords::getGridJ() const {
	return grid_j;
}

/**
* @brief Returns the grid index along the Z axis
* @return grid index k
*/
int Coords::getGridK() const {
	return grid_k;
}

/**
* @brief Returns the projection parameters
* @param[out] proj_type projection type
* @param[out] proj_args optional arguments
*/
void Coords::getProj(std::string& proj_type, std::string& proj_args) const {
	proj_type = coordsystem;
	if(coordsystem=="LOCAL") {
		std::stringstream dms;
		dms << "(" << decimal_to_dms(ref_latitude) << " , " << decimal_to_dms(ref_longitude) << ")";
		proj_args=dms.str();
	} else {
		proj_args = coordparam;
	}
}

/**
* @brief Print a nicely formatted lat/lon in degrees, minutes, seconds
* @return lat/lon
*/
std::string Coords::printLatLon() const {
	std::stringstream dms;
	dms << "(" << decimal_to_dms(latitude) << " , " << decimal_to_dms(longitude) << ")";

	return dms.str();
}

/**
* @brief Set latitude and longitude
* The automatic update of the easting/northing can be turned off so that
* both lat/lon and east/north coordinates can be provided in order to thereafter check the
* coordinates by calling check().
* @param[in] _coordinates string containing the lat/lon to read
* @param[in] _altitude altitude to set (optional)
* @param[in] _update should the easting/northing be updated? (default=true)
*/
void Coords::setLatLon(const std::string& _coordinates, const double _altitude, const bool _update) {
	double lat, lon;
	parseLatLon(_coordinates, lat, lon);
	setLatLon(lat, lon, _altitude, _update);
}

/**
* @brief Set latitude and longitude
* The automatic update of the easting/northing can be turned off so that
* both lat/lon and east/north coordinates can be provided in order to thereafter check the
* coordinates by calling check().
* @param[in] _latitude latitude to set
* @param[in] _longitude longitude to set
* @param[in] _altitude altitude to set
* @param[in] _update should the easting/northing be updated? (default=true)
*/
void Coords::setLatLon(const double _latitude, const double _longitude, const double _altitude, const bool _update) {
	latitude = _latitude;
	longitude = _longitude;
	if(_altitude!=IOUtils::nodata) {
		altitude = _altitude;
	}
	if(coordsystem!="NULL" && _update==true) {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
	grid_i = grid_j = grid_k = IOUtils::inodata;
}

/**
* @brief Set easting and northing
* The automatic update of the latitude/longitude can be turned off so that
* both lat/lon and east/north coordinates can be provided in order to thereafter check the
* coordinates by calling check().
* @param[in] _easting easting to set
* @param[in] _northing northing to set
* @param[in] _altitude altitude to set
* @param[in] _update should the easting/northing be updated? (default=true)
*/
void Coords::setXY(const double _easting, const double _northing, const double _altitude, const bool _update) {
	easting = _easting;
	northing = _northing;
	if(_altitude!=IOUtils::nodata) {
		altitude = _altitude;
	}
	if(coordsystem!="NULL" && _update==true) {
		convert_to_WGS84(easting, northing, latitude, longitude);
	}
	grid_i = grid_j = grid_k = IOUtils::inodata;
}

/**
* @brief Set grid indices
* This index represent the position in a cartesian grid. It can not be automatically matched with
* a set of geographic coordinates because it needs the information about the said grid.
* Therefore, the coordinate object needs to be given to a grid object that will either set (i,j) or
* (lat,lon)/(easting,northing) as well as the grid's matching altitude if it is not already set.
* Any subsequent change of either (lat,lon) or (easting,northing) will reset these indexes to IOUtils::inodata.
* By default, setting (i,j) will invalidate (ie: delete) ALL geographic coordinates for the object (since we can not
* convert from grid indices to/from geographic coordinates in the current object by lack of information).
* Finally, the given indices are <b>NOT checked for validity</b>: such check must be done
* by calling Grid2DObject::gridify or Grid3DObject::gridify .
*
* @note To make it short: <b>use this method with caution!!</b>
* @param[in] _grid_i grid index along the X direction
* @param[in] _grid_j grid index along the Y direction
* @param[in] _grid_k grid index along the Z direction
* @param[in] _invalidate should the geographic coordinates be invalidated? (default=true, this flag should be false ONLY when calling from Grid2/3DObject)
*/
void Coords::setGridIndex(const int _grid_i, const int _grid_j, const int _grid_k, const bool _invalidate) {
//HACK TODO make grid_i,j,k friends of Grid2/3DObject -> remove _invalidate and ALWAYS invalidate
	grid_i = _grid_i;
	grid_j = _grid_j;
	grid_k = _grid_k;
	if(_invalidate==true) {
		latitude = IOUtils::nodata;
		longitude = IOUtils::nodata;
		easting = IOUtils::nodata;
		northing = IOUtils::nodata;
		altitude = IOUtils::nodata;
	}
}

/**
* \anchor Coordinate_types
* @brief Set projection to use
* This projection will be used for converting between lat/lon and East/North
* @param[in] _coordinatesystem string identifying the coordinate system to use
* @param[in] _parameters string giving some additional parameters for the projection (optional)
*
* The coordinate system can be any of the following:
* - CH1903 for coordinates in the Swiss Grid [ref: http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf]
* - UTM for UTM coordinates (the zone must be specified in the parameters, for example 31T) [ref: http://www.oc.nps.edu/oc2902w/maps/utmups.pdf]
* - PROJ4 for coordinate conversion relying on the Proj4 library [ref: http://trac.osgeo.org/proj/]
* - LOCAL for local coordinate system (using the horizontal and vertical distance from a reference point, see Coords::geo_distances for the available choice of distance algorithms)
*/
void Coords::setProj(const std::string& _coordinatesystem, const std::string& _parameters) {
	//the latitude/longitude had not been calculated, so we do it first in order to have our reference
	//before further conversions (usage scenario: giving a x,y and then converting to anyother x,y in another system
	if ((coordsystem != "NULL") && ((latitude==IOUtils::nodata) || (longitude==IOUtils::nodata))) {
		convert_to_WGS84(easting, northing, latitude, longitude);
	}

	if(_coordinatesystem == "") {
		coordsystem = std::string("NULL");
	} else {
		coordsystem = _coordinatesystem;
	}
	coordparam  = _parameters;
	setFunctionPointers();

	//since lat/long is our reference, we refresh x,y (only if lat/lon exist)
	if(latitude!=IOUtils::nodata && longitude!=IOUtils::nodata) {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Set the local projection reference coordinates
* This projection will be used for converting between lat/lon and East/North
* @param[in] _ref_latitude latitude of the local origin
* @param[in] _ref_longitude longitude of the local origin
*/
void Coords::setLocalRef(const double _ref_latitude, const double _ref_longitude) {
	if(_ref_latitude==IOUtils::nodata || _ref_longitude==IOUtils::nodata) {
		throw InvalidArgumentException("For LOCAL projection, please provide both reference latitude and longitude!", AT);
	}
	ref_latitude = _ref_latitude;
	ref_longitude = _ref_longitude;
	if(coordsystem=="LOCAL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Set the local projection reference coordinates
* This projection will be used for converting between lat/lon and East/North
* @param[in] _coordparam string containing the (lat,lon) of the local origin
*/
void Coords::setLocalRef(const std::string _coordparam) {
	coordparam = _coordparam;
	parseLatLon(coordparam, ref_latitude, ref_longitude);
	if(coordsystem=="LOCAL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Set the algorithm to use to compute distances
* Various algorithm exist that offer various precision/complexity tradeoffs.
* @param[in] _algo enum giving the algorithm to be used (see documentation for geo_distances)
*/
void Coords::setDistances(const geo_distances _algo) {
	distance_algo = _algo;
	if(coordsystem=="LOCAL") {
		convert_from_WGS84(latitude, longitude, easting, northing);
	}
}

/**
* @brief Check consistency of coordinates
* When both latitude/longitude and easting/northing are given, there is
* a risk of inconsistency between these two sets of coordinates.
* This method checks that enough information is available (ie: at least one set
* of coordinates is present) and if more than one is present, that it is consistent (within 5 meters)
* It throws and exception if something is not right.
*/
void Coords::check()
{
	//calculate/check coordinates if necessary
	if(coordsystem=="LOCAL" && (ref_latitude==IOUtils::nodata || ref_longitude==IOUtils::nodata)) {
		throw InvalidArgumentException("please define a reference point for LOCAL coordinate system", AT);
	}

	if(latitude==IOUtils::nodata || longitude==IOUtils::nodata) {
		if(easting==IOUtils::nodata || northing==IOUtils::nodata) {
			throw InvalidArgumentException("missing positional parameters (easting,northing) or (lat,long) for coordinate", AT);
		}
		convert_to_WGS84(easting, northing, latitude, longitude);
	} else {
		if(easting==IOUtils::nodata || northing==IOUtils::nodata) {
			convert_from_WGS84(latitude, longitude, easting, northing);
		} else {
			double tmp_lat, tmp_lon;
			convert_to_WGS84(easting, northing, tmp_lat, tmp_lon);

			if(!IOUtils::checkEpsilonEquality(latitude, tmp_lat, IOUtils::lat_epsilon) || !IOUtils::checkEpsilonEquality(longitude, tmp_lon, IOUtils::lon_epsilon)) {
				throw InvalidArgumentException("Latitude/longitude and xllcorner/yllcorner don't match for coordinate", AT);
			}
		}
	}
}

/**
* @brief Calculate the distance between two points
* @param destination destination coordinate
* @return distance in meters
*/
double Coords::distance(const Coords& destination) const {
	double dist, bearing;
	distance(destination, dist, bearing);
	return dist;
}

/**
* @brief Check if two Coords object are using the same projection
* @param target coordinate to compare to
* @return true or false
*/
bool Coords::isSameProj(const Coords& target) const {
	if(coordsystem=="LOCAL") {
		return ( target.coordsystem=="LOCAL" && ref_latitude==target.ref_latitude && ref_longitude==target.ref_longitude );
	} else {
		return ( coordsystem==target.coordsystem && coordparam==target.coordparam );
	}
}

/**
* @brief Copy the projection parameters of another Coords object
* @param source source object to copy the projection from
* @param _update should the necessary coordinates be updated? (default=true)
*/
void Coords::copyProj(const Coords& source, const bool _update) {
	if(!isSameProj(source)) {
		//we only do a copy if we are not already using the same projection
		if(source.coordsystem=="LOCAL") {
			coordsystem="LOCAL";
			coordparam=source.coordparam;
			ref_latitude=source.ref_latitude;
			ref_longitude=source.ref_longitude;
		} else {
			coordsystem=source.coordsystem;
			coordparam=source.coordparam;
		}
		setFunctionPointers();

		if(_update==true) {
			if((latitude!=IOUtils::nodata) && (longitude!=IOUtils::nodata)) {
				convert_from_WGS84(latitude, longitude, easting, northing);
			} else {
				convert_to_WGS84(easting, northing, latitude, longitude);
			}
		}
	}
}

/////////////////////////////////////////////////////private methods
/**
* @brief Method converting towards WGS84
* @param[in] easting easting of the coordinate to convert
* @param[in] northing northing of the coordinate to convert
* @param[out] latitude converted latitude
* @param[out] longitude converted longitude
*/
void Coords::convert_to_WGS84(double easting, double northing, double& latitude, double& longitude) const
{
	if((easting!=IOUtils::nodata) && (northing!=IOUtils::nodata)) {
		(this->*convToWGS84)(easting, northing, latitude, longitude);
	} else {
		latitude = IOUtils::nodata;
		longitude = IOUtils::nodata;
	}
}

/**
* @brief Method converting towards WGS84
* @param[in] latitude latitude of the coordinate to convert
* @param[in] longitude longitude of the coordinate to convert
* @param[out] easting converted easting
* @param[out] northing converted northing
*/
void Coords::convert_from_WGS84(double latitude, double longitude, double& easting, double& northing) const
{
	if((latitude!=IOUtils::nodata) && (longitude!=IOUtils::nodata)) {
		(this->*convFromWGS84)(latitude, longitude, easting, northing);
	} else {
		easting = IOUtils::nodata;
		northing = IOUtils::nodata;
	}
}

/**
* @brief Parse a latitude or longitude
* It can be formatted as any of the following examples:
* - 46° 48' 03" (with or without spaces, decimal or integer numbers)
* - 46d 48' 03" (with or without spaces, decimal or integer numbers)
* - 46 48' 03" (with spaces, decimal or integer numbers)
* - 46° 48.02'(with or without spaces, decimal or integer numbers)
* - 46d 48.02'(with or without spaces, decimal or integer numbers)
* - 46 48.02'(with spaces, decimal or integer numbers)
* - 46.802°
* - 46.802d
* - 46.802
* @param[in] dms string containing the coordinate
* @return coordinate in decimal
*/
double Coords::dms_to_decimal(const std::string& dms) {
	double d=IOUtils::nodata, m=IOUtils::nodata, s=IOUtils::nodata, decimal=IOUtils::nodata;

	if 	((sscanf(dms.c_str(), "%lf°%lf'%lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lf° %lf' %lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lfd%lf'%lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lfd %lf' %lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lf %lf' %lf\"", &d, &m ,&s) < 3) &&
		(sscanf(dms.c_str(), "%lf° %lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lf°%lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lfd %lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lfd%lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lf %lf'", &d, &m) < 2) &&
		(sscanf(dms.c_str(), "%lf°", &d) < 1) &&
		(sscanf(dms.c_str(), "%lfd", &d) < 1) &&
		(sscanf(dms.c_str(), "%lf", &d) < 1)) {
			throw InvalidFormatException("Can not parse given latitude or longitude: "+dms,AT);
	}

	decimal = d;
	if(m!=IOUtils::nodata) {
		decimal += m/60.;
	}
	if(s!=IOUtils::nodata) {
		decimal += s/3600.;
	}

	return decimal;
}

/**
* @brief Parse a latitude-longitude pair
* It can be formatted as any of the following examples:
* - lat lon (without any spaces in the latitude or longitude string)
* - lat/lon
* - (lat;lon)
* - (lat,lon)
* @param[in] coordinates string containing the coordinates
* @param[out] lat parsed latitude
* @param[out] lon parsed longitude
*/
void Coords::parseLatLon(const std::string& coordinates, double&
lat, double& lon) {
	char lat_str[32]="";
	char lon_str[32]="";

	if(coordinates.size()>(32+32)) {
		throw InvalidFormatException("Given lat/lon string is too long! ",AT);
	}

	if 	((sscanf(coordinates.c_str(), "%[0-9.,°d'\"-] %[0-9.,°d'\"-]", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "%[0-9.,°d'\"- ]/%[0-9.,°d'\"- ]", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "(%[0-9.,°d'\"- ];%[0-9.,°d'\"- ])", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "(%[0-9.°d'\"- ],%[0-9.°d'\"- ])", lat_str, lon_str) < 2)) {
			throw InvalidFormatException("Can not parse given lat/lon: "+coordinates,AT);
	}

	lat = dms_to_decimal(std::string(lat_str));
	lon = dms_to_decimal(std::string(lon_str));
}

/**
* @brief Converts a decimal latitude or longitude to degrees, minutes, seconds
* It formats its arguments as in the following example: 46°48'03"
* @param[in] decimal decimal coordinate to convert
* @return string containing the formatted coordinate
*/
std::string Coords::decimal_to_dms(const double& decimal) {
	std::stringstream dms;
	const int d = (int)floor(decimal);
	const int m = (int)floor( (decimal - (double)d)*60. );
	const int s = (int)floor( decimal - (double)d - (double)m/60. );

	dms << d << "°" << m << "'" << s << "\"";
	return dms.str();
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to Swiss grid
* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
* @param[in] lat_in Decimal Latitude
* @param[in] long_in Decimal Longitude
* @param[out] east_out easting coordinate (Swiss system)
* @param[out] north_out northing coordinate (Swiss system)
*/
void Coords::WGS84_to_CH1903(double lat_in, double long_in, double& east_out, double& north_out) const
{
	//converts WGS84 coordinates (lat,long) to the Swiss coordinates. See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf
	//The elevation is supposed to be above sea level, so it does not require any conversion
	//lat and long must be decimal (and they will be converted to seconds)
	const double phi_p = (lat_in*3600. - 169028.66) / 10000.;
	const double lambda_p = (long_in*3600. - 26782.5) / 10000.;

	east_out = 600072.37
		+ 211455.93	* lambda_p
		- 10938.51	* lambda_p * phi_p
		- 0.36		* lambda_p * (phi_p*phi_p)
		- 44.54		* (lambda_p*lambda_p*lambda_p);

	north_out = 200147.07
		+ 308807.95	* phi_p
		+ 3745.25	* (lambda_p*lambda_p)
		+ 76.63		* (phi_p*phi_p)
		- 194.56	* (lambda_p*lambda_p) * phi_p
		+ 119.79	* (phi_p*phi_p*phi_p);

	/*// if necessary for the elevation, uncomment this block
	h_out = h_in - 49.55
		+ 2.73		* lambda_p
		+ 6.94		* phi_p;
	*/
}

/**
* @brief Coordinate conversion: from Swiss grid to WGS84 Lat/Long
* See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf for more.
* @param[in] east_in easting coordinate (Swiss system)
* @param[in] north_in northing coordinate (Swiss system)
* @param[out] lat_out Decimal Latitude
* @param[out] long_out Decimal Longitude
*/
void Coords::CH1903_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
	//converts Swiss coordinates to WGS84 coordinates (lat,long). See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf
	//The elevation is supposed to be above sea level, so it does not require any conversion
	//lat and long are decimal
	const double y_p = (east_in - 600000.) / 1000000.;
	const double x_p = (north_in - 200000.) / 1000000.;

	const double lambda_p = 2.6779094
		+ 4.728982	* y_p
		+ 0.791484	* y_p * x_p
		+ 0.1306	* y_p * (x_p*x_p)
		- 0.0436	* (y_p*y_p*y_p);

	const double phi_p = 16.9023892
		+ 3.238272	* x_p
		- 0.270978	* (y_p*y_p)
		- 0.002528	* (x_p*x_p)
		- 0.0447	* (y_p*y_p) * x_p
		- 0.0140	* (x_p*x_p*x_p);

	lat_out = phi_p * 100./36.;
	long_out = lambda_p * 100./36.;

	/*// if necessary for the elevation, uncomment this block
	h_out = h_in + 49.55
		- 12.60		* y_p
		- 22.64		* x_p;
	*/
}

int Coords::getUTMZone(const double _latitude, const double _longitude, std::string& zone_out) const
{//This routine determines the correct UTM letter designator for the given latitude
//UTM limits its coverage to [80S , 84N], outside of this, returns Y/Z/A/B for the zone

	//computing zone number, assuming longitude in [-180. ; 180[
	int ZoneNumber = int((_longitude + 180.)/6.) + 1;

	// Special zones for Scandinavia
	if( _latitude >= 72.0 && _latitude < 84.0 ) {
		if(      _longitude >= 0.0  && _longitude <  9.0 ) ZoneNumber = 31;
		else if( _longitude >= 9.0  && _longitude < 21.0 ) ZoneNumber = 33;
		else if( _longitude >= 21.0 && _longitude < 33.0 ) ZoneNumber = 35;
		else if( _longitude >= 33.0 && _longitude < 42.0 ) ZoneNumber = 37;
	 }
	if( latitude >= 56.0 && _latitude < 64.0 && _longitude >= 3.0 && _longitude < 12.0 ) {
		ZoneNumber = 32;
	}

	//getting zone letter
	char zoneLetter='Z';
	if     ((0 >= _longitude) && (_latitude >  84)) zoneLetter = 'Y';
	else if((0 <  _longitude) && (_latitude >  84)) zoneLetter = 'Z';
	else if((84 >= _latitude) && (_latitude >= 72)) zoneLetter = 'X';
	else if((72 > _latitude) && (_latitude >= 64)) zoneLetter = 'W';
	else if((64 > _latitude) && (_latitude >= 56)) zoneLetter = 'V';
	else if((56 > _latitude) && (_latitude >= 48)) zoneLetter = 'U';
	else if((48 > _latitude) && (_latitude >= 40)) zoneLetter = 'T';
	else if((40 > _latitude) && (_latitude >= 32)) zoneLetter = 'S';
	else if((32 > _latitude) && (_latitude >= 24)) zoneLetter = 'R';
	else if((24 > _latitude) && (_latitude >= 16)) zoneLetter = 'Q';
	else if((16 > _latitude) && (_latitude >= 8)) zoneLetter = 'P';
	else if(( 8 > _latitude) && (_latitude >= 0)) zoneLetter = 'N';
	else if(( 0 > _latitude) && (_latitude >= -8)) zoneLetter = 'M';
	else if((-8 > _latitude) && (_latitude >= -16)) zoneLetter = 'L';
	else if((-16 > _latitude) && (_latitude >= -24)) zoneLetter = 'K';
	else if((-24 > _latitude) && (_latitude >= -32)) zoneLetter = 'J';
	else if((-32 > _latitude) && (_latitude >= -40)) zoneLetter = 'H';
	else if((-40 > _latitude) && (_latitude >= -48)) zoneLetter = 'G';
	else if((-48 > _latitude) && (_latitude >= -56)) zoneLetter = 'F';
	else if((-56 > _latitude) && (_latitude >= -64)) zoneLetter = 'E';
	else if((-64 > _latitude) && (_latitude >= -72)) zoneLetter = 'D';
	else if((-72 > _latitude) && (_latitude >= -80)) zoneLetter = 'C';
	else if((0 >= _longitude) && (latitude <= -80)) zoneLetter = 'A';
	else if((0 <  _longitude) && (latitude <= -80)) zoneLetter = 'B';

	std::stringstream zone;
	zone << ZoneNumber << zoneLetter;
	zone_out = zone.str();

	return ZoneNumber;
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to UTM grid
* See http://www.oc.nps.edu/oc2902w/maps/utmups.pdf for more.
* @param[in] lat_in Decimal Latitude
* @param[in] long_in Decimal Longitude
* @param[out] east_out easting coordinate (Swiss system)
* @param[out] north_out northing coordinate (Swiss system)
*/
void Coords::WGS84_to_UTM(double lat_in, double long_in, double& east_out, double& north_out) const
{//converts WGS84 coordinates (lat,long) to UTM coordinates.
//See USGS Bulletin 1532 or http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf
//also http://www.uwgb.edu/dutchs/usefuldata/UTMFormulas.HTM
//also http://www.oc.nps.edu/oc2902w/maps/utmups.pdf or Chuck Gantz (http://www.gpsy.com/gpsinfo/geotoutm/)
	//Geometric constants
	const double to_rad = M_PI / 180.0;
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double eP2 = e2 / (1.-e2);	//second ellispoid eccentricity, squared (=(a²-b²)/b²)
	const double k0 = 0.9996;	//scale factor for the projection

	//getting posistion parameters
	std::string zone;
	long_in = fmod(long_in+360.+180., 360.) - 180.; //normalized to [-180. ; 180.[
	const double Long = long_in * to_rad;
	const double Lat = lat_in * to_rad;
	const int zoneNumber = getUTMZone(lat_in, long_in, zone);
	const double long0 = (double)((zoneNumber - 1)*6 - 180 + 3) * to_rad; //+3 puts origin in middle of zone

	//Geometrical parameters
	const double nu = a / sqrt(1.-e2*IOUtils::pow2(sin(Lat))); //radius of curvature of the earth perpendicular to the meridian plane
	const double p = (Long-long0);

	//calculating first the coefficients of the series, then the Meridional Arc M itself
	const double n = (a-b)/(a+b);
	const double n2=n*n, n3=n*n*n, n4=n*n*n*n, n5=n*n*n*n*n, n6=n*n*n*n*n*n;
	const double A = a           * (1. - n + 5./4.*(n2 - n3) + 81./64.*(n4 - n5));
	const double B = (3./2.*a)   * (n - n2 + 7./8.*(n3 - n4) + 55./64.*(n5 - n6));
	const double C = (15./16.*a) * (n2 - n3 + 3./4.*(n4 - n5));
	const double D = (35./48.*a) * (n3 - n4 + 11./16.*(n5 - n6));
	const double E = (315./512.*a) * (n4 - n5); //correctrion of ~0.03mm
	const double M = A*Lat - B*sin(2.*Lat) + C*sin(4.*Lat) - D*sin(6.*Lat) + E*sin(8.*Lat);

	//calculating the coefficients for the series
	const double K1 = M*k0;
	const double K2 = 1./4.*k0*nu*sin(2.*Lat);
	const double K3 = (k0*nu*sin(Lat)*IOUtils::pow3(cos(Lat))*1./24.) * (5. - IOUtils::pow2(tan(Lat)) + 9.*eP2*IOUtils::pow2(cos(Lat)) + 4.*eP2*eP2*IOUtils::pow4(cos(Lat)));
	const double K4 = k0*nu*cos(Lat);
	const double K5 = (k0*nu*IOUtils::pow3(cos(Lat))*1./6.) * (1. - IOUtils::pow2(tan(Lat)) + eP2*IOUtils::pow2(cos(Lat)));

	north_out = K1 + K2*p*p + K3*p*p*p*p;
	east_out = K4*p + K5*p*p*p + 500000.0;

	if(Lat < 0)
		north_out += 10000000.0; //offset for southern hemisphere
}

/**
* @brief Coordinate conversion: from UTM grid to WGS84 Lat/Long
* See http://www.oc.nps.edu/oc2902w/maps/utmups.pdf for more.
* @param[in] east_in easting coordinate (UTM)
* @param[in] north_in northing coordinate (UTM)
* @param[out] lat_out Decimal Latitude
* @param[out] long_out Decimal Longitude
*/
void Coords::UTM_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{//converts UTM coordinates to WGS84 coordinates (lat,long).
//See USGS Bulletin 1532 or http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf
//also http://www.uwgb.edu/dutchs/usefuldata/UTMFormulas.HTM
//also http://www.oc.nps.edu/oc2902w/maps/utmups.pdf or Chuck Gantz (http://www.gpsy.com/gpsinfo/geotoutm/)
	//Geometric constants
	const double to_deg = 180.0 / M_PI;
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double eP2 = e2 / (1.-e2);	//second ellispoid eccentricity, squared (=(a²-b²)/b²)
	const double k0 = 0.9996;		//scale factor for the projection

	//UTM Zone information
	int zoneNumber;
	char zoneLetter;
	if 	((sscanf(coordparam.c_str(), "%d%c", &zoneNumber, &zoneLetter) < 2) &&
		(sscanf(coordparam.c_str(), "%d %c)", &zoneNumber, &zoneLetter) < 2)) {
			throw InvalidFormatException("Can not parse given UTM zone: "+coordparam,AT);
	}
	zoneLetter = toupper(zoneLetter); //just in case... (sorry for the pun!)
	if(zoneLetter=='Y' || zoneLetter=='Z' || zoneLetter=='A' || zoneLetter=='B') {
			//Special zones for the poles: we should NOT use UTM in these regions!
			throw InvalidFormatException("Invalid UTM zone: "+coordparam+" (trying to use UTM in polar regions)",AT);
	}

	//set reference parameters: central meridian of the zone, true northing and easting
	//please note that the special zones still use the reference meridian as given by their zone number (ie: even if it might not be central anymore)
	const int long0 = (zoneNumber - 1)*6 - 180 + 3;  //+3 puts origin in "middle" of zone as required for the projection meridian (might not be the middle for special zones)
	if(zoneLetter<='N') {
		north_in -= 10000000.0; //offset used for southern hemisphere
	}
	east_in -= 500000.0; //longitude offset: x coordinate is relative to central meridian

	//calculating footprint latitude fp
	const double arc = north_in/k0; //Meridional arc
	const double mu = arc / (a*(1.-e2/4.-3.*e2*e2/64.-5.*e2*e2*e2/256.));
	const double e1 = (1.-b/a) / (1.+b/a); //simplification of [1 - (1 - e2)1/2]/[1 + (1 - e2)1/2]
	const double J1 = (3./2.*e1 - 27./32.*e1*e1*e1);
	const double J2 = (21./16.*e1*e1 - 55./32.*e1*e1*e1*e1);
	const double J3 = (151./96.*e1*e1*e1);
	const double J4 = (1097./512.*e1*e1*e1*e1);
	const double fp = mu + J1*sin(2.*mu) + J2*sin(4.*mu) + J3*sin(6.*mu) + J4*sin(8.*mu);

	//calculating the parameters
	const double C1 = eP2 * IOUtils::pow2(cos(fp));
	const double T1 = IOUtils::pow2( tan(fp) );
	const double R1 = a*(1.-e2) / pow((1.-e2*IOUtils::pow2(sin(fp))), 1.5);
	const double N1 = a / sqrt(1.-e2*IOUtils::pow2(sin(fp)));
	const double D = east_in / (N1*k0);

	//calculating the coefficients of the series for latitude and longitude
	const double Q1 = N1*tan(fp)/R1;
	const double Q2 = 0.5*D*D;
	const double Q3 = (5. + 3.*T1 + 10.*C1 - 4.*C1*C1 - 9.*eP2) * 1./24.*D*D*D*D;
	const double Q4 = (61. + 90.*T1 + 298.*C1 + 45.*T1*T1 - 3.*C1*C1 - 252.*eP2) * 1./720.*D*D*D*D*D*D;
	const double Q5 = D;
	const double Q6 = (1. + 2.*T1 + C1) * 1./6.*D*D*D;
	const double Q7 = (5. - 2.*C1 + 28.*T1 - 3.*C1*C1 + 8.*eP2 + 24.*T1*T1) * 1./120.*D*D*D*D*D;

	lat_out = (fp - Q1 * (Q2 - Q3 + Q4))*to_deg;
	long_out = (double)long0 + ((Q5 - Q6 + Q7)/cos(fp))*to_deg;
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to proj4 parameters
* @param lat_in Decimal Latitude
* @param long_in Decimal Longitude
* @param east_out easting coordinate (target system)
* @param north_out northing coordinate (target system)
*/
void Coords::WGS84_to_PROJ4(double lat_in, double long_in, double& east_out, double& north_out) const
{
#ifdef PROJ4
	const std::string src_param="+proj=latlong +datum=WGS84 +ellps=WGS84";
	projPJ pj_latlong, pj_dest;
	double x=long_in*DEG_TO_RAD, y=lat_in*DEG_TO_RAD;

	if ( !(pj_dest = pj_init_plus(coordparam.c_str())) ) {
		pj_free(pj_dest);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+coordparam, AT);
	}
	if ( !(pj_latlong = pj_init_plus(src_param.c_str())) ) {
		pj_free(pj_latlong);
		pj_free(pj_dest);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+src_param, AT);
	}

	const int p = pj_transform(pj_latlong, pj_dest, 1, 1, &x, &y, NULL );
	if(p!=0) {
		pj_free(pj_latlong);
		pj_free(pj_dest);
		throw ConversionFailedException("PROJ4 conversion failed: "+p, AT);
	}
	east_out = x;
	north_out = y;
	pj_free(pj_latlong);
	pj_free(pj_dest);
#else
	(void)lat_in;
	(void)long_in;
	(void)east_out;
	(void)north_out;
	throw IOException("Not compiled with PROJ4 support", AT);
#endif
}

/**
* @brief Coordinate conversion: from proj4 parameters to WGS84 Lat/Long
* @param east_in easting coordinate (Swiss system)
* @param north_in northing coordinate (Swiss system)
* @param lat_out Decimal Latitude
* @param long_out Decimal Longitude
*/
void Coords::PROJ4_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
#ifdef PROJ4
	const std::string dest_param="+proj=latlong +datum=WGS84 +ellps=WGS84";
	projPJ pj_latlong, pj_src;
	double x=east_in, y=north_in;

	if ( !(pj_src = pj_init_plus(coordparam.c_str())) ) {
		pj_free(pj_src);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+coordparam, AT);
	}
	if ( !(pj_latlong = pj_init_plus(dest_param.c_str())) ) {
		pj_free(pj_latlong);
		pj_free(pj_src);
		throw InvalidArgumentException("Failed to initalize Proj4 with given arguments: "+dest_param, AT);
	}

	const int p = pj_transform(pj_src, pj_latlong, 1, 1, &x, &y, NULL );
	if(p!=0) {
		pj_free(pj_latlong);
		pj_free(pj_src);
		throw ConversionFailedException("PROJ4 conversion failed: "+p, AT);
	}
	long_out = x*RAD_TO_DEG;
	lat_out = y*RAD_TO_DEG;
	pj_free(pj_latlong);
	pj_free(pj_src);
#else
	(void)east_in;
	(void)north_in;
	(void)lat_out;
	(void)long_out;
	throw IOException("Not compiled with PROJ4 support", AT);
#endif
}

void Coords::distance(const Coords& destination, double& distance, double& bearing) const {
//HACK: this is the 2D distance, it does not work in 3D!!
	if(isSameProj(destination)) {
		//we can use simple cartesian grid arithmetic
		const double to_deg = 180.0 / M_PI;
		distance = sqrt( IOUtils::pow2(easting - destination.getEasting()) + IOUtils::pow2(northing - destination.getNorthing()) );
		bearing = atan2( northing - destination.getNorthing() , easting - destination.getEasting() );
		bearing = fmod( bearing*to_deg+360. , 360. );
	} else {
		switch(distance_algo) {
			case GEO_COSINE:
				distance = cosineDistance(latitude, longitude, destination.getLat(), destination.getLon(), bearing);
				break;
			case GEO_VINCENTY:
				distance = VincentyDistance(latitude, longitude, destination.getLat(), destination.getLon(), bearing);
				break;
			default:
				throw InvalidArgumentException("Unrecognized geodesic distance algorithm selected", AT);
		}
	}
}

/**
* @brief Coordinate conversion: from WGS84 Lat/Long to local grid as given in coordparam
* @param lat_in Decimal Latitude
* @param long_in Decimal Longitude
* @param east_out easting coordinate (target system)
* @param north_out northing coordinate (target system)
*/
void Coords::WGS84_to_local(double lat_in, double long_in, double& east_out, double& north_out) const
{
	double alpha;
	const double to_rad = M_PI / 180.0;
	double distance;

	if((ref_latitude==IOUtils::nodata) || (ref_longitude==IOUtils::nodata)) {
		east_out = IOUtils::nodata;
		north_out = IOUtils::nodata;
		//throw InvalidArgumentException("No reference coordinate provided for LOCAL projection", AT);
	} else {
		switch(distance_algo) {
			case GEO_COSINE:
				distance = cosineDistance(ref_latitude, ref_longitude, lat_in, long_in, alpha);
				break;
			case GEO_VINCENTY:
				distance = VincentyDistance(ref_latitude, ref_longitude, lat_in, long_in, alpha);
				break;
			default:
				throw InvalidArgumentException("Unrecognized geodesic distance algorithm selected", AT);
		}

		east_out = distance*sin(alpha*to_rad);
		north_out = distance*cos(alpha*to_rad);
	}
}


/**
* @brief Coordinate conversion: from local grid as given in coordparam to WGS84 Lat/Long
* @param east_in easting coordinate (Swiss system)
* @param north_in northing coordinate (Swiss system)
* @param lat_out Decimal Latitude
* @param long_out Decimal Longitude
*/
void Coords::local_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
	const double to_deg = 180.0 / M_PI;
	const double distance = sqrt( IOUtils::pow2(east_in) + IOUtils::pow2(north_in) );
	const double bearing = fmod( atan2(east_in, north_in)*to_deg+360. , 360.);

	if((ref_latitude==IOUtils::nodata) || (ref_longitude==IOUtils::nodata)) {
		lat_out = IOUtils::nodata;
		long_out = IOUtils::nodata;
		//throw InvalidArgumentException("No reference coordinate provided for LOCAL projection", AT);
	} else {
		switch(distance_algo) {
			case GEO_COSINE:
				cosineInverse(ref_latitude, ref_longitude, distance, bearing, lat_out, long_out);
				break;
			case GEO_VINCENTY:
				VincentyInverse(ref_latitude, ref_longitude, distance, bearing, lat_out, long_out);
				break;
			default:
				throw InvalidArgumentException("Unrecognized geodesic distance algorithm selected", AT);
		}
	}
}

void Coords::NULL_to_WGS84(double /*east_in*/, double /*north_in*/, double& /*lat_out*/, double& /*long_out*/) const
{
	throw InvalidArgumentException("The projection has not been initialized!", AT);
}

void Coords::WGS84_to_NULL(double /*lat_in*/, double /*long_in*/, double& /*east_out*/, double& /*north_out*/) const
{
	throw InvalidArgumentException("The projection has not been initialized!", AT);
}

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
void Coords::cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon) const
{
	const double Rearth = 6371.e3;
	const double to_rad = M_PI / 180.0;
	const double to_deg = 180.0 / M_PI;
	const double lat_ref_rad = lat_ref*to_rad;
	const double bearing_rad = bearing*to_rad;

	lat = asin( sin(lat_ref_rad)*cos(distance/Rearth) +
				cos(lat_ref_rad)*sin(distance/Rearth)*cos(bearing_rad) );
	lon = lon_ref*to_rad + atan2( sin(bearing_rad)*sin(distance/Rearth)*cos(lat_ref_rad) ,
					cos(distance/Rearth) - sin(lat_ref_rad)*sin(lat) );
	lon = fmod(lon+M_PI, 2.*M_PI) - M_PI;

	lat *= to_deg;
	lon *= to_deg;
}

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
double Coords::cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha) const
{
	const double Rearth = 6371.e3;
	const double to_rad = M_PI / 180.0;
	const double d = acos(
		sin(lat1*to_rad) * sin(lat2*to_rad)
		+ cos(lat1*to_rad) * cos(lat2*to_rad) * cos((lon2-lon1)*to_rad)
		) * Rearth;

	alpha = atan2( sin((lon2-lon1)*to_rad)*cos(lat2*to_rad) ,
			cos(lat1*to_rad)*sin(lat2*to_rad) - sin(lat1*to_rad)*cos(lat2*to_rad)*cos((lon2-lon1)*to_rad)
 		) / to_rad;
	alpha = fmod((alpha+360.), 360.);

	return d;
}

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
double Coords::VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha) const
{
	const double thresh = 1.e-12;	//convergence absolute threshold
	const int n_max = 100;		//maximum number of iterations
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = M_PI / 180.0;

	const double L = (lon1 - lon2)*to_rad;
	const double U1 = atan( (1.-f)*tan(lat1*to_rad) );
	const double U2 = atan( (1.-f)*tan(lat2*to_rad) );

	double lambda = L, lambda_p=0., delta_sigma;
	double sin_sigma, cos_sigma, sigma, sin_alpha, cos_alpha2, cos_2sigma_m;
	double C, u2, A, B, s;
	int n=0;
	do {
		sin_sigma = sqrt( IOUtils::pow2(cos(U2)*sin(lambda)) + IOUtils::pow2(cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)) );
		if(sin_sigma==0.) {
			//co-incident points
			return 0.;
		}
		cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
		sigma = atan2(sin_sigma,cos_sigma);
		sin_alpha = cos(U1)*cos(U2)*sin(lambda) / sin_sigma;
		cos_alpha2 = 1. - IOUtils::pow2(sin_alpha);
		if(lat1==0. && lat2==0.) {
			cos_2sigma_m = 0.;
		} else {
			cos_2sigma_m = cos_sigma - 2.*sin(U1)*sin(U2)/cos_alpha2;
		}
		C = f/16. * cos_alpha2*(4.+f*(4.-3.*cos_alpha2));
		lambda_p = lambda;
		lambda = L + (1.-C)*f*sin_alpha*(
			sigma + C*sin_sigma*( cos_2sigma_m + C * cos_sigma * (-1.+2.*IOUtils::pow2(cos_2sigma_m)) )
			);
		n++;
	} while ( (n<n_max) && (fabs(lambda - lambda_p) > thresh) );

	if(n>n_max) {
		throw IOException("Distance calculation not converging", AT);
	}

	u2 = cos_alpha2 * (a*a - b*b) / (b*b);
	A = 1. + u2/16384. * ( 4096.+u2*(-768.+u2*(320.-175.*u2)) );
	B = u2/1024. * ( 256.+u2*(-128.+u2*(74.-47.*u2)) );
	delta_sigma = B*sin_sigma*( cos_2sigma_m+B/4.*( cos_sigma*(-1.+2.*IOUtils::pow2(cos_2sigma_m)) - B/6.*(cos_2sigma_m*(-3.+4.*IOUtils::pow2(sin_sigma))*(-3.+4.*IOUtils::pow2(cos_2sigma_m))) ) );

	s = b*A*(sigma - delta_sigma);	//distance between the two points

	//computation of the average forward bearing
	double alpha1 = atan2(cos(U2)*sin(lambda), cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda)) / to_rad; //forward azimuth
	//double alpha2 = atan2(cos(U1)*sin(lambda), sin(U1)*cos(U2)-cos(U1)*sin(U2)*cos(lambda)) / to_rad; //reverse azimuth

	//trying to get a normal compass bearing... TODO: make sure it works and understand why
	alpha1 = fmod(-alpha1+360., 360.);
	//alpha2 = fmod(alpha2+180., 360.);

	//we only keep the forward bearing, otherwise the reverse projection will not produce the initial point
	alpha = alpha1;
	return s;
}

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
void Coords::VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon) const
{//well, actually this is the DIRECT Vincenty formula
	const double thresh = 1.e-12;	//convergence absolute threshold
	const double a = ellipsoids[E_WGS84].a;	//major ellipsoid semi-axis, value for wgs84
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis, value for wgs84
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = M_PI / 180.0;
	const double to_deg = 180.0 / M_PI;

	const double alpha1 = bearing*to_rad;
	const double tanU1 = (1.-f)*tan(lat_ref*to_rad);
	const double cosU1 = 1./sqrt(1.+tanU1*tanU1);
	const double sinU1 = tanU1*cosU1;
	const double sigma1 = atan2(tanU1,cos(alpha1));
	const double sinAlpha = cosU1*sin(alpha1);
	const double cos2alpha = 1. - sinAlpha*sinAlpha;
	const double u2 = cos2alpha * (a*a - b*b) / (b*b);
	const double A = 1. + u2/16384. * (4096. + u2*(-768.+u2*(320.-175.*u2)) );
	const double B = u2/1024. * (256. + u2*(-128.+u2*(74.-47.*u2)));

	double sigma = distance / (b*A);
	double sigma_p = 2.*M_PI;
	double cos2sigma_m = cos( 2.*sigma1 + sigma ); //required to avoid uninitialized value

	while (fabs(sigma - sigma_p) > thresh) {
		cos2sigma_m = cos( 2.*sigma1 + sigma );
		double delta_sigma = B*sin(sigma) * ( cos2sigma_m + B/4. * (
			cos(sigma)*(-1.+2.*cos2sigma_m*cos2sigma_m)
			-B/6. * cos2sigma_m * (-3.+4.*IOUtils::pow2(sin(sigma))) * (-3.+4.*cos2sigma_m*cos2sigma_m)
			) );
		sigma_p = sigma;
		sigma = distance / (b*A) + delta_sigma;
	}

	lat = atan2( sinU1*cos(sigma) + cosU1*sin(sigma)*cos(alpha1),
		     (1.-f) * sqrt( sinAlpha*sinAlpha + IOUtils::pow2(sinU1*sin(sigma) - cosU1*cos(sigma)*cos(alpha1)) )
		   );
	const double lambda = atan2( sin(sigma)*sin(alpha1), cosU1*cos(sigma) - sinU1*sin(sigma)*cos(alpha1) );
	const double C = f/16. * cos2alpha * (4.+f*(4.-3.*cos2alpha));
	const double L = lambda - (1.-C) * f * sinAlpha * (
				sigma + C * sin(sigma) * ( cos2sigma_m+C*cos(sigma) * (-1.+2.*cos2sigma_m*cos2sigma_m) )
				);

	lat = lat * to_deg;
	lon = lon_ref + (L*to_deg);
	//const double alpha2 = atan2( sinAlpha, -(sinU1*sin(sigma)-cosU1*cos(sigma)*cos(alpha1)) ); //reverse azimuth
}

void Coords::initializeMaps() {
//Please don't forget to mirror the keywords here in the documentation in Coords.h!!
	to_wgs84["CH1903"]   = &Coords::CH1903_to_WGS84;
	from_wgs84["CH1903"] = &Coords::WGS84_to_CH1903;
	to_wgs84["UTM"]   = &Coords::UTM_to_WGS84;
	from_wgs84["UTM"] = &Coords::WGS84_to_UTM;
	to_wgs84["PROJ4"]   = &Coords::PROJ4_to_WGS84;
	from_wgs84["PROJ4"] = &Coords::WGS84_to_PROJ4;
	to_wgs84["LOCAL"]   = &Coords::local_to_WGS84;
	from_wgs84["LOCAL"] = &Coords::WGS84_to_local;
	to_wgs84["NULL"]   = &Coords::NULL_to_WGS84;
	from_wgs84["NULL"] = &Coords::WGS84_to_NULL;
}

void Coords::setFunctionPointers() {
	//check whether there exists a tranformation for the given coordinatesystem
	//init function pointers
	std::map<std::string, convfunc>::iterator mapitTo;
	std::map<std::string, convfunc>::iterator mapitFrom;
	mapitTo   = to_wgs84.find(coordsystem);
	mapitFrom = from_wgs84.find(coordsystem);

	if ((mapitTo == to_wgs84.end()) || (mapitFrom == from_wgs84.end()))
		throw IOException("No known conversions exist for coordinate system " + coordsystem, AT);

	convToWGS84   = mapitTo->second;
	convFromWGS84 = mapitFrom->second;

	if(coordsystem=="LOCAL" && coordparam!="") {
		parseLatLon(coordparam, ref_latitude, ref_longitude);
	}
}

void Coords::clearCoordinates() {
//sets safe defaults for all internal variables (except function pointers and maps)
	latitude = longitude = IOUtils::nodata;
	altitude = IOUtils::nodata;
	easting = northing = IOUtils::nodata;
	grid_i = grid_j = grid_k = IOUtils::inodata;
}

void Coords::setDefaultValues() {
//sets safe defaults for all internal variables (except function pointers and maps)
	clearCoordinates();
	ref_latitude = ref_longitude = IOUtils::nodata;
	distance_algo = GEO_COSINE;
}

#ifdef _POPC_
#include "marshal_meteoio.h"
void Coords::Serialize(POPBuffer &buf, bool pack)
{
	if (pack){
		buf.Pack(&coordsystem, 1);
		buf.Pack(&coordparam, 1);
		buf.Pack(&ref_latitude, 1);
		buf.Pack(&ref_longitude, 1);
		buf.Pack(&latitude, 1);
		buf.Pack(&longitude, 1);
		buf.Pack(&altitude, 1);
		buf.Pack(&easting, 1);
		buf.Pack(&northing, 1);
		buf.Pack(&grid_i, 1);
		buf.Pack(&grid_j, 1);
		buf.Pack(&grid_k, 1);
		marshal_geo_distances(buf, distance_algo, 0, FLAG_MARSHAL, NULL);
	}else{
		buf.UnPack(&coordsystem, 1);
		buf.UnPack(&coordparam, 1);
		buf.UnPack(&ref_latitude, 1);
		buf.UnPack(&ref_longitude, 1);
		buf.UnPack(&latitude, 1);
		buf.UnPack(&longitude, 1);
		buf.UnPack(&altitude, 1);
		buf.UnPack(&easting, 1);
		buf.UnPack(&northing, 1);
		buf.UnPack(&grid_i, 1);
		buf.UnPack(&grid_j, 1);
		buf.UnPack(&grid_k, 1);
		marshal_geo_distances(buf, distance_algo, 0, !FLAG_MARSHAL, NULL);
		initializeMaps();
		setFunctionPointers();
	}
}
#endif
