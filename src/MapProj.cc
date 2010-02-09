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
#include "MapProj.h"

#ifdef PROJ4
	#include <proj_api.h>
#endif
#ifndef PI
	#define PI 3.141592653589
#endif

const struct MapProj::ELLIPSOID MapProj::ellipsoids[] = {
		{ 6378137.,	6356752.3142 }, //E_WGS84
		{ 6378137.,	6356752.3141 }, //E_GRS80
		{ 6377563.396,	6356256.909 }, //E_AIRY
		{ 6378388.,	6356911.946 }, //E_INTL1924
		{ 6378249.145,	6356514.86955 }, //E_CLARKE1880
		{ 6378160.,	6356774.719 } //E_GRS67
};

void MapProj::initializeMaps()
{	//Please don't forget to mirror the keywords here in the documentation in MapProj.h!!
	to_wgs84["CH1903"]   = &MapProj::CH1903_to_WGS84;
	from_wgs84["CH1903"] = &MapProj::WGS84_to_CH1903;
	to_wgs84["UTM"]   = &MapProj::UTM_to_WGS84;
	from_wgs84["UTM"] = &MapProj::WGS84_to_UTM;
	to_wgs84["PROJ4"]   = &MapProj::PROJ4_to_WGS84;
	from_wgs84["PROJ4"] = &MapProj::WGS84_to_PROJ4;
	to_wgs84["LOCAL"]   = &MapProj::local_to_WGS84;
	from_wgs84["LOCAL"] = &MapProj::WGS84_to_local;
}

void MapProj::setFunctionPointers()
{
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

	if(coordsystem=="LOCAL") {
		//if ref_latitude has not yet been set (ie: we have not come from the specific
		//local constructor
		if(ref_latitude==IOUtils::nodata) {
			parseLatLon(coordparam, ref_latitude, ref_longitude);
		}
	}
}

MapProj::MapProj()
{
	coordsystem = "NULL";
	coordparam  = "NULL";
}

MapProj::MapProj(const std::string& _coordinatesystem, const std::string& _parameters)
{
	initializeMaps();
	coordsystem = _coordinatesystem;
	coordparam  = _parameters;
	ref_latitude = ref_longitude = IOUtils::nodata;
	setFunctionPointers();
}

MapProj::MapProj(const double& _lat_ref, const double& _long_ref)
{
	initializeMaps();
	coordsystem = std::string("LOCAL");
	
	if(_lat_ref==IOUtils::nodata || _long_ref==IOUtils::nodata) {
		throw InvalidArgumentException("For LOCAL projection, please provide both reference latitude and longitude!", AT);
	}
	ref_latitude = _lat_ref;
	ref_longitude = _long_ref;
	coordparam = "";
	setFunctionPointers();
}

void MapProj::convert_to_WGS84(double easting, double northing, double& latitude, double& longitude) const
{
	(this->*convToWGS84)(easting, northing, latitude, longitude);
}

void MapProj::convert_from_WGS84(double latitude, double longitude, double& easting, double& northing) const
{
	(this->*convFromWGS84)(latitude, longitude, easting, northing);
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
double MapProj::dms_to_decimal(const std::string& dms) {
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
void MapProj::parseLatLon(const std::string& coordinates, double&
lat, double& lon) {
	char lat_str[32]="";
	char lon_str[32]="";

	if(coordinates.size()>(32+32)) {
		throw InvalidFormatException("Given lat/lon string is too long! ",AT);
	}

	if 	((sscanf(coordinates.c_str(), "%[0-9.,°d'\"] %[0-9.,°d'\"]", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "%[0-9.,°d'\" ]/%[0-9.,°d'\" ]", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "(%[0-9.,°d'\" ];%[0-9.,°d'\" ])", lat_str, lon_str) < 2) &&
		(sscanf(coordinates.c_str(), "(%[0-9.°d'\" ],%[0-9.°d'\" ])", lat_str, lon_str) < 2)) {
			throw InvalidFormatException("Can not parse given lat/lon: "+coordinates,AT);
	}

	lat = dms_to_decimal(string(lat_str));
	lon = dms_to_decimal(string(lon_str));
}

/**
* @brief Converts a decimal latitude or longitude to degrees, minutes, seconds
* It formats its arguments as in the following example: 46°48'03"
* @param[in] decimal decimal coordinate to convert
* @return string containing the formatted coordinate
*/
std::string MapProj::decimal_to_dms(const double& decimal) {
	std::stringstream dms;
	int d = (int)floor(decimal);
	int m = (int)floor( (decimal - (double)d)*60. );
	int s = (int)floor( decimal - (double)d - (double)m/60. );

	dms << d << "°" << m << "'" << s << "\"";
	return dms.str();
}


void MapProj::WGS84_to_CH1903(double lat_in, double long_in, double& east_out, double& north_out) const
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

void MapProj::CH1903_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
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

int MapProj::getUTMZone(const double latitude, const double longitude, std::string& zone_out) const
{//This routine determines the correct UTM letter designator for the given latitude
//UTM limits its coverage to [80S , 84N], outside of this, returns Z for the zone

	//computing zone number, assuming longitude in [-180. ; 180[
	int ZoneNumber = int((longitude + 180.)/6.) + 1;

	// Special zones for Scandinavia
	if( latitude >= 72.0 && latitude < 84.0 ) {
		if(      longitude >= 0.0  && longitude <  9.0 ) ZoneNumber = 31;
		else if( longitude >= 9.0  && longitude < 21.0 ) ZoneNumber = 33;
		else if( longitude >= 21.0 && longitude < 33.0 ) ZoneNumber = 35;
		else if( longitude >= 33.0 && longitude < 42.0 ) ZoneNumber = 37;
	 }
	if( latitude >= 56.0 && latitude < 64.0 && longitude >= 3.0 && longitude < 12.0 ) {
		ZoneNumber = 32;
	}

	//getting zone letter
	char zoneLetter = 'Z';
	if     ((84 >= latitude) && (latitude >= 72)) zoneLetter = 'X';
	else if((72 > latitude) && (latitude >= 64)) zoneLetter = 'W';
	else if((64 > latitude) && (latitude >= 56)) zoneLetter = 'V';
	else if((56 > latitude) && (latitude >= 48)) zoneLetter = 'U';
	else if((48 > latitude) && (latitude >= 40)) zoneLetter = 'T';
	else if((40 > latitude) && (latitude >= 32)) zoneLetter = 'S';
	else if((32 > latitude) && (latitude >= 24)) zoneLetter = 'R';
	else if((24 > latitude) && (latitude >= 16)) zoneLetter = 'Q';
	else if((16 > latitude) && (latitude >= 8)) zoneLetter = 'P';
	else if(( 8 > latitude) && (latitude >= 0)) zoneLetter = 'N';
	else if(( 0 > latitude) && (latitude >= -8)) zoneLetter = 'M';
	else if((-8 > latitude) && (latitude >= -16)) zoneLetter = 'L';
	else if((-16 > latitude) && (latitude >= -24)) zoneLetter = 'K';
	else if((-24 > latitude) && (latitude >= -32)) zoneLetter = 'J';
	else if((-32 > latitude) && (latitude >= -40)) zoneLetter = 'H';
	else if((-40 > latitude) && (latitude >= -48)) zoneLetter = 'G';
	else if((-48 > latitude) && (latitude >= -56)) zoneLetter = 'F';
	else if((-56 > latitude) && (latitude >= -64)) zoneLetter = 'E';
	else if((-64 > latitude) && (latitude >= -72)) zoneLetter = 'D';
	else if((-72 > latitude) && (latitude >= -80)) zoneLetter = 'C';

	std::stringstream zone;
	zone << ZoneNumber << zoneLetter;
	zone_out = zone.str();

	return ZoneNumber;
}

void MapProj::WGS84_to_UTM(double lat_in, double long_in, double& east_out, double& north_out) const
{//converts WGS84 coordinates (lat,long) to UTM coordinates.
//See USGS Bulletin 1532 or http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf
//also http://www.uwgb.edu/dutchs/usefuldata/UTMFormulas.HTM
//also http://www.oc.nps.edu/oc2902w/maps/utmups.pdf or Chuck Gantz (http://www.gpsy.com/gpsinfo/geotoutm/)
	//Geometric constants
	const double to_rad = PI / 180.0;
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
	const double nu = a / sqrt(1.-e2*pow2(sin(Lat))); //radius of curvature of the earth perpendicular to the meridian plane
	const double p = (Long-long0);

	//calculating first the coefficients of the series, then the Meridional Arc M itself
	const double n = (a-b)/(a+b);
	const double n2=n*n, n3=n*n*n, n4=n*n*n*n, n5=n*n*n*n*n, n6=n*n*n*n*n*n;
	const double A = a           * (1. - n + 5./4.*(n2 - n3) + 81./64.*(n4 - n5));
	const double B = (3./2.*a)   * (n - n2 + 7./8.*(n3 - n4) + 55./64.*(n5 - n6));
	const double C = (15./16.*a) * (n2 - n3 + 3./4.*(n4 - n5));
	const double D = (35./48.*a) * (n3 - n4 + 11./16.*(n5 - n6));
	//const double E = (315./512.*a) * (n4 - n5); //correctrion of ~0.03mm
	const double M = A*Lat - B*sin(2.*Lat) + C*sin(4.*Lat) - D*sin(6.*Lat) /*+ E*sin(8.*Lat)*/;

	//calculating the coefficients for the series
	const double K1 = M*k0;
	const double K2 = 1./4.*k0*nu*sin(2.*Lat);
	const double K3 = (k0*nu*sin(Lat)*pow3(cos(Lat))*1./24.) * (5. - pow2(tan(Lat)) + 9.*eP2*pow2(cos(Lat)) + 4.*eP2*eP2*pow4(cos(Lat)));
	const double K4 = k0*nu*cos(Lat);
	const double K5 = (k0*nu*pow3(cos(Lat))*1./6.) * (1. - pow2(tan(Lat)) + eP2*pow2(cos(Lat)));

	north_out = K1 + K2*p*p + K3*p*p*p*p;
	east_out = K4*p + K5*p*p*p + 500000.0;

	if(Lat < 0)
		north_out += 10000000.0; //offset for southern hemisphere
}

void MapProj::UTM_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{//converts UTM coordinates to WGS84 coordinates (lat,long).
//See USGS Bulletin 1532 or http://earth-info.nga.mil/GandG/publications/tm8358.2/TM8358_2.pdf
//also http://www.uwgb.edu/dutchs/usefuldata/UTMFormulas.HTM
//also http://www.oc.nps.edu/oc2902w/maps/utmups.pdf or Chuck Gantz (http://www.gpsy.com/gpsinfo/geotoutm/)
	//Geometric constants
	const double to_deg = 180.0 / PI;
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double e2 = (a*a-b*b) / (a*a);	//ellispoid eccentricity, squared
	const double eP2 = e2 / (1.-e2);	//second ellispoid eccentricity, squared (=(a²-b²)/b²)
	const double k0 = 0.9996;		//scale factor for the projection

	//UTM Zone information
	char zoneLetter;
	int zoneNumber;
	if 	((sscanf(coordparam.c_str(), "%d%c", &zoneNumber, &zoneLetter) < 2) &&
		(sscanf(coordparam.c_str(), "%d %c)", &zoneNumber, &zoneLetter) < 2)) {
			throw InvalidFormatException("Can not parse given UTM zone: "+coordparam,AT);
	}
	if(zoneLetter<='N') {
		north_in -= 10000000.0; //offset used for southern hemisphere
	}
	east_in -= 500000.0; //longitude offset: x coordinate is relative to central meridian
	const int long0 = (zoneNumber - 1)*6 - 180 + 3;  //+3 puts origin in middle of zone; HACK: this does not account for the Scandinavian irregular zones

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
	const double C1 = eP2 * pow2(cos(fp));
	const double T1 = pow2( tan(fp) );
	const double R1 = a*(1.-e2) / pow((1.-e2*pow2(sin(fp))), 1.5);
	const double N1 = a / sqrt(1.-e2*pow2(sin(fp)));
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
	long_out = long0 + ((Q5 - Q6 + Q7)/cos(fp))*to_deg;
}

void MapProj::WGS84_to_PROJ4(double lat_in, double long_in, double& east_out, double& north_out) const
{
#ifdef PROJ4
	const string src_param="+proj=latlong +datum=WGS84 +ellps=WGS84";
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

void MapProj::PROJ4_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
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

void MapProj::WGS84_to_local(double lat_ref, double lon_ref, const double& lat, const double& lon, double& easting, double& northing, const enum GEO_DISTANCES algo)
{
	double alpha;
	const double to_rad = PI / 180.0;
	double distance;

	switch(algo) {
		case GEO_COSINE:
			distance = cosineDistance(lat_ref, lon_ref, lat, lon, alpha);
			break;
		case GEO_VINCENTY:
			distance = VincentyDistance(lat_ref, lon_ref, lat, lon, alpha);
			break;
		default:
			throw InvalidArgumentException("Unrecognized geodesic distance algorithm selected", AT);
	}
	
	easting = distance*sin(alpha*to_rad);
	northing = distance*cos(alpha*to_rad);
}

void MapProj::local_to_WGS84(double lat_ref, double lon_ref, const double& easting, const double& northing, double& lat, double& lon, const enum GEO_DISTANCES algo)
{
	const double to_deg = 180.0 / PI;
	const double distance = sqrt( IOUtils::pow2(easting) + IOUtils::pow2(northing) );
	const double bearing = fmod( atan2(easting, northing)*to_deg+360. , 360.);

	switch(algo) {
		case GEO_COSINE:
			cosineInverse(lat_ref, lon_ref, distance, bearing, lat, lon);
			break;
		case GEO_VINCENTY:
			VincentyInverse(lat_ref, lon_ref, distance, bearing, lat, lon);
			break;
		default:
			throw InvalidArgumentException("Unrecognized geodesic distance algorithm selected", AT);
	}
}

void MapProj::WGS84_to_local(double lat_in, double long_in, double& east_out, double& north_out) const
{
	WGS84_to_local(ref_latitude, ref_longitude, lat_in, long_in, east_out, north_out, GEO_COSINE);
}

void MapProj::local_to_WGS84(double east_in, double north_in, double& lat_out, double& long_out) const
{
	local_to_WGS84(ref_latitude, ref_longitude, east_in, north_in, lat_out, long_out, GEO_COSINE);
}

void MapProj::cosineInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon)
{
	const double Rearth = 6371.e3;
	const double to_rad = PI / 180.0;
	const double to_deg = 180.0 / PI;
	const double lat_ref_rad = lat_ref*to_rad;
	const double bearing_rad = bearing*to_rad;

	lat = asin( sin(lat_ref_rad)*cos(distance/Rearth) + 
				cos(lat_ref_rad)*sin(distance/Rearth)*cos(bearing_rad) );
	lon = lon_ref*to_rad + atan2( sin(bearing_rad)*sin(distance/Rearth)*cos(lat_ref_rad) , 
					cos(distance/Rearth) - sin(lat_ref_rad)*sin(lat) );
	lon = fmod(lon+PI, 2.*PI) - PI;

	lat *= to_deg;
	lon *= to_deg;
}

double MapProj::cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2)
{
	double alpha;
	return cosineDistance(lat1, lon1, lat2, lon2, alpha);
}

double MapProj::cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha)
{
	const double Rearth = 6371.e3;
	const double to_rad = PI / 180.0;
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

double MapProj::VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2)
{
	double alpha;
	return VincentyDistance(lat1, lon1, lat2, lon2, alpha);
}

double MapProj::VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha)
{
	const double thresh = 1.e-12;	//convergence absolute threshold
	const int n_max = 100;		//maximum number of iterations
	const double a = ellipsoids[E_WGS84].a; //major ellipsoid semi-axis
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = PI / 180.0;
	
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

void MapProj::VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon)
{//well, actually this is the DIRECT Vincenty formula
	const double thresh = 1.e-12;	//convergence absolute threshold
	const double a = ellipsoids[E_WGS84].a;	//major ellipsoid semi-axis, value for wgs84
	const double b = ellipsoids[E_WGS84].b;	//minor ellipsoid semi-axis, value for wgs84
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = PI / 180.0;
	const double to_deg = 180.0 / PI;

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
	double sigma_p = 2.*PI;
	double cos2sigma_m = cos( 2.*sigma1 + sigma ); //required to avoid uninitialized value

	while (fabs(sigma - sigma_p) > thresh) {
		cos2sigma_m = cos( 2.*sigma1 + sigma );
		double delta_sigma = B*sin(sigma) * ( cos2sigma_m + B/4. * ( 
			cos(sigma)*(-1.+2.*cos2sigma_m*cos2sigma_m) 
			-B/6. * cos2sigma_m * (-3.+4.*pow2(sin(sigma))) * (-3.+4.*cos2sigma_m*cos2sigma_m)
			) );
		sigma_p = sigma;
		sigma = distance / (b*A) + delta_sigma;
	}
	
	lat = atan2( sinU1*cos(sigma) + cosU1*sin(sigma)*cos(alpha1),
		     (1.-f) * sqrt( sinAlpha*sinAlpha + pow2(sinU1*sin(sigma) - cosU1*cos(sigma)*cos(alpha1)) )
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

bool MapProj::operator==(const MapProj& in) const
{//we only compare on coordsystem/coordparam since the function pointers might be different 
 //when comparing between different objects. Obviously, this assumes that the function pointers are
 //always correctly initialized...

	return ((coordsystem==in.coordsystem) && (coordparam==in.coordparam));
}

bool MapProj::operator!=(const MapProj& in) const
{
	return !(*this==in);
}

#ifdef _POPC_
void MapProj::Serialize(POPBuffer &buf, bool pack)
{
	if (pack){
		buf.Pack(&coordsystem, 1);
		buf.Pack(&coordparam, 1);
		buf.Pack(&ref_latitude, 1);
		buf.Pack(&ref_longitude, 1);
	}else{
		buf.UnPack(&coordsystem, 1);
		buf.UnPack(&coordparam, 1);
		buf.UnPack(&ref_latitude, 1);
		buf.UnPack(&ref_longitude, 1);
		initializeMaps();
		setFunctionPointers();
	}
}
#endif
