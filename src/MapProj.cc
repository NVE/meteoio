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
{
	to_wgs84["CH1903"]   = &MapProj::CH1903_to_WGS84;
	from_wgs84["CH1903"] = &MapProj::WGS84_to_CH1903;
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
		if(ref_latitude==IOUtils::nodata) {
			parseLocalParameters(ref_latitude, ref_longitude);
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

MapProj::MapProj(const std::string& _coordinatesystem, const double& _lat_ref, const double& _long_ref)
{
	initializeMaps();
	coordsystem = _coordinatesystem;
	if(coordsystem!="LOCAL") {
		throw InvalidArgumentException("Improperly using MapProj constructor that is ONLY suitable for LOCAL conversion", AT);
	}
	if(_lat_ref==IOUtils::nodata || _long_ref==IOUtils::nodata) {
		throw InvalidArgumentException("For LOCAL projection, please provide both reference latitude and longitude!", AT);
	}
	ref_latitude = _lat_ref;
	ref_longitude = _long_ref;
	coordparam = "";
	setFunctionPointers();
}

void MapProj::convert_to_WGS84(const double& easting, const double& northing, double& latitude, double& longitude) const
{
	(this->*convToWGS84)(easting, northing, latitude, longitude);
}

void MapProj::convert_from_WGS84(const double& latitude, const double& longitude, double& easting, double& northing) const
{
	(this->*convFromWGS84)(latitude, longitude, easting, northing);
}

void MapProj::parseLocalParameters(double& lat_ref, double& lon_ref) const
{//lat/lon for the reference point of local grid is read from projection parameter coordparam
	if 	((sscanf(coordparam.c_str(), "%lf %lf", &lat_ref, &lon_ref) < 2) && 
		(sscanf(coordparam.c_str(), "(%lf,%lf)", &lat_ref, &lon_ref) < 2) &&
		(sscanf(coordparam.c_str(), "%lf/%lf", &lat_ref, &lon_ref) < 2)) {
			throw InvalidFormatException("Can not parse given lat/lon: "+coordparam,AT);
	}
}

double MapProj::normalizeBearing(double angle)
{//TODO: use fmod
	if(angle<0.) angle = 360.0 + angle;
	if(angle>360.) angle = angle - 360.;
	return angle;
}


void MapProj::WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out) const
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

void MapProj::CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out) const
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

void MapProj::WGS84_to_PROJ4(const double& lat_in, const double& long_in, double& east_out, double& north_out) const
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

void MapProj::PROJ4_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out) const
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

void MapProj::WGS84_to_local(const double& lat_ref, const double& lon_ref, const double& lat, const double& lon, double& easting, double& northing, const enum GEO_DISTANCES algo)
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

void MapProj::local_to_WGS84(const double& lat_ref, const double& lon_ref, const double& easting, const double& northing, double& lat, double& lon, const enum GEO_DISTANCES algo)
{
	const double to_deg = 180.0 / PI;
	double bearing = atan2(northing, easting)*to_deg;
	const double distance = sqrt( IOUtils::pow2(easting) + IOUtils::pow2(northing) );

	bearing = normalizeBearing(bearing); //TODO: are we really using compass bearing here?

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

void MapProj::WGS84_to_local(const double& lat_in, const double& long_in, double& east_out, double& north_out) const
{
	WGS84_to_local(ref_latitude, ref_longitude, lat_in, long_in, east_out, north_out, GEO_COSINE);
}

void MapProj::local_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out) const
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
	double alpha2 = atan2(cos(U1)*sin(lambda), sin(U1)*cos(U2)-cos(U1)*sin(U2)*cos(lambda)) / to_rad; //reverse azimuth

	//trying to get a normal compass bearing... TODO: make sure it works and understand why
	alpha1 = fmod(-alpha1+360., 360.);
	alpha2 = fmod(alpha2+180., 360.);

	alpha = (alpha1+alpha2)/2.;
	return s;
}

void MapProj::VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon)
{
	const double thresh = 1.e-12;	//convergence absolute threshold
	const double a = 6378137.;	//major ellipsoid semi-axis, value for wgs84
	const double b = 6356752.3142;	//minor ellipsoid semi-axis, value for wgs84
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = PI / 180.0;
	const double to_deg = 180.0 / PI;

	const double alpha1 = bearing*to_rad;
	const double tanU1 = (1.-f)*tan(lat_ref*to_rad);
	const double cosU1 = 1./sqrt(1.+IOUtils::pow2(tanU1));
	const double sinU1 = tanU1*cosU1;
	const double sigma1 = atan2(tanU1,cos(alpha1));
	const double sinAlpha = cosU1*sin(alpha1);
	const double cos2alpha = 1. - IOUtils::pow2(sinAlpha);
	const double u2 = cos2alpha * (a*a - b*b) / (b*b);
	const double A = 1. + u2/16384. * (4096. + u2*(-768.+u2*(320.-175.*u2)) );
	const double B = u2/1024. * (256. + u2*(-128.+u2*(74.-47.*u2)));
	
	double sigma = distance / (b*A);
	double sigma_p = 2.*PI;
	double cos2sigma_m = cos( 2.*sigma1 + sigma ); //required to avoid uninitialized value

	while (fabs(sigma - sigma_p) > thresh) {
		cos2sigma_m = cos( 2.*sigma1 + sigma );
		double delta_sigma = B*sin(sigma) * ( cos2sigma_m + B/4. * ( 
			cos(sigma)*(-1.+2.*IOUtils::pow2(cos2sigma_m)) 
			-B/6. * cos2sigma_m * (-3.+4.*IOUtils::pow2(sin(sigma))) * (-3.+4.*IOUtils::pow2(cos2sigma_m))
			) );
		sigma_p = sigma;
		sigma = distance / (b*A) + delta_sigma;
	}
	
	lat = atan2( sinU1*cos(sigma) + cosU1*cos(sigma)*cos(alpha1),
		     (1.-f) * sqrt( IOUtils::pow2(sinAlpha) + IOUtils::pow2(sinU1*sin(sigma) - cosU1*cos(sigma)*cos(alpha1)) )
		   );
	const double lambda = atan2( sin(sigma)*sin(alpha1), cosU1*cos(sigma) - sinU1*sin(sigma)*cos(alpha1) );
	const double C = f/16. * cos2alpha * (4.+f*(4.-3.*cos2alpha));
	const double L = lambda - (1.-C) * f * sinAlpha * (
				sigma + C * sin(sigma) * ( cos2sigma_m+C*cos(sigma) * (-1.+2.*IOUtils::pow2(cos2sigma_m)) )
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
