#include "MapProj.h"

void MapProj::initializeMaps()
{
	to_wgs84["CH1903"]   = &CH1903_to_WGS84;
	from_wgs84["CH1903"] = &WGS84_to_CH1903;
}


MapProj::MapProj(const std::string& _coordinatesystem, const std::string& _parameters)
{
	initializeMaps();
	
	//check whether there exists a tranformation for the given coordinatesystem
	//init function pointers
	std::map<std::string, convfunc>::iterator mapitTo;
	std::map<std::string, convfunc>::iterator mapitFrom;
	mapitTo   = to_wgs84.find(_coordinatesystem);	
	mapitFrom = from_wgs84.find(_coordinatesystem);	
	
	if ((mapitTo == to_wgs84.end()) || (mapitFrom == from_wgs84.end()))
		throw IOException("No known conversions exist for coordinate system " + _coordinatesystem, AT);
	
	convToWGS84   = mapitTo->second;
	convFromWGS84 = mapitFrom->second;

	coordsystem = _coordinatesystem;
	coordparam  = _parameters;
}

void MapProj::convert_to_WGS84(const double& easting, const double& northing, double& latitude, double& longitude)
{
	convToWGS84(easting, northing, latitude, longitude);
}

void MapProj::convert_from_WGS84(const double& latitude, const double& longitude, double& easting, double& northing)
{
	convFromWGS84(latitude, longitude, easting, northing);
}


void MapProj::WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out)
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

void MapProj::CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out)
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
