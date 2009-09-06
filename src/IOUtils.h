#ifndef __IOUTILS_H__
#define __IOUTILS_H__

#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <cctype>

#include "IOExceptions.h"
#include "Date_IO.h"

#ifndef MAX
#define MAX(x,y)    (((x) < (y)) ? (y) : (x))
#endif
#ifndef MIN
#define MIN(x,y)    (((x) < (y)) ? (x) : (y))
#endif

namespace IOUtils {

	const double nodata = -999.0;	//HACK: we should define the same nodata everywhere...
	const unsigned int unodata = (unsigned int)-1;
	const unsigned int npos    = (unsigned int)-1;

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
	* @brief Coordinate conversion: from WGS84 Lat/Long to local metric grid
	* @param lat_ref Decimal Latitude of origin (const double&)
	* @param lon_ref Decimal Longitude of origin (const double&)
	* @param lat Decimal Latitude (const double&)
	* @param lon Decimal Longitude(const double&)
	* @param easting easting coordinate in local grid (double&)
	* @param northing northing coordinate in local grid (double&)
	*/
	void WGS84_to_local(const double& lat_ref, const double& lon_ref, const double& lat, const double& lon, double& easting, double& northing);

	/**
	* @brief Coordinate conversion: from local metric grid to WGS84 Lat/Long
	* @param lat_ref Decimal Latitude of origin (const double&)
	* @param lon_ref Decimal Longitude of origin (const double&)
	* @param easting easting coordinate in local grid (double&)
	* @param northing northing coordinate in local grid (double&)
	* @param lat Decimal Latitude (double&)
	* @param lon Decimal Longitude(double&)
	*/
	void local_to_WGS84(const double& lat_ref, const double& lon_ref, const double& easting, const double& northing, double& lat, double& lon);

	/**
	* @brief Spherical law of cosine Distance calculation between points in WGS84 (decimal Lat/Long)
	* See http://www.movable-type.co.uk/scripts/latlong.html for more
	* @param lat1 Decimal Latitude (const double&)
	* @param lon1 Decimal Longitude (const double&)
	* @param lat2 Decimal Latitude (const double&)
	* @param lon2 Decimal Longitude (const double&)
	* @return distance (double)
	*/
	double cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2);

	/**
	* @brief Vincenty Distance calculation between points in WGS84 (decimal Lat/Long)
	* See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems", 
	* Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599, 
	* see http://www.springerlink.com/content/y7108u6862473583 for more
	* @param lat1 Decimal Latitude (const double&)
	* @param lon1 Decimal Longitude (const double&)
	* @param lat2 Decimal Latitude (const double&)
	* @param lon2 Decimal Longitude (const double&)
	* @param alpha average bearing (double&)
	* @return distance (double)
	*/
	double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha);
	double VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2);

	/**
	* @brief Vincenty Inverse calculation giving WGS84 (decimal Lat/Long) position
	* given a start location (lat,lon) a distance and a bearing
	* See T. Vincenty, "Closed formulas for the direct and reverse geodetic problems", 
	* Journal of Geodesy, 51, 3, 1977, DOI:10.1007/BF02521599, 
	* see http://www.springerlink.com/content/y7108u6862473583 for more
	* @param lat_ref Decimal Latitude (const double&)
	* @param lon_ref Decimal Longitude (const double&)
	* @param distance Distance in meters (const double&)
	* @param bearing bearing in degrees, 0 being north (const double&)
	* @param lat Decimal latitude of target point (double&)
	* @param lon Decimal longitude of target point (double&)
	*/
	void VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon);

	/**
	* @brief Check whether two values are equal regarding a certain epsilon environment (within certain radius of each other)
	* @param val1
	* @param val2
	* @param epsilon is a radius around val1
	* @return true if val2 is within the radius around val1, false otherwise.
	*/
	bool checkEpsilonEquality(double val1, double val2, double epsilon);

	double pow2(const double val);

	/**
	* @brief Converts an angle to a bearing, that is between 0 and 360 degrees
	* @param angle angle in degrees (double&)
	* @return bearing between 0 and 360 degrees (double)
	*/
	double normalizeBearing(double angle);

	void readDirectory(const string& path, list<string>& dirlist, const string& pattern = "");

	bool validFileName(const std::string& filename);

	bool fileExists(const std::string& filename);

	/**
	* @brief Removes trailing and leading whitespaces, tabs and newlines from a string. 
	* @param s The reference of the string to trim (in/out parameter)
	*/
	void trim(std::string &s);

	char getEoln(std::istream& fin);

	void skipLines(std::istream& fin, unsigned int nbLines, char eoln='\n');

	/**
	* @brief read a string line, parse it and save it into a map object, that is passed by reference
	* @param in_line (const string&) string to parse
	* @param delimiter (const string&) delimiter to use for the parsing
	* @param out_map (map\<string,string\>&) map after parsing
	* @param keyprefix this string is prefixed before the key, defaults to no prefix: ""
	* @return (bool) true when line is empty
	*/
	bool readKeyValuePair(const std::string& in_line, 
					  const std::string& delimiter, 
					  std::map<std::string,std::string>& out_map,
					  const std::string& keyprefix="");

	void toUpper(std::string& str);
	unsigned int readLineToVec(const std::string& line_in, std::vector<string>& vecString);
	unsigned int readLineToVec(const std::string& line_in, std::vector<string>& vecString, const char& delim);
	void readKeyValueHeader(std::map<std::string, std::string>& headermap, 
				    std::istream& bs,
				    const unsigned int& linecount=1, 
				    const std::string& delimiter="=");


	/** 
	* @brief Convert a string to the requested type (template function). 
	* @tparam T   [in] The type wanted for the return value (template type parameter). 
	* @param t   [out] The value converted to the requested type. 
	* @param str   [in] The input string to convert; trailling whitespaces are ignored, 
	*              comment after non-string values are allowed, but multiple values are not allowed. 
	* @param f   [in] The radix for reading numbers, such as std::dec or std::oct; default is std::dec. 
	* @return true if everything went fine, false otherwise
	*/
	template <class T> bool convertString(T& t, const std::string str, std::ios_base& (*f)(std::ios_base&) = std::dec) {
		std::string s = str; 
		trim(s); //delete trailing and leading whitespaces and tabs
		if (s.size() == 0) {
			t = static_cast<T> (nodata);
			return true;
		} else {
			std::istringstream iss(s);
			iss >> f >> t; //Convert first part of stream with the formatter (e.g. std::dec, std::oct)
		
			if (iss.fail()) {
				//Conversion failed
				return false;
			}
			std::string tmp="";
			getline(iss,  tmp); //get rest of line, if any
			trim(tmp);
			if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {
				//if line holds more than one value it's invalid
				return false;
			}
			return true;
		}
	}
	// fully specialized template functions (implementation must not be in header)
	template<> bool convertString<std::string>(std::string& t, const std::string str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<bool>(bool& t, const std::string str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<Date_IO>(Date_IO& t, const std::string str, std::ios_base& (*f)(std::ios_base&));

	/**
	* @brief Returns, with the requested type, the value associated to a key (template function).
	* @tparam T   [in] The type wanted for the return value (template type parameter).
	* @param properties   [in] A map containing all the parameters.
	* @param key   [in] The key of the parameter to retrieve.
	* @param t   [out] The value associated to the key, converted to the requested type
	*/
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties, const std::string& key, T& t) {
		if (key == "") {
			throw InvalidArgumentException("Empty key", AT);
		}
		const std::string value = (const_cast<std::map<std::string,std::string>&>(properties))[key];

		if (value == "") {
			throw UnknownValueException("No value for key " + key, AT);
		}

		if(!convertString<T>(t, value, std::dec)) {
			throw ConversionFailedException(value, AT);
		}
	}
  
}; //end namespace IOUtils

#endif
