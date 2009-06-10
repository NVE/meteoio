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

#include "IOExceptions.h"
#include "Date_IO.h"

using namespace std;

namespace IOUtils {

	const double nodata = -999.0;	//HACK: we should define the same nodata everywhere...

	/**
	* @brief Coordinate conversion: from WGS84 Lat/Long to Swiss grid
	* @param lat_in Decimal Latitude 
	* @param long_in Decimal Longitude
	* @param east_out easting coordinate (Swiss system)
	* @param north_out northing coordinate (Swiss system)
	*/
	void WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out);

	/**
	* @brief Coordinate conversion: from Swiss grid to WGS84 Lat/Long
	* @param east_in easting coordinate (Swiss system)
	* @param north_in northing coordinate (Swiss system)
	* @param lat_out Decimal Latitude 
	* @param long_out Decimal Longitude
	*/
	void CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out);

	/**
	* @brief Check whether two values are equal regarding a certain epsilon environment (within certain radius of each other)
	* @param val1
	* @param val2
	* @param epsilon is a radius around val1
	* @return true if val2 is within the radius around val1, false otherwise.
	*/
	bool checkEpsilonEquality(double val1, double val2, double epsilon);

	void readDirectory(const string& path, list<string>& dirlist, const string& pattern = "");

	bool validFileName(const std::string& filename);

	bool fileExists(const std::string& filename);

	/**
	* Removes trailing and leading whitespaces, tabs and newlines from a string. 
	* @param string   The reference of the string to trim (in/out parameter)
	*/
	void trim(string &s);

	char getEoln(std::istream& fin);

	void skipLines(std::istream& fin, unsigned int nbLines, char eoln='\n');

	bool readKeyValuePair(const string& in_line, const string& delimiter, map<string,string>& out_map);

	unsigned int readLineToVec(const string& line_in, vector<string>& vecString);

	void readKeyValueHeader(map<string,string>& headermap, 
				    std::istream& bs,
				    const unsigned int& linecount=1, 
				    const std::string& delimiter="=");


	/** 
	* Convert a string to the requested type (template function). 
	* @param class T   [in] The type wanted for the return value (template type parameter). 
	* @param t   [out] The value converted to the requested type. 
	* @param str   [in] The input string to convert; trailling whitespaces are ignored, 
	*              comment after non-string values are allowed, but multiple values are not allowed. 
	* @param f   [in] The radix for reading numbers, such as std::dec or std::oct; default is std::dec. 
	* @return true if everything went fine, false otherwise
	*/
	template <class T> bool convertString(T& t, const std::string str, std::ios_base& (*f)(std::ios_base&) = std::dec) {
		string s = str; 
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
			string tmp="";
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
	template<> bool convertString<string>(string& t, const std::string str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<bool>(bool& t, const std::string str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<Date_IO>(Date_IO& t, const std::string str, std::ios_base& (*f)(std::ios_base&));

	/**
	* Returns, with the requested type, the value associated to a key (template function). 
	* @param class T   [in] The type wanted for the return value (template type parameter). 
	* @param properties   [in] A map containing all the parameters. 
	* @param key   [in] The key of the parameter to retrieve.
	* @param t   [out] The value associated to the key, converted to the requested type
	*/
	template <class T> void getValueForKey(const map<string,string>& properties, const string& key, T& t) {
		if (key == "") {
			THROW InvalidArgumentException("Empty key", AT);
		}
		const string value = (const_cast<map<string,string>&>(properties))[key];  

		if (value == "") {
			THROW UnknownValueException("No value for key " + key, AT);
		}

		if(!convertString<T>(t, value, std::dec)) {
			THROW ConversionFailedException(value, AT);
		}
	}
  
}; //end namespace IOUtils

#endif
