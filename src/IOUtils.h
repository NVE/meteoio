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
#ifndef __IOUTILS_H__
#define __IOUTILS_H__


#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>
#include <cctype>
#include <limits>

#include "IOExceptions.h"
#include "Date_IO.h"
#include "Coords.h"
class Coords;

#ifndef MAX
#define MAX(x,y)    (((x) < (y)) ? (y) : (x))
#endif
#ifndef MIN
#define MIN(x,y)    (((x) < (y)) ? (x) : (y))
#endif

#ifndef C_TO_K
#define C_TO_K( T ) ( T + 273.15 )	  // Celsius to Kelvin
#endif
#ifndef K_TO_C
#define K_TO_C( T ) ( T - 273.15 )	  // Kelvin to Celsius
#endif

namespace IOUtils {
	const double nodata = -999.0; ///<This is the internal nodata value
	//const double not_set = std::numeric_limits<double>::max()-2.;
	const unsigned int unodata = (unsigned int)-1;
	const int inodata = -999;
	const unsigned int npos    = (unsigned int)-1;

	const double earth_radius = 6371e3; ///<Earth radius in meters
	const double grid_epsilon = 5.; ///<What is an acceptable small distance on a grid, in meters
	const double lon_epsilon = grid_epsilon / earth_radius; ///<in degrees. Small angle for longitudes, so sin(x)=x
	const double lat_epsilon = lon_epsilon/2.; ///<in degrees. Small angle for latitudes. Since for longitudes it is for 360deg, it has to be 180deg for latitudes

	typedef enum NODATA_HANLDING {
		RAW_NODATA, ///< no special handling of nodata
		PARSE_NODATA ///< process nodata as "no data"
	} nodata_handling;

	/**
	* @brief Check whether two values are equal regarding a certain epsilon environment (within certain radius of each other)
	* @param val1
	* @param val2
	* @param epsilon is a radius around val1
	* @return true if val2 is within the radius around val1, false otherwise.
	*/
	bool checkEpsilonEquality(const double& val1, const double& val2, const double& epsilon);

	double pow2(const double& val);
	double pow3(const double& val);
	double pow4(const double& val);

	void readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern = "");

	bool validFileName(const std::string& filename);

	bool fileExists(const std::string& filename);

	/**
	* @brief Removes trailing and leading whitespaces, tabs and newlines from a string. 
	* @param s The reference of the string to trim (in/out parameter)
	*/
	void trim(std::string &s);

	char getEoln(std::istream& fin);

	void skipLines(std::istream& fin, const unsigned int& nbLines, const char& eoln='\n');

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
	unsigned int readLineToVec(const std::string& line_in, std::vector<std::string>& vecString);
	unsigned int readLineToVec(const std::string& line_in, std::vector<std::string>& vecString, const char& delim);
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
	template <class T> bool convertString(T& t, const std::string& str, std::ios_base& (*f)(std::ios_base&) = std::dec) {
		std::string s = str; 
		trim(s); //delete trailing and leading whitespaces and tabs
		if (s.size() == 0) {
			t = static_cast<T> (nodata);
			return true;
		} else {
			std::istringstream iss(s);
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10); //try to read values with maximum precision
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
	template<> bool convertString<std::string>(std::string& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<bool>(bool& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<Date_IO>(Date_IO& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<Coords>(Coords& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));

	/**
	* @brief Returns, with the requested type, the value associated to a key (template function).
	* @tparam T   [in] The type wanted for the return value (template type parameter).
	* @param[in] properties   A map containing all the parameters.
	* @param[in] key   The key of the parameter to retrieve.
	* @param[out] t   The value associated to the key, converted to the requested type
	*/
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties,
						const std::string& key, T& t) {
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

	/**
	* @brief Returns, with the requested type, the value associated to a key (template function).
	* @tparam T           [in] The type wanted for the return value (template type parameter).
	* @param[in] properties  A map containing all the parameters.
	* @param[in] key         The key of the parameter to retrieve.
	* @param[out] vecT        The vector of values associated to the key, each value is converted to the requested type
	*/
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties, 
							const std::string& key, std::vector<T>& vecT) {
		if (key == "") {
			throw InvalidArgumentException("Empty key", AT);
		}
		const std::string value = (const_cast<std::map<std::string,std::string>&>(properties))[key];

		if (value == "") {
			throw UnknownValueException("No value for key " + key, AT);
		}

		//split value string 
		std::vector<std::string> vecUnconvertedValues;
		unsigned int counter = readLineToVec(value, vecUnconvertedValues);

		for (unsigned int ii=0; ii<counter; ii++){
			T myvar;
			if(!convertString<T>(myvar, vecUnconvertedValues.at(ii), std::dec)){
				throw ConversionFailedException(vecUnconvertedValues.at(ii), AT);
			}
			vecT.push_back(myvar);
		}
	}

	/**
	* @brief Standardize a given value to use MeteoIO's internal nodata value (if applicable)
	* @tparam T           [in] The type wanted for the return value (template type parameter).
	* @param[in] value  the value to check/convert
	* @param[in] plugin_nodata plugin-specific nodata value
	* @return checked/converted value
	*/
	template <class T> T standardizeNodata(const T& value, const double& plugin_nodata) {
		if(value==plugin_nodata) return static_cast<T> (nodata);
		else return value;
	}

} //end namespace IOUtils

#endif
