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
#ifndef __CONFIGREADER_H__
#define __CONFIGREADER_H__

#include "IOUtils.h"
#include "IOExceptions.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <map>

/**
 * @class ConfigReader
 * @brief A class that reads a key/value file. These files (typically named *.ini) have the following structure:
 * - they consist of 0 or more explicitly indicated sections, which start with a sectionname enclosed in square brackets
 *   e.g. [General] or [Filter]
 * - within each section there are 0 or more key value pairs defined: KEY = VALUE 
 * - in this implementation each key is unique within its section
 * - lines that start with a semicolon ';' or a hash sign '#' are ignored (regarded as comments)
 * - empty lines are ignored
 * - if there is no section name given in a file, the default section called "GENERAL" is assumed 
 * - a VALUE for a KEY can consist of multiple whitespace separated values (e.g. MYNUMBERS = 17.77 -18.55 8888 99.99)
 *
 * @author Thomas Egger
 * @date   2008-11-30
 */
namespace mio {

namespace IOUtils {
	void toUpper(std::string& str);
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties,
								    const std::string& key, T& t);
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties, 
								    const std::string& key, std::vector<T>& vecT);
}

#ifdef _POPC_
class ConfigReader : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class ConfigReader {
#endif
	public:
		/**
		 * @brief Empty constructor. The user MUST later one fill the internal key/value map object
		*/
		ConfigReader();

		/**
		 * @brief Main constructor. The file is parsed and a key/value map object is internally created
		 * @param[in] filename_in string representing the absolute filename of the key/value file
		*/
		ConfigReader(const std::string& filename_in);

		/**
		 * @brief Add the content of a file to the internal key/value map object
		 * @param[in] filename_in string representing the absolute filename of the key/value file
		*/
		void addFile(const std::string& filename_in);

		/**
		 * @brief Add the content of the given command line to the internal key/value map object
		 * @param[in] cmd_line string representing the command line to be parsed for key/value pairs or switches
		*/
		void addCmdLine(const std::string& cmd_line);

		/**
		 * @brief Add a specific key/value pair to the internal key/value map object
		 * @param[in] key string representing the key to be added
		 * @param[in] value string representing the matching value to be added
		*/
		void addKey(const std::string& key, const std::string& value);

		/**
		 * @brief Add a specific key/value pair to the internal key/value map object
		 * @param[in] key string representing the key to be added
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[in] value string representing the matching value to be added
		*/
		void addKey(const std::string& key, const std::string& section, const std::string& value);

		/**
		 * @brief Returns the filename that the ConfigReader object was constructed with.
		 * @return std::string The absolute filename of the key/value file.
		 */
		std::string getSourceName();

		/**
		 * @brief Template function to retrieve a vector of values of class T for a certain key
		 * @code
		 * algorithms = lsm linregres idw_kriging\n
		 * MYNUMBERS = 17.87 19.89 -19.89 +7777.007
		 * @endcode
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[out] vecT a variable of class vector<T> into which the values for the corresponding key are saved
		 * @param[in] options indicating whether an exception should be raised, when key is not present
		 */
		template <class T> void getValue(const std::string& key, 
								   std::vector<T>& vecT, 
								   const unsigned int& options=0) const {
			getValue(key, "GENERAL", vecT, options);
		}
		
		/**
		 * @brief Template function to retrieve a vector of values of class T for a certain key
		 * @code
		 * algorithms = lsm linregres idw_kriging\n
		 * MYNUMBERS = 17.87 19.89 -19.89 +7777.007
		 * @endcode
		 * @param[in] key std::string representing a KEY in the key/value file
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[out] vecT a variable of class vector<T> into which the values for the corresponding key are saved 
		 * @param[in] options indicating whether an exception should be raised, when key is not present
		 */
		template <class T> void getValue(const std::string& key, 
								   const std::string& section, 
								   std::vector<T>& vecT, 
								   const unsigned int& options=0) const {
			try {
				vecT.clear();
				std::string _section = section;
				IOUtils::toUpper(_section);
				IOUtils::getValueForKey<T>(properties, _section + "::" + key, vecT);
			} catch(std::exception& e){
				if (options != ConfigReader::nothrow) {
					std::cerr << "[E] Error for ConfigReader of " << sourcename << ": " << e.what() << std::endl;
					throw;
				}
			}
		}

		/**
		 * @brief Template function to retrieve a value of class T for a certain key
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[out] t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		 * @param[in] options indicating whether an exception should be raised, when key is not present
		 */
		template <class T> void getValue(const std::string& key, T& t, const unsigned int& options=0) const {
			getValue(key, "GENERAL", t, options);
		}

		/**
		 * @brief Template function to retrieve a value of class T for a certain key
		 * @param[in] key std::string representing a KEY in the key/value file
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[out] t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		 * @param[in] options indicating whether an exception should be raised, when key is not present
		 */
		template <class T> void getValue(const std::string& key, const std::string& section, 
								   T& t, 
								   const unsigned int& options=0) const {
			try {
				std::string _section = section;
				IOUtils::toUpper(_section);
				IOUtils::getValueForKey<T>(properties, _section + "::" + key, t);
			} catch(std::exception& e){
				if (options != ConfigReader::nothrow) {
					std::cerr << "[E] Error for ConfigReader of " << sourcename << ": " << e.what() << std::endl;
					throw;
				}
			}
		}
		

		/**
		 * @brief Function that searches for a given string within the keys of section (default: GENERAL)
		 *         it returns the number of matches (partial matches are considered) and writes all the keys
		 *         into a vector\<string\> that is handed to the function as reference
		 * @param[out] vecResult A vector that will hold all keys that partially match keystart
		 * @param[in] keystart A string representing the beginning of a key to search for
		 * @param[in] section A string defining which section to search through (default: GENERAL)
		 * @code
		 *  vector<string> myVec;
		 *  unsigned int nrOfMatches = getKeys(myVec, "TA::", "Filters");
		 * @endcode
		 */
		unsigned int findKeys(std::vector<std::string>& vecResult, 
						  const std::string keystart, const std::string section="GENERAL") const;

		static const unsigned int nothrow;

	private:
		void parseCmdLine(const std::string& cmd_line);
		void parseFile(const std::string& filename);
		void parseLine(const unsigned int& linenr, std::string& line, std::string& section);

		std::map<std::string, std::string> properties; //Save key value pairs
		std::string sourcename; //description of the data source for the key/value pair

}; //end class definition ConfigReader

} //end namespace mio

#endif
