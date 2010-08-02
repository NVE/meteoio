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

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace mio {


/**
 * @class Config
 * @brief A class that reads a key/value file. These files (typically named *.ini) follow the INI file format standard (see http://en.wikipedia.org/wiki/INI_file) and have the following structure:
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

class ConfigProxy;

#ifdef _POPC_
class Config : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Config {
#endif
	public:
		enum Options { dothrow, nothrow };

		/**
		 * @brief Empty constructor. The user MUST later one fill the internal key/value map object
		*/
		Config();

		/**
		 * @brief Main constructor. The file is parsed and a key/value map object is internally created
		 * @param[in] filename_in string representing the absolute filename of the key/value file
		*/
		Config(const std::string& filename_in);

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
		 * @brief Add a specific key/value pair to the internal key/value map object.
		 *        key is case insensitive
		 * @param[in] key string representing the key to be added
		 * @param[in] value string representing the matching value to be added
		*/
		void addKey(const std::string& key, const std::string& value);

		/**
		 * @brief Add a specific key/value pair to the internal key/value map object. 
		 *        key and section are case insensitive
		 * @param[in] key string representing the key to be added
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[in] value string representing the matching value to be added
		*/
		void addKey(std::string key, std::string section, const std::string& value);

		/**
		 * @brief Returns the filename that the Config object was constructed with.
		 * @return std::string The absolute filename of the key/value file.
		 */
		std::string getSourceName();

		/**
		* @brief Print the content of the Config object (usefull for debugging)
		* The Config is bound by "<Config>" and "</Config>" on separate lines
		*/
		friend std::ostream& operator<<(std::ostream& os, const Config& cfg);

		template <typename T> std::vector<T> getValue(const std::string& key, const Options& opt=Config::dothrow) const
		{
			std::vector<T> tmp;
			getValue(key, Config::defaultSection, tmp, opt);
			return tmp;
		}

		template <typename T> std::vector<T> getValue(const std::string& key, const std::string& section,
                                                        const Options& opt=Config::dothrow) const
		{
			std::vector<T> tmp;
			getValue(key, section, tmp, opt);
			return tmp;
		}

		/**
		 * @brief Template function to retrieve a vector of values of class T for a certain key
		 * @code
		 * algorithms = lsm linregres idw_kriging\n
		 * MYNUMBERS = 17.87 19.89 -19.89 +7777.007
		 * @endcode
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[out] vecT a variable of class vector<T> into which the values for the corresponding key are saved
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key,
								   std::vector<T>& vecT,
								   const Options& opt=Config::dothrow) const
		{
			getValue(key, "GENERAL", vecT, opt);
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
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key, const std::string& section,
								      std::vector<T>& vecT, const Options& opt=Config::dothrow) const 
		{
			try {
				vecT.clear();
				std::string _key(key);
				std::string _section(section);
				IOUtils::toUpper(_key);
				IOUtils::toUpper(_section);
				IOUtils::getValueForKey<T>(properties, _section + "::" + _key, vecT);
			} catch(std::exception& e){
				if (opt != Config::nothrow) {
					std::stringstream ss;
					ss << "[E] Error for Config of " << sourcename << ": " << e.what();
					throw UnknownValueException(ss.str(), AT);
				}
			}
		}

		/**
		 * @ brief A function that allows to retrieve a value for a key as return parameter (vectors of values too)
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 * @return A value of type T
		 */
		ConfigProxy get(const std::string& key, const Options& opt=Config::dothrow) const;

		/**
		 * @ brief A function that allows to retrieve a value for a key as return parameter (vectors of values too)
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 * @return A value of type T
		 *
		 * Example Usage:
		 * @code
		 * Config cfg("io.ini");
		 * vector<int> = cfg.get("DEPTHS", "INPUT", Config::nothrow);
		 * string mystr = cfg.get("PATH", "OUTPUT");
		 * @endcode
		 */
		ConfigProxy get(const std::string& key, const std::string& section, const Options& opt=Config::dothrow) const;

		/**
		 * @brief Template function to retrieve a value of class T for a certain key
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[out] t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key, T& t, const Options& opt=Config::dothrow) const
		{
			getValue(key, "GENERAL", t, opt);
		}

		/**
		 * @brief Template function to retrieve a value of class T for a certain key
		 * @param[in] key std::string representing a KEY in the key/value file
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[out] t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key, const std::string& section, T& t,
                                              const Options& opt=Config::dothrow) const
		{
			try {
				std::string _key(key);
				std::string _section(section);
				IOUtils::toUpper(_key);
				IOUtils::toUpper(_section);
				IOUtils::getValueForKey<T>(properties, _section + "::" + _key, t);
			} catch(std::exception& e){
				if (opt != Config::nothrow) {
					std::stringstream ss;
					ss << "[E] Error for Config of " << sourcename << ": " << e.what();
					throw UnknownValueException(ss.str(), AT);
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
						  std::string keystart, std::string section="GENERAL") const;

		static const std::string defaultSection;

	private:
		void parseCmdLine(const std::string& cmd_line);
		void parseFile(const std::string& filename);
		void parseLine(const unsigned int& linenr, std::string& line, std::string& section);

		std::map<std::string, std::string> properties; //Save key value pairs
		std::string sourcename; //description of the data source for the key/value pair

}; //end class definition Config

class ConfigProxy {
 public:
	const Config& proxycfg;
	const std::string& key;
	const std::string& section;
     const Config::Options& opt;

	ConfigProxy(const Config& i_cfg, const std::string& i_key, 
			  const std::string& i_section, const Config::Options& i_opt) 
			: proxycfg(i_cfg), key(i_key),section(i_section), opt(i_opt) { }
			
	template<typename T> operator T()
		{ 
			T tmp;
			proxycfg.getValue(key, section, tmp, opt);
			return tmp; 
		}
};

} //end namespace mio

#endif
