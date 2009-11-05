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
#ifdef _POPC_
class ConfigReader : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
		ConfigReader(){};
#else
class ConfigReader {
#endif
	public:
		/**
		 * @brief Main constructor. The file is parsed and a key/value map object is internally created
		 * @param[in] filename_in string representing the absolute filename of the key/value file
		*/
		ConfigReader(const std::string& filename_in);

		/**
		 * @brief Returns the filename that the ConfigReader object was constructed with.
		 * @return std::string The absolute filename of the key/value file.
		 */
		std::string getFileName();

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
								   vector<T>& vecT, 
								   const unsigned int& options=0) const {
			try {
				vecT.clear();
				IOUtils::getValueForKey<T>(properties, section + "::" + key, vecT);
			} catch(std::exception& e){
				if (options != ConfigReader::nothrow)
					throw;
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
				IOUtils::getValueForKey<T>(properties, section + "::" + key, t);
			} catch(std::exception& e){
				if (options != ConfigReader::nothrow)
					throw;
			}
		}
		

		/**
		 * @brief Function that searches for a given string within the keys of section (default: GENERAL)
		 *         it returns the number of matches (partial matches are considered) and writes all the keys
		 *         into a vector<string> that is handed to the function as reference
		 * @param[out] vecResult A vector that will hold all keys that partially match keystart
		 * @param[in] keystart A string representing the beginning of a key to search for
		 * @param[in] section A string defining which section to search through (default: GENERAL)
		 * @code
		 *  vector<string> myVec;
		 *  unsigned int nrOfMatches = getKeys(myVec, "TA::", "Filters");
		 * @endcode
		 */
		unsigned int findKeys(vector<string>& vecResult, const std::string keystart, const std::string section="GENERAL");

		//LEGACY
		static bool readConfigLine(std::istream& fin, int lineNb, int& lineType, string& str1, string& str2);

		/** [Constant for lineType] The config line is the end of the file */
		static const int CfgLineEOF = 0;
		/** [Constant for lineType] The config line is a comment (or a blank line) */
		static const int CfgLineComment = 1;
		/** [Constant for lineType] The config line is the start of a section */
		static const int CfgLineSection = 2;
		/** [Constant for lineType] The config line is a key-value pair */
		static const int CfgLineKeyValue = 3;
		/** [Constant for lineType] The config line is of an unknown type (in fact, an invalid format) */
		static const int CfgLineUnknown = 9;

		static const unsigned int nothrow;

	private:
		void parseFile();
		void parseLine(const unsigned int& linenr, std::string& line, std::string& section);

		std::map<string, string> properties; //Save key value pairs
		std::string filename; //Absolute filename of the key/value file

}; //end class definition ConfigReader

#endif
