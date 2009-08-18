#ifndef __CONFIGREADER_H__
#define __CONFIGREADER_H__

#include "IOUtils.h"
#include "IOExceptions.h"

#include <fstream>
#include <string>
#include <map>

/**
 * @class ConfigReader
 * @brief A class that reads a key/value file. These files (typically named *.ini) have the following structure:
 * - line 1 contains the string "[Parameters]" (without the quotation marks)
 * - all following lines are either blank, comments or contain a "KEY = VALUE" assignment
 * - in this implementation each key is unique
 * - lines that don't have a delimiter in them are ignored
 * - lines that start with a semicolon ';' or a hash sign '#' are ignored (regarded as comments)
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
		* @param filename_in string representing the absolute filename of the key/value file
		*/
		ConfigReader(const std::string& filename_in);

		/**
		* @brief Explicit copy constructor. It is required, because the private ifstream fin cannot be copied.
		* @param configreader ConfigReader object to be copied
		*/
		ConfigReader(const ConfigReader& configreader);

		/**
		* @brief Default destructor. Closes the key/value file if it is open.
		*/
		~ConfigReader() throw(){};


		/**
		* @brief Returns the filename that the ConfigReader object was constructed with.
		* @return std::string The absolute filename of the key/value file.
		*/
		std::string getFileName();


		/**
		* @brief Template function to retrieve a value of class T for a certain key
		* @param key std::string representing a KEY in the key/value file
		* @param t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		*/
		template <class T> void getValue(const std::string& key, T& t) const {
			IOUtils::getValueForKey<T>(properties, key, t);
		}

		/** 
		* Read a line of a .ini config file
		* @param fin   [in] The file input stream to read from.
		* @param lineNb   [in] The line number of the file (used for user-friendly error messages).
		* @param lineType   [out] Returns what sort of line it is (section, key/value, comment, eof).
		* @param str1   [out] Returns a string: section name for section, key for key/value. 
		* @param str2   [out] Returns a string: value for key/value, unset otherwise. 
		* @return true if succesfully read, false if end-of-line or problems encountered
		* 
		*/
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

	private:
		void parseFile();

		std::map<string, string> properties; //Save key value pairs
		std::string filename; //Absolute filename of the key/value file

}; //end class definition ConfigReader

#endif
