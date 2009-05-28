#include "ConfigReader.h"

using namespace std;

//DECONSTRUCTOR
ConfigReader::~ConfigReader() throw()
{
	//properties.clear();
	//file.~string();

	if (fin.is_open()) {//close fin if open
		fin.close();
	}
}

//CONSTRUCTORS
ConfigReader::ConfigReader(const ConfigReader& configreader) : properties(configreader.properties), file(configreader.file){}

ConfigReader::ConfigReader(const std::string& filename/*, const std::string& delimiter*/)
{
	//Check whether file exists and is accessible -> throw FileNotFound or FileNotAccessible Exception
	//Read first line -> if not [Parameters] -> throw WrongFormatException
	//Read Key Value Pairs and put them into map 
	file = filename;
	if (!IOUtils::validFileName(filename)) {
		THROW InvalidFileNameException(filename,AT);
	}

	int linenr = 0;

	//Check whether file exists
	if (!IOUtils::fileExists(filename)) {
		THROW FileNotFoundException(filename, AT);
	}

	//Open file
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		THROW FileAccessException(filename, AT);
	}

	int lineType = ConfigReader::CfgLineComment /*dummy non-problematic initial value*/;
	string str1, str2;
	stringstream tmpStringStream;

	//Read header, should be [Parameters]
	ConfigReader::readConfigLine(fin, ++linenr, lineType, str1, str2);
	if (lineType != ConfigReader::CfgLineSection || str1 != "Parameters") {
		THROW InvalidFormatException("Missing header [Parameters] in file " + filename, AT);    
	}
	//Go through file, save key value pairs
	while (ConfigReader::readConfigLine(fin, ++linenr, lineType, str1, str2)) {
		if (lineType == ConfigReader::CfgLineKeyValue) {
			properties[str1] = str2;
		} else if (lineType != ConfigReader::CfgLineComment) {
			tmpStringStream << "Expected key-value pair at line "<<linenr<<".";
			THROW InvalidFormatException(tmpStringStream.str(), AT);
		}
	}

	fin.close();
} // end ConfigReader::ConfigReader

  
bool ConfigReader::readConfigLine(std::istream& fin, int lineNb, int& lineType, string& str1, string& str2)
{
	string line="";
	unsigned int pos = 0;
	stringstream tmpStringStream;
	bool isSuccess = false;

	str1 = "";
	str2 = "";
	if (! std::getline(fin, line)) {
		lineType = CfgLineEOF;
		//cout << "[D] readConfigLine lineNb="<<lineNb<<" lineType="<<lineType<<" str1='"<<str1<<"' str2='"<<str2<<"'."<<endl;
		return isSuccess;
	}
	IOUtils::trim(line);

	if (line[0] == ';' || line[0] == '#') {
		// handle comment
		lineType = CfgLineComment;
		str1 = line;
		isSuccess = true;

	} else if (line[0] == '[') {
		// handle sections
		lineType = CfgLineSection;
		pos = line.find_first_of(']');
		if (pos != line.npos) {
			str1 = line.substr(1, pos-1);
			line = line.substr(pos+1);
			IOUtils::trim(line);
			if (line.length() > 0 && line[0] != ';' && line[0] != '#') {
				tmpStringStream << "Invalid non-comment data at line "<<lineNb<<" after section ["<<str1<<"].";
				THROW InvalidFormatException(tmpStringStream.str(), AT);
			}
			isSuccess = true;
		} else {
			tmpStringStream << "Config section name at line "<<lineNb<<" not properly ended with ']'.";
			THROW InvalidFormatException(tmpStringStream.str(), AT);
		}

	} else if (line.length() > 0) {
		pos = line.find_first_of('=');
		if (pos != line.npos) {
			// handle key/value pair
			lineType = CfgLineKeyValue;
			str1 = line.substr(0, pos);
			str2 = line.substr(pos + 1);
			IOUtils::trim(str1);
			IOUtils::trim(str2);
			if (str1 == "") {
				tmpStringStream << "Unallowed empty key at line "<<lineNb<<" of config file.";
				THROW InvalidFormatException(tmpStringStream.str(), AT);
			}
			if (str2 == "") {
				tmpStringStream << "Unallowed empty value at line "<<lineNb<<" of config file.";
				THROW InvalidFormatException(tmpStringStream.str(), AT);
			}
			isSuccess = true;

		} else {
			// handle unknown
			lineType = CfgLineUnknown;
			str1 = line;
			tmpStringStream << "Unexpected content at line "<<lineNb<<" of config file.";
			THROW InvalidFormatException(tmpStringStream.str(), AT);
		}

	} else {
		// handle empty line
		lineType = CfgLineComment;
		isSuccess = true;

	}
	//cout << "[D] readConfigLine lineNb="<<lineNb<<" lineType="<<lineType<<" str1='"<<str1<<"' str2='"<<str2<<"'."<<endl;
	return isSuccess;
}

//Return key/value filename
std::string ConfigReader::getFileName()
{
	return file;
}
