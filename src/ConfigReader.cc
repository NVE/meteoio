#include "ConfigReader.h"

using namespace std;

//CONSTRUCTORS
#ifdef _POPC_
ConfigReader::ConfigReader(const ConfigReader& configreader) : paroc_base(), properties(configreader.properties), filename(configreader.filename){}
#else
ConfigReader::ConfigReader(const ConfigReader& configreader) : properties(configreader.properties), filename(configreader.filename){}
#endif

ConfigReader::ConfigReader(const std::string& filename_in/*, const std::string& delimiter*/)
{
	//Check whether file exists and is accessible -> throw FileNotFound or FileNotAccessible Exception
	//Read first line -> if not [Parameters] -> throw WrongFormatException
	//Read Key Value Pairs and put them into map 
	filename = filename_in;
	parseFile();

} // end ConfigReader::ConfigReader

void ConfigReader::parseFile()
{
	std::ifstream fin; //Input file streams
	int linenr = 0;

	if (!IOUtils::validFileName(filename)) {
		throw InvalidFileNameException(filename,AT);
	}

	//Check whether file exists
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}

	//Open file
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}

	try {
		int lineType = ConfigReader::CfgLineComment /*dummy non-problematic initial value*/;
		string str1, str2;
		stringstream tmpStringStream;
	
		//Read header, should be [Parameters]
		ConfigReader::readConfigLine(fin, ++linenr, lineType, str1, str2);
		if (lineType != ConfigReader::CfgLineSection || str1 != "Parameters") {
			throw InvalidFormatException("Missing header [Parameters] in file " + filename, AT);    
		}
		//Go through file, save key value pairs
		while (ConfigReader::readConfigLine(fin, ++linenr, lineType, str1, str2)) {
			if (lineType == ConfigReader::CfgLineKeyValue) {
				properties[str1] = str2;
			} else if (lineType != ConfigReader::CfgLineComment) {
				tmpStringStream << "Expected key-value pair at line "<<linenr<<".";
				throw InvalidFormatException(tmpStringStream.str(), AT);
			}
		}

		//done reading, so closing the file
		fin.close();
	} catch (...) {
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
}

bool ConfigReader::readConfigLine(std::istream& fin, int lineNb, int& lineType, string& str1, string& str2)
{
	string line="";
	std::string::size_type pos = 0;
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
				throw InvalidFormatException(tmpStringStream.str(), AT);
			}
			isSuccess = true;
		} else {
			tmpStringStream << "Config section name at line "<<lineNb<<" not properly ended with ']'.";
			throw InvalidFormatException(tmpStringStream.str(), AT);
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
				throw InvalidFormatException(tmpStringStream.str(), AT);
			}
			if (str2 == "") {
				tmpStringStream << "Unallowed empty value at line "<<lineNb<<" of config file.";
				throw InvalidFormatException(tmpStringStream.str(), AT);
			}
			isSuccess = true;

		} else {
			// handle unknown
			lineType = CfgLineUnknown;
			str1 = line;
			tmpStringStream << "Unexpected content at line "<<lineNb<<" of config file.";
			throw InvalidFormatException(tmpStringStream.str(), AT);
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
	return filename;
}

#ifdef _POPC_
#include "marshal_meteoio.h"
void ConfigReader::Serialize(POPBuffer &buf, bool pack)
{
	if (pack)
	{
		buf.Pack(&filename,1);
		marshal_map_str_str(buf, properties, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&filename,1);
		marshal_map_str_str(buf, properties, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif
