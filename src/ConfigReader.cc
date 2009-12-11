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
#include "ConfigReader.h"

using namespace std;

const unsigned int ConfigReader::nothrow = 666;

ConfigReader::ConfigReader(const std::string& filename_in/*, const std::string& delimiter*/)
{
	//Check whether file exists and is accessible -> throw FileNotFound or FileNotAccessible Exception
	//Read Key Value Pairs and put them into map, depending on the section they are in
	//If there is no section definition, the default section "General" is assumed
	filename = filename_in;
	parseFile();
} 

void ConfigReader::parseFile()
{
	std::ifstream fin; //Input file streams
	unsigned int linenr = 0;
	std::string line="", section="GENERAL";
	
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

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
	
	try {
		do {
			getline(fin, line, eoln); //read complete line
			parseLine(linenr++, line, section);
		} while(!fin.eof());
		fin.close();
	} catch(std::exception& e){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
}

void ConfigReader::parseLine(const unsigned int& linenr, std::string& line, std::string& section)
{
	stringstream tmp;

	IOUtils::trim(line); //delete leading and trailing whitespace characters
	if (line.length() == 0) //ignore empty lines
		return;
	
	if ((line[0]=='#') || (line[0]==';')) //check whether line is a comment
		return;

	if (line[0] == '['){
		size_t endpos = line.find_last_of(']');
		if ((endpos == string::npos) || (endpos < 2) || (endpos != (line.length()-1))){
			tmp << linenr;
			throw IOException("Section header corrupt in line " + tmp.str(), AT);
		} else {
			section = line.substr(1,endpos-1);
			IOUtils::toUpper(section);
			return;
		}
	}

	//At this point line can only be a key value pair
	if (!IOUtils::readKeyValuePair(line, "=", properties, section+"::")){
		tmp << linenr;
		throw InvalidFormatException("Error reading key value pair in " + filename + " line:" + tmp.str(), AT);    
	}
}

//Return key/value filename
std::string ConfigReader::getFileName()
{
	return filename;
}

unsigned int ConfigReader::findKeys(vector<string>& vecResult, 
							const std::string keystart, 
							const std::string section)
{
	vecResult.clear();
	string _section = section;
	if (_section.length() == 0) //enforce the default section if user tries to give empty section string
		_section = "GENERAL";

	IOUtils::toUpper(_section);
	string _keystart = _section + "::" + keystart;

	//Loop through keys, look for substring match - push it into vecResult
	map<string,string>::iterator it;
	for (it=properties.begin(); it != properties.end(); it++){
		string tmp = (*it).first;
		tmp = tmp.substr(0,_keystart.length());

		int matchcount = _keystart.compare(tmp);

		if (matchcount == 0){ //perfect match
			string tmp2 = it->first;
			tmp2 = tmp2.substr(_section.length() + 2);
			vecResult.push_back(tmp2);
		}
	}

	return vecResult.size();
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
