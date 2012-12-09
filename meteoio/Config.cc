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
#include <meteoio/Config.h>

using namespace std;

namespace mio {

const std::string Config::defaultSection = "GENERAL";

//Constructors
Config::Config() : properties(), sourcename()
{
	//nothing is even put in the property map, the user will have to fill it by himself
}

Config::Config(const std::string& i_filename) : properties(), sourcename(i_filename)
{
	addFile(i_filename);
}

ConfigProxy Config::get(const std::string& key, const Options& opt) const
{
	return ConfigProxy(*this, key, Config::defaultSection, opt);
}

ConfigProxy Config::get(const std::string& key, const std::string& section, const Options& opt) const
{
	return ConfigProxy(*this, key, section, opt);
}

//Populating the property map
void Config::addFile(const std::string& i_filename)
{
	sourcename = i_filename;
	parseFile(i_filename);
}

void Config::addCmdLine(const std::string& cmd_line)
{
	sourcename = std::string("Command line");
	parseCmdLine(cmd_line);
}

void Config::addKey(const std::string& key, const std::string& value)
{
	const std::string section=defaultSection;
	addKey(key, section, value);
}

void Config::addKey(std::string key, std::string section, const std::string& value)
{
	IOUtils::toUpper(key);
	IOUtils::toUpper(section);
	properties[section + "::" + key] = value;
}

void Config::deleteKey(std::string key, std::string section)
{
	IOUtils::toUpper(key);
	IOUtils::toUpper(section);
	properties.erase(section + "::" + key);
}

std::ostream& operator<<(std::ostream &os, const Config& cfg)
{
	os << "<Config>\n";
	map<string,string>::const_iterator it;
	for (it=cfg.properties.begin(); it != cfg.properties.end(); it++){
		os << (*it).first << " -> " << (*it).second << "\n";
	}
	os << "</Config>\n";
	return os;
}

//Parsing
void Config::parseCmdLine(const std::string& cmd_line)
{
	(void)cmd_line;
	throw IOException("Nothing implemented here", AT);
}

void Config::parseFile(const std::string& filename)
{
	std::ifstream fin; //Input file streams

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

	std::string line="", section=defaultSection;
	const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
	unsigned int linenr = 0;

	try {
		do {
			getline(fin, line, eoln); //read complete line
			parseLine(linenr++, line, section);
		} while(!fin.eof());
		fin.close();
	} catch(const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
}

void Config::parseLine(const unsigned int& linenr, std::string& line, std::string& section)
{
	//First thing cut away any possible comments (may start with "#" or ";")
	IOUtils::stripComments(line);

	IOUtils::trim(line);    //delete leading and trailing whitespace characters
	if (line.length() == 0) //ignore empty lines
		return;

	stringstream tmp;       //stringstream to convert the unsigned int linenr into a string
	if (line[0] == '['){
		const size_t endpos = line.find_last_of(']');
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
	if (!IOUtils::readKeyValuePair(line, "=", properties, section+"::", true)){
		tmp << linenr;
		throw InvalidFormatException("Error reading key value pair in \"" + sourcename + "\" at line " + tmp.str(), AT);
	}
}

//Return key/value filename
std::string Config::getSourceName() const
{
	return sourcename;
}

size_t Config::findKeys(std::vector<std::string>& vecResult, const std::string& keystart,
                        std::string section) const
{
	vecResult.clear();

	if (section.length() == 0) //enforce the default section if user tries to give empty section string
		section = defaultSection;

	const string tmp_keystart = IOUtils::strToUpper(section) + "::" + IOUtils::strToUpper(keystart);
	//Loop through keys, look for substring match - push it into vecResult
	for (map<string,string>::const_iterator it=properties.begin(); it != properties.end(); it++){
		const string tmp = (it->first).substr(0, tmp_keystart.length());
		const int matchcount = tmp_keystart.compare(tmp);

		if (matchcount == 0){ //perfect match
			const string tmp2 = (it->first).substr(section.length() + 2);
			vecResult.push_back(tmp2);
		}
	}

	return vecResult.size();
}

std::string Config::extract_section(std::string& key)
{
	const string::size_type pos = key.find("::");

	if (pos != string::npos){
		const string sectionname = key.substr(0, pos);
		key.erase(key.begin(), key.begin() + pos + 2); //delete section name
		return sectionname;
	}
	return defaultSection;
}

void Config::write(const std::string& filename)
{

	std::ofstream fout;
	fout.open(filename.c_str());
	if (fout.fail())
		throw FileAccessException(filename, AT);

	try {
		string current_section;
		unsigned int sectioncount = 0;
		for (map<string,string>::const_iterator it=properties.begin(); it != properties.end(); it++){
			string tmp = it->first;
			const string section = extract_section(tmp);

			if (current_section != section){
				current_section = section;
				if (sectioncount != 0)
					fout << endl;
				sectioncount++;
				fout << "[" << section << "]" << endl;
			}

			fout << tmp << " = " << it->second << endl;
		}
	} catch(...) {
		if (fout.is_open()) //close fout if open
			fout.close();

		throw;
	}

	if (fout.is_open()) //close fout if open
		fout.close();
}

} //end namespace

#ifdef _POPC_
#include "marshal_meteoio.h"
using namespace mio; //HACK for POPC
void Config::Serialize(POPBuffer &buf, bool pack)
{
	if (pack)
	{
		buf.Pack(&sourcename,1);
		marshal_map_str_str(buf, properties, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		buf.UnPack(&sourcename,1);
		marshal_map_str_str(buf, properties, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

//} //namespace //HACK for POPC
