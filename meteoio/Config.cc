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
#include <algorithm>

using namespace std;

namespace mio {

const std::string Config::defaultSection = "GENERAL";

//Constructors
Config::Config() : properties(), imported(), sourcename()
{
	//nothing is even put in the property map, the user will have to fill it by himself
}

Config::Config(const std::string& i_filename) : properties(), imported(), sourcename(i_filename)
{
	addFile(i_filename);
}

ConfigProxy Config::get(const std::string& key, const IOUtils::ThrowOptions& opt) const
{
	return ConfigProxy(*this, key, Config::defaultSection, opt);
}

ConfigProxy Config::get(const std::string& key, const std::string& section, const IOUtils::ThrowOptions& opt) const
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

void Config::addKey(const std::string& key, const std::string& section, const std::string& value)
{
	properties[ IOUtils::strToUpper(section) + "::" + IOUtils::strToUpper(key) ] = value;
}

void Config::deleteKey(const std::string& key, const std::string& section)
{
	properties.erase( IOUtils::strToUpper(section) + "::" + IOUtils::strToUpper(key) );
}

bool Config::keyExists(const std::string& key, const std::string& section) const
{
	const string full_key = IOUtils::strToUpper(section) + "::" + IOUtils::strToUpper(key);
	const map<string,string>::const_iterator it = properties.find(full_key);
	return (it!=properties.end());
}

const std::string Config::toString() const {
	std::stringstream os;
	os << "<Config>\n";
	os << "Source: " << sourcename << "\n";
	for (map<string,string>::const_iterator it = properties.begin(); it != properties.end(); ++it){
		os << (*it).first << " -> " << (*it).second << "\n";
	}
	os << "</Config>\n";
	return os.str();
}

std::iostream& operator<<(std::iostream& os, const Config& cfg) {
	const size_t s_source = cfg.sourcename.size();
	os.write(reinterpret_cast<const char*>(&s_source), sizeof(size_t));
	os.write(reinterpret_cast<const char*>(&cfg.sourcename[0]), s_source*sizeof(cfg.sourcename[0]));

	const size_t s_map = cfg.properties.size();
	os.write(reinterpret_cast<const char*>(&s_map), sizeof(size_t));
	for (map<string,string>::const_iterator it = cfg.properties.begin(); it != cfg.properties.end(); ++it){
		const string& key = (*it).first;
		const size_t s_key = key.size();
		os.write(reinterpret_cast<const char*>(&s_key), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&key[0]), s_key*sizeof(key[0]));

		const string& value = (*it).second;
		const size_t s_value = value.size();
		os.write(reinterpret_cast<const char*>(&s_value), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&value[0]), s_value*sizeof(value[0]));
	}

	return os;
}

std::iostream& operator>>(std::iostream& is, Config& cfg) {
	size_t s_source;
	is.read(reinterpret_cast<char*>(&s_source), sizeof(size_t));
	cfg.sourcename.resize(s_source);
	is.read(reinterpret_cast<char*>(&cfg.sourcename[0]), s_source*sizeof(cfg.sourcename[0]));

	cfg.properties.clear();
	size_t s_map;
	is.read(reinterpret_cast<char*>(&s_map), sizeof(size_t));
	for(size_t ii=0; ii<s_map; ii++) {
		size_t s_key, s_value;
		is.read(reinterpret_cast<char*>(&s_key), sizeof(size_t));
		string key;
		key.resize(s_key);
		is.read(reinterpret_cast<char*>(&key[0]), s_key*sizeof(key[0]));

		is.read(reinterpret_cast<char*>(&s_value), sizeof(size_t));
		string value;
		value.resize(s_value);
		is.read(reinterpret_cast<char*>(&value[0]), s_value*sizeof(value[0]));

		cfg.properties[key] = value;
	}
	return is;
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

	std::string section=defaultSection;
	const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
	unsigned int linenr = 0;
	std::vector<std::string> import_after; //files to import after the current one
	bool accept_import_before = true;
	imported.push_back(filename);

	try {
		do {
			std::string line;
			getline(fin, line, eoln); //read complete line
			parseLine(linenr++, import_after, accept_import_before, line, section);
		} while(!fin.eof());
		fin.close();
	} catch(const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}

	std::reverse(import_after.begin(), import_after.end());
	while(!import_after.empty()) {
		const string file_name = import_after.back();
		addFile(file_name);
		import_after.pop_back();
	}
}

void Config::parseLine(const unsigned int& linenr, std::vector<std::string> &import_after, bool &accept_import_before, std::string &line, std::string &section)
{
	//First thing cut away any possible comments (may start with "#" or ";")
	IOUtils::stripComments(line);
	IOUtils::trim(line);    //delete leading and trailing whitespace characters
	if (line.empty()) //ignore empty lines
		return;

	//if this is a section header, read it
	if(line[0] == '[') {
		const size_t endpos = line.find_last_of(']');
		if ((endpos == string::npos) || (endpos < 2) || (endpos != (line.length()-1))) {
			stringstream tmp;
			tmp << linenr;
			throw IOException("Section header corrupt in line " + tmp.str(), AT);
		} else {
			section = line.substr(1,endpos-1);
			IOUtils::toUpper(section);
			return;
		}
	}

	//this can only be a key value pair...
	string key, value;
	if(IOUtils::readKeyValuePair(line, "=", key, value, true)) {
		if(key=="IMPORT_BEFORE") {
			const std::string file_and_path = clean_import_path(value);
			if(!accept_import_before)
				throw IOException("Error in \""+sourcename+"\": IMPORT_BEFORE key MUST occur before any other key!", AT);
			if(std::find(imported.begin(), imported.end(), file_and_path)!=imported.end())
				throw IOException("Can not import \"" + value + "\" again: it has already been imported!", AT);
			parseFile(file_and_path);
			return;
		}
		if(key=="IMPORT_AFTER") {
			const std::string file_and_path = clean_import_path(value);
			if(std::find(imported.begin(), imported.end(), file_and_path)!=imported.end())
				throw IOException("Can not import \"" + value + "\" again: it has already been imported!", AT);
			import_after.push_back(file_and_path);
			return;
		}

		properties[section+"::"+key] = value; //save the key/value pair
		accept_import_before = false; //this is not an import, so no further import_before allowed
	} else {
		stringstream tmp;
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

	if (section.empty()) //enforce the default section if user tries to give empty section string
		section = defaultSection;

	const string tmp_keystart = IOUtils::strToUpper(section) + "::" + IOUtils::strToUpper(keystart);
	//Loop through keys, look for substring match - push it into vecResult
	for (map<string,string>::const_iterator it=properties.begin(); it != properties.end(); ++it){
		const string tmp = (it->first).substr(0, tmp_keystart.length());
		const int matchcount = tmp_keystart.compare(tmp);

		if (matchcount == 0){ //perfect match
			const string tmp2 = (it->first).substr(section.length() + 2);
			vecResult.push_back(tmp2);
		}
	}

	return vecResult.size();
}

std::string Config::extract_section(std::string key) const
{
	const string::size_type pos = key.find("::");

	if (pos != string::npos){
		const string sectionname = key.substr(0, pos);
		key.erase(key.begin(), key.begin() + pos + 2); //delete section name
		return sectionname;
	}
	return defaultSection;
}

std::string Config::clean_import_path(const std::string& in_path) const
{
#ifdef WIN32
	const bool is_absolute = (in_path.size()>=2 && in_path[1]==':');
#else
	const bool is_absolute = (in_path.size()>=1 && in_path[0]=='/');
#endif
	//if this is a relative path, prefix the import path with the current path
	const std::string prefix = (!is_absolute)? IOUtils::getPath(sourcename, true)+"/" : "";
	const std::string path = IOUtils::getPath(prefix+in_path, true);  //clean & resolve path
	const std::string filename = IOUtils::getFilename(in_path);

	return path + "/" + filename;
}

void Config::write(const std::string& filename) const
{

	std::ofstream fout;
	fout.open(filename.c_str());
	if (fout.fail())
		throw FileAccessException(filename, AT);

	try {
		string current_section;
		unsigned int sectioncount = 0;
		for (map<string,string>::const_iterator it=properties.begin(); it != properties.end(); ++it){
			const string key = it->first;
			const string section = extract_section(key);

			if (current_section != section){
				current_section = section;
				if (sectioncount != 0)
					fout << endl;
				sectioncount++;
				fout << "[" << section << "]" << endl;
			}

			fout << key << " = " << it->second << endl;
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
