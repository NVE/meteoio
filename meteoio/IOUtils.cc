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
#include <cmath>
#include <cstring>
#include <algorithm>

#include <meteoio/IOUtils.h>
#include <meteoio/meteolaws/Meteoconst.h> //for math constants
#include <meteoio/Config.h>    // to avoid forward declaration hell
#include <meteoio/MeteoData.h> // to avoid forward declaration hell

#ifdef WIN32
	#include <windows.h>
	#include "Shlwapi.h"
	//removing two macros defined in windows.h
	#undef max
	#undef min
#else
	#include <dirent.h>
	#include <sys/stat.h>
#endif

namespace mio {

#ifdef _MSC_VER
//This is C99, Microsoft should move on and suppport it, it is almost 15 years old!!
double round(const double& x) {
	//middle value point test
	if (ceil(x+0.5) == floor(x+0.5)) {
		const int a = (int)ceil(x);
		if (a%2 == 0) {
			return ceil(x);
		} else {
			return floor(x);
		}
	} else {
		return floor(x+0.5);
	}
}
#endif

std::string getLibVersion() {
	std::stringstream ss;
	ss << _VERSION << " compiled on " << __DATE__ << " " << __TIME__;
	return ss.str();
}

namespace IOUtils {

double bearing_to_angle(const double& bearing) {
	return (fmod(360.-bearing+90., 360.)*Cst::to_rad);
}

double angle_to_bearing(const double& angle) {
	return (fmod( 90.-angle*Cst::to_deg+360. , 360. ));
}

double bearing(std::string bearing_str)
{
	trim(bearing_str);
	toUpper(bearing_str);

	if(bearing_str=="N") return 0.;
	if(bearing_str=="NNE") return 22.5;
	if(bearing_str=="NE") return 45.;
	if(bearing_str=="ENE") return 67.5;
	if(bearing_str=="E") return 90.;
	if(bearing_str=="ESE") return 112.5;
	if(bearing_str=="SE") return 135.;
	if(bearing_str=="SSE") return 157.5;
	if(bearing_str=="S") return 180.;
	if(bearing_str=="SSW") return 202.5;
	if(bearing_str=="SW") return 225.;
	if(bearing_str=="WSW") return 247.5;
	if(bearing_str=="W") return 270.;
	if(bearing_str=="WNW") return 292.5;
	if(bearing_str=="NW") return 315.;
	if(bearing_str=="NNW") return 337.5;

	//no match
	return nodata;
}

void stripComments(std::string& str)
{
	const size_t found = str.find_first_of("#;");
	if (found != std::string::npos){
		str.erase(found); //rest of line disregarded
	}
}

void copy_file(const std::string& src, const std::string& dest)
{
	if (src == dest) return; //copying to the same file doesn't make sense, but is no crime either

	std::ifstream fin(src.c_str(), std::ios::binary);
	if (fin.fail()) throw FileAccessException(src, AT);

	std::ofstream fout(dest.c_str(), std::ios::binary);
	if (fout.fail()) {
		fin.close();
		throw FileAccessException(dest, AT);
	}

	fout << fin.rdbuf();

	fin.close();
	fout.close();
}

std::string cleanPath(const std::string& in_path, const bool& resolve)
{
	if(!resolve) { //do not resolve links, relative paths, etc
		std::string out_path(in_path);
		std::replace(out_path.begin(), out_path.end(), '\\', '/');
		return out_path;
	} else {
	#ifdef WIN32
		//if this would not suffice, see http://pdh11.blogspot.ch/2009/05/pathcanonicalize-versus-what-it-says-on.html
		char **ptr = NULL;
		char *out_buff = (char*)calloc(MAX_PATH, sizeof(char));
		const DWORD status = GetFullPathName(in_path.c_str(), MAX_PATH, out_buff, ptr);
		std::string out_path = (status!=0 && status<=MAX_PATH)? out_buff : in_path;
		free(out_buff);

		std::replace(out_path.begin(), out_path.end(), '\\', '/');
		return out_path;
	#else //POSIX
		std::string out_path(in_path);
		std::replace(out_path.begin(), out_path.end(), '\\', '/');

		char *real_path = realpath(out_path.c_str(), NULL); //POSIX
		if(real_path!=NULL) {
			const std::string tmp(real_path);
			free(real_path);
			return tmp;
		} else
			return out_path; //something failed in realpath, keep it as it is
	#endif
	}
}

std::string getExtension(const std::string& filename)
{
	const size_t start_basename = filename.find_last_of("/\\"); //we will skip the path
	const size_t startpos = filename.find_last_of('.');
	if( startpos==std::string::npos ) return std::string();
	if( start_basename!=std::string::npos && startpos<start_basename ) return std::string();

	const std::string whitespaces(" \t\f\v\n\r");
	const size_t endpos = filename.find_last_not_of(whitespaces); // Find the first character position from reverse af

	return filename.substr(startpos+1, endpos-startpos);
}

std::string removeExtension(const std::string& filename)
{
	const size_t start_basename = filename.find_last_of("/\\"); //we will skip the path
	const size_t startpos = filename.find_last_of('.');
	if( startpos==std::string::npos ) return filename;
	if( start_basename!=std::string::npos && startpos<start_basename ) return filename;

	return filename.substr(0, startpos);
}

std::string getPath(const std::string& filename, const bool& resolve)
{
	const std::string clean_filename = cleanPath(filename, resolve);
	const size_t end_path = clean_filename.find_last_of("/");
	if(end_path!=std::string::npos) {
		return clean_filename.substr(0, end_path);
	} else {
		return cleanPath("./", resolve);
	}
}

std::string getFilename(const std::string& path)
{
	const size_t start_basename = path.find_last_of("/\\");
	if(start_basename!=std::string::npos)
		return path.substr(start_basename+1, std::string::npos);
	else
		return path;
}

void trim(std::string& str)
{
	const std::string whitespaces(" \t\f\v\n\r");
	const size_t startpos = str.find_first_not_of(whitespaces); // Find the first character position after excluding leading blank spaces
	const size_t endpos = str.find_last_not_of(whitespaces); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(startpos!=std::string::npos && endpos!=std::string::npos) {
		str.erase(endpos+1); //right trim
		str.erase(0, startpos); //left trim
	} else {
		str.clear();
	}
}

void toUpper(std::string& str) {
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

void toLower(std::string& str) {
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

std::string strToUpper(const std::string &str) {
    std::string strcopy(str.size(), 0);
    std::transform(str.begin(), str.end(), strcopy.begin(), ::toupper);
    return strcopy;
}

std::string strToLower(const std::string &str) {
    std::string strcopy(str.size(), 0);
    std::transform(str.begin(), str.end(), strcopy.begin(), ::tolower);
    return strcopy;
}

bool isNumeric(std::string str, const unsigned int& nBase)
{
	trim(str); //delete trailing and leading whitespaces and tabs
	std::istringstream iss(str);

	if( nBase == 10 ) {
		double tmp;
		iss >> tmp;
	} else if( nBase == 8 || nBase == 16 ) {
		int tmp;
		iss >> ( ( nBase == 8 ) ? std::oct : std::hex ) >> tmp;
	} else
		return false;

	if( !iss ) //the conversion failed
		return false;

	return ( iss.rdbuf()->in_avail() == 0 ); //true if nothing was left after conversion
}

bool readKeyValuePair(const std::string& in_line, const std::string& delimiter, std::string &key, std::string &value, const bool& setToUpperCase)
{
	const size_t pos = ((delimiter==" ") || (delimiter=="\t"))? in_line.find_first_of(" \t", 0) : in_line.find(delimiter); //first occurence of delimiter

	if(pos != std::string::npos) { //ignore in_lines that are empty or without '='
		key = in_line.substr(0, pos);
		value = in_line.substr(pos + 1);

		trim(key);
		trim(value);

		if (key.empty() || value.empty()) {
			return false;
		}

		if (setToUpperCase)
			toUpper(key);
	} else {
		key="";
		value="";
		return false;
	}

	return true;
}

bool validFileName(const std::string& filename)
{
	const size_t startpos = filename.find_first_not_of(" \t\n"); // Find the first character position after excluding leading blank spaces
	if((startpos!=0) || (filename==".") || (filename=="..")) {
		return false;
	}

	return true;
}

#ifdef WIN32
bool fileExists(const std::string& filename)
{
	return ( GetFileAttributes( filename.c_str() ) != INVALID_FILE_ATTRIBUTES );
}

void readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	const size_t path_length = path.length();
	if (path_length > (MAX_PATH - 3)) {
		std::cerr << "Path " << path << "is too long (" << path_length << " characters)" << std::endl;
		throw FileAccessException("Error opening directory " + path, AT);
	}

	const std::string filepath = path+"\\"+pattern;
	WIN32_FIND_DATA ffd;
	const HANDLE hFind = FindFirstFile(filepath.c_str(), &ffd);
	if (INVALID_HANDLE_VALUE == hFind) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	do {
		if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			//this is a directory -> do nothing
		} else {
			const std::string filename(ffd.cFileName);
			dirlist.push_back(filename);
		}
	}
	while (FindNextFile(hFind, &ffd) != 0);

	const DWORD dwError = GetLastError();
	if (dwError != ERROR_NO_MORE_FILES) {
		throw FileAccessException("Error listing files in directory " + path, AT);
	}

	FindClose(hFind);
}
#else
bool fileExists(const std::string& filename)
{
	struct stat buffer ;

	if ((stat( filename.c_str(), &buffer))==0) {//File exists if stat returns 0
		return true ;
	}

	return false;
}

void readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	DIR *dp = opendir(path.c_str());
	if(dp == NULL) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	struct dirent *dirp;
	while ((dirp = readdir(dp)) != NULL) {
		const std::string tmp(dirp->d_name);
		if( tmp.compare(".")!=0 && tmp.compare("..")!=0 ) { //skip "." and ".."
			if (pattern.empty()) {
				dirlist.push_back(tmp);
			} else {
				const size_t pos = tmp.find(pattern);
				if (pos!=std::string::npos) {
					dirlist.push_back(tmp);
				}
			}
		}
	}
	closedir(dp);
}
#endif

std::string getLogName() {
	char *tmp;

	if((tmp=getenv("USERNAME"))==NULL) { //Windows & Unix
		if((tmp=getenv("LOGNAME"))==NULL) { //Unix
			tmp=getenv("USER"); //Windows & Unix
		}
	}

	if(tmp==NULL) return std::string("N/A");
	return std::string(tmp);
}

void readKeyValueHeader(std::map<std::string,std::string>& headermap,
                                 std::istream& fin,
                                 const size_t& linecount,
                                 const std::string& delimiter)
{
	size_t linenr = 0;
	std::string line;

	//make a test for end of line encoding:
	const char eol = getEoln(fin);

	for (size_t ii=0; ii< linecount; ii++){
		if (std::getline(fin, line, eol)) {
			std::string key, value;
			linenr++;
			const bool result = readKeyValuePair(line, delimiter, key, value);
			if(result) {
				headermap[key] = value;
			} else { //  means if ((key == "") || (value==""))
				std::stringstream out;
				out << "Invalid key value pair in line: " << linenr << " of header";
				throw IOException(out.str(), AT);
			}
		} else {
			throw InvalidFormatException("Premature EOF while reading Header", AT);
		}
	}
}

char getEoln(std::istream& fin)
{
	std::streambuf* pbuf;
	char tmp = '0';
	int chars = 0;

	const std::streampos file_start = fin.tellg();

	do {
		fin.get(tmp);
		chars++;

		if ((tmp == '\r') || (tmp == '\n')) {
			char peekc = tmp;
			while ((!fin.eof() && ((peekc=='\r') || (peekc=='\n')))) {
				tmp = peekc;
				fin.get(peekc);
				chars++;
			}
			pbuf = fin.rdbuf();
			pbuf->pubseekpos(file_start); //rewind
			fin.clear(); //reset eof flag, etc
			return tmp;
		}
	} while ((chars < 3000) && (!fin.eof()));

	pbuf = fin.rdbuf();
	pbuf->pubseekpos(file_start); //rewind
	fin.clear(); //reset eof flag, etc

	return '\n';
}

void skipLines(std::istream& fin, const size_t& nbLines, const char& eoln)
{
	std::string dummy;
	for (size_t ii=0; ii<nbLines; ii++) {
		if(!getline(fin, dummy, eoln)) {
			throw InvalidFormatException("Premature EOF while skipping lines", AT);
		}
	}
}

size_t readLineToVec(const std::string& line_in, std::vector<double>& vec_data)
{
	vec_data.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input
	iss.setf(std::ios::fixed);
	iss.precision(std::numeric_limits<double>::digits10);

	double tmp;
	while (!iss.eof()) {
		iss >> std::skipws >> tmp;

		if (iss.fail()) {
			std::stringstream ss;
			ss << "Can not read column " << vec_data.size()+1 << " in data line \"" << line_in << "\"";
			throw InvalidFormatException(ss.str(), AT);
		}
		vec_data.push_back(tmp);
	}

	return vec_data.size();
}

size_t readLineToVec(const std::string& line_in, std::vector<std::string>& vecString)
{
	vecString.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input

	std::string tmp_string;
	while (!iss.eof()) {
		iss >> std::skipws >> tmp_string;

		if (!tmp_string.empty()) {
			vecString.push_back(tmp_string);
			tmp_string.clear();
		}
	}

	return vecString.size();
}

size_t readLineToVec(const std::string& line_in, std::vector<std::string>& vecString, const char& delim)
{
	vecString.clear();
	std::string tmp_string;
	std::istringstream iss(line_in);

	while (getline(iss, tmp_string, delim)){
		vecString.push_back(tmp_string);
	}

	return vecString.size();
}

// generic template function convertString must be defined in the header

const char ALPHANUM[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
const char NUM[] = "0123456789";

template<> bool convertString<std::string>(std::string& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	(void)f;
	t = str;
	trim(t); //delete trailing and leading whitespaces and tabs
	return true;
}

template<> bool convertString<bool>(bool& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s(str);
	trim(s); //delete trailing and leading whitespaces and tabs

	if (toupper(s[0])=='T' || toupper(s[0])=='Y') {
		t = true;
	} else if (toupper(s[0])=='F' || toupper(s[0])=='N') {
		t = false;
	} else {
		std::istringstream iss(s);
		int i;
		iss >> f >> i; //Convert first part of stream with the formatter (e.g. std::dec, std::oct)
		if (iss.fail()) {//Conversion failed
			return false;
		}
		t = (i != 0);
	}

	const std::string::size_type pos = s.find_first_not_of(ALPHANUM);
	if (pos != std::string::npos) {
		std::string tmp = s.substr(pos);
		trim(tmp);
		if (!tmp.empty() && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
			return false;
		}
	}

	return true;
}

template<> bool convertString<unsigned int>(unsigned int& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s(str);
	trim(s); //delete trailing and leading whitespaces and tabs
	if (s.empty()) {
		t = unodata;
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
		std::string tmp;
		getline(iss,  tmp); //get rest of line, if any
		trim(tmp);
		if (!tmp.empty() && tmp[0] != '#' && tmp[0] != ';') {
			//if line holds more than one value it's invalid
			return false;
		}
		return true;
	}
}

bool convertString(Date& t, const std::string& str, const double& time_zone, std::ios_base& (*f)(std::ios_base&))
{
	std::string s(str);
	trim(s); //delete trailing and leading whitespaces and tabs

	(void)f;
	int year;
	unsigned int month, day, hour, minute;
	double second;
	char rest[32] = "";

	//HACK: we read the seconds, but we ignore them...
	const char *c_str = s.c_str();
	if (sscanf(c_str, "%d-%u-%u %u:%u:%lg%31s", &year, &month, &day, &hour, &minute, &second, rest) >= 6) {
		std::string timezone_iso(rest);
		stripComments(timezone_iso);
		const double tz = (timezone_iso.empty())? time_zone : Date::parseTimeZone(timezone_iso);
		if(tz==nodata) return false;
		t.setDate(year, month, day, hour, minute, tz);
		return true;

	} else if (sscanf(c_str, "%d-%u-%uT%u:%u:%lg%31s", &year, &month, &day, &hour, &minute, &second, rest) >= 6) { //ISO
		std::string timezone_iso(rest);
		stripComments(timezone_iso);
		const double tz = (timezone_iso.empty())? time_zone : Date::parseTimeZone(timezone_iso);
		if(tz==nodata) return false;
		t.setDate(year, month, day, hour, minute, tz);
		return true;

	} else if (sscanf(c_str, "%d-%u-%u %u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
		std::string timezone_iso(rest);
		stripComments(timezone_iso);
		const double tz = (timezone_iso.empty())? time_zone : Date::parseTimeZone(timezone_iso);
		if(tz==nodata) return false;
		t.setDate(year, month, day, hour, minute, tz);
		return true;

	} else if (sscanf(c_str, "%d-%u-%uT%u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
		std::string timezone_iso(rest);
		stripComments(timezone_iso);
		const double tz = (timezone_iso.empty())? time_zone : Date::parseTimeZone(timezone_iso);
		if(tz==nodata) return false;
		t.setDate(year, month, day, hour, minute, tz);
		return true;

	} else if (sscanf(c_str, "%d-%u-%u%31s", &year, &month, &day, rest) >= 3) {
		std::string timezone_iso(rest);
		stripComments(timezone_iso);
		const double tz = (timezone_iso.empty())? time_zone : Date::parseTimeZone(timezone_iso);
		if(tz==nodata) return false;
		t.setDate(year, month, day, (unsigned)0, (unsigned)0, tz);
		return true;

	} else if (sscanf(c_str, "%u:%u%31s", &hour, &minute, rest) >= 2) {
		std::string timezone_iso(rest);
		stripComments(timezone_iso);
		const double tz = (timezone_iso.empty())? time_zone : Date::parseTimeZone(timezone_iso);
		if(tz==nodata) return false;
		t.setDate( ((double)hour)/24. + ((double)minute)/24./60. , tz);
		return true;

	} else {
		//try to read purely numerical date, potentially surrounded by other chars
		//and potentially containing an ISO time zone string
		const size_t in_len = s.length();

		//extract date/time
		const size_t date_beg = s.find_first_of(NUM);
		if(date_beg==npos || date_beg==in_len) return false;
		size_t date_end = s.find_first_not_of(NUM, date_beg+1);
		if(date_end==npos) date_end = in_len;
		const std::string date = s.substr(date_beg, date_end-date_beg);

		//parse date/time
		const size_t date_len = date.length();
		if(date_len<10 || date_len>14) return false;
		if( convertString(year,date.substr(0,4))==false ) return false;
		if( convertString(month,date.substr(4,2))==false ) return false;
		if( convertString(day,date.substr(6,2))==false ) return false;
		if( convertString(hour,date.substr(8,2))==false ) return false;
		if(date_len==10)
			minute=0;
		else {
			if(date_len>=12) {
				if( convertString(minute,date.substr(10,2))==false ) return false;
			} else
				return false;
			if(date_len==12)
				second=0;
			else {
				if(date_len==14) {
					if( convertString(second,date.substr(12,2))==false ) return false;
				} else
					return false;
			}
		}

		//extract potential ISO time zone string
		double tz = time_zone;
		const size_t tz_beg = s.find_first_of("+-", date_end);
		if(tz_beg!=npos && tz_beg!=in_len) {
			size_t tz_end = s.find_first_not_of("0123456789:", date_end+1);
			if(tz_end==npos) tz_end = in_len;
			const std::string timezone_iso = s.substr(tz_beg, tz_end-tz_beg);
			if(!timezone_iso.empty()) tz = Date::parseTimeZone(timezone_iso);
		}

		t.setDate( year, month, day, hour, minute, tz );
	}

	return true;
}

template<> bool convertString<Coords>(Coords& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s(str);
	trim(s); //delete trailing and leading whitespaces and tabs

	(void)f;
	double lat, lon;
	try {
		Coords::parseLatLon(s, lat, lon);
	} catch(const IOException&) {
		return false;
	}
	t.setLatLon(lat, lon, nodata);

	return true;
}


void getProjectionParameters(const Config& cfg, std::string& coordin, std::string& coordinparam,
                                      std::string& coordout, std::string& coordoutparam) {
	cfg.getValue("COORDSYS", "Input", coordin);
	cfg.getValue("COORDPARAM", "Input", coordinparam, IOUtils::nothrow);
	cfg.getValue("COORDSYS", "Output", coordout, IOUtils::nothrow);
	cfg.getValue("COORDPARAM", "Output", coordoutparam, IOUtils::nothrow);
}

void getTimeZoneParameters(const Config& cfg, double& tz_in, double& tz_out) {
	cfg.getValue("TIME_ZONE", "Input", tz_in, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", tz_out, IOUtils::nothrow);
}

//returns index of element, if element does not exist it returns closest index after soughtdate
//the element needs to be an exact hit or embedded between two measurments
size_t seek(const Date& soughtdate, const std::vector<MeteoData>& vecM, const bool& exactmatch)
{
	if (vecM.empty() || soughtdate > vecM.back().date || soughtdate < vecM.front().date) {
		//the sought date is not contained in the vector, return npos
		return npos;
	}

	const size_t max_idx = vecM.size()-1; //obviously, the index must be <= max_idx

	//since usually the sampling rate is quite constant, try to guess where our point
	//should be and provide a much smaller search interval around it
	const double start_date = vecM.front().date.getJulian(true);
	const double end_date = vecM.back().date.getJulian(true);
	const double curr_date = soughtdate.getJulian(true);
	const double raw_pos = (curr_date-start_date) / (end_date-start_date) * static_cast<double>(max_idx); //always >=0
	const size_t start_idx = (size_t)floor(raw_pos*.9);
	const size_t end_idx = MIN( (size_t)ceil(raw_pos*1.1), max_idx);

	//first and last index of the search interval, either using our initial guess or the full vector
	size_t first = (curr_date >= vecM[start_idx].date.getJulian(true))? start_idx : 0;
	size_t last = (curr_date <= vecM[end_idx].date.getJulian(true))? end_idx : max_idx;

	//if we reach this point: the date is spanned by the buffer and there are at least two elements
	if (exactmatch){
		//perform binary search
		while (first <= last) {
			const size_t mid = (first + last) / 2;  // compute mid point
			if (soughtdate > vecM[mid].date)
				first = mid + 1;                   // repeat search in top half
			else if (soughtdate < vecM[mid].date)
				last = mid - 1;                    // repeat search in bottom half
			else
				return mid;                        // found it. return position
		}
	} else {
		//perform binary search
		while (first <= last) {
			const size_t mid = (first + last) / 2;  // compute mid point

			if (mid < max_idx) {
				if ((soughtdate > vecM[mid].date) && (soughtdate < vecM[mid+1].date))
					return mid+1;
			}

			if (soughtdate > vecM[mid].date)
				first = mid + 1;                   // repeat search in top half
			else if (soughtdate < vecM[mid].date)
				last = mid - 1;                    // repeat search in bottom half
			else
				return mid;                        // found it. return position
		}
	}

	return npos;
}

std::string printFractionalDay(const double& fractional) {
	const double hours=floor(fractional*24.);
	const double minutes=floor((fractional*24.-hours)*60.);
	const double seconds=fractional*24.*3600.-hours*3600.-minutes*60.;

	std::stringstream tmp;
	tmp << std::fixed << std::setfill('0') << std::setprecision(0);
	tmp << std::setw(2) << hours << ":";
	tmp << std::setw(2) << minutes << ":";
	tmp << std::setw(2) << seconds;

	return tmp.str();
}

void getArraySliceParams(const unsigned int& dimx, const unsigned int& nbworkers, const unsigned int &wk, unsigned int& startx, unsigned int& nx)
{
	if(nbworkers>dimx) {
		std::stringstream ss;
		ss << "Can not split " << dimx << " columns in " << nbworkers << " bands!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	const unsigned int nx_avg = dimx / nbworkers;
	const unsigned int remainder = dimx % nbworkers;

	if(wk<=remainder) { //distribute remainder, 1 extra column per worker, on first workers
		nx = nx_avg+1;
		startx = (wk-1)*nx;
	} else { //all remainder has been distributed, we now attribute a normal number of columns
		nx = nx_avg;
		startx = remainder*(nx+1) + (wk-1-remainder)*nx;
	}
}

void FileIndexer::setIndex(const Date& i_date, const std::streampos& i_pos)
{
	const file_index elem(i_date, i_pos);

	//check if we can simply append the new index
	if(vecIndex.empty() || elem>vecIndex.back()) {
		vecIndex.push_back(elem);
		return;
	}

	//look for the proper position for insertion of the new index
	const std::vector< struct file_index >::iterator it = std::upper_bound(vecIndex.begin(), vecIndex.end(), elem);
	if(it>vecIndex.begin() && (it-1)->date!=elem.date) { //check that we don't try to insert a duplicate
		vecIndex.insert(it, elem); //insertion is at the proper place -> remains ordered
		return;
	}
}

void FileIndexer::setIndex(const std::string& i_date, const std::streampos& i_pos)
{
	Date tmpdate;
	convertString(tmpdate, i_date, 0.);
	setIndex(tmpdate, i_pos);
}

void FileIndexer::setIndex(const double& i_date, const std::streampos& i_pos)
{
	const Date tmpdate(i_date, 0.);
	setIndex(tmpdate, i_pos);
}

std::streampos FileIndexer::getIndex(const Date& i_date) const
{
	const size_t foundIdx = binarySearch(i_date);
	if(foundIdx==static_cast<size_t>(-1)) return static_cast<std::streampos>(-1);
	else return vecIndex[foundIdx].pos;
}

std::streampos FileIndexer::getIndex(const std::string& i_date) const
{
	Date tmpdate;
	convertString(tmpdate, i_date, 0.);
	return getIndex(tmpdate);
}

std::streampos FileIndexer::getIndex(const double& i_date) const
{
	const Date tmpdate(i_date, 0.);
	return getIndex(tmpdate);
}

size_t FileIndexer::binarySearch(const Date& soughtdate) const
{//perform binary search, return the first element that is GREATER than the provided value
	if(vecIndex.empty()) return static_cast<size_t>(-1);
	if(soughtdate<vecIndex.front().date) return static_cast<size_t>(-1);
	if(soughtdate>=vecIndex.back().date) return vecIndex.size()-1;

	const file_index elem(soughtdate, 0);
	//returns the first element that is GREATER than the provided value
	const std::vector< struct file_index >::const_iterator it = std::upper_bound(vecIndex.begin(), vecIndex.end(), elem);
	if(it>vecIndex.begin()) return it-vecIndex.begin()-1;
	else return static_cast<size_t>(-1);
}

const std::string FileIndexer::toString() const
{
	std::stringstream os;
	os << "<FileIndexer>\n";
	for(size_t ii=0; ii<vecIndex.size(); ii++)
		os << "\t" << "[" << ii << "] - " << vecIndex[ii].date.toString(Date::ISO) << " -> #" << std::hex << vecIndex[ii].pos << std::dec << "\n";
	os << "</FileIndexer>\n";
	return os.str();
}

} //namespace IOUtils
} //namespace
