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

#include <meteoio/IOUtils.h>
#include <meteoio/meteolaws/Meteoconst.h> //for math constants
#include <meteoio/Config.h>    // to avoid forward declaration hell
#include <meteoio/MeteoData.h> // to avoid forward declaration hell

#ifdef _WIN32
	#include <windows.h>
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

bool IOUtils::checkEpsilonEquality(const double& val1, const double& val2, const double& epsilon)
{
	if (((val1-epsilon) < val2) && ((val1+epsilon) > val2)) {
		return true;
	}

	return false;
}

double IOUtils::pow2(const double& val)
{
	return (val*val);
}

double IOUtils::pow3(const double& val)
{
	return (val*val*val);
}

double IOUtils::pow4(const double& val)
{
	return (val*val*val*val);
}

double IOUtils::bearing_to_angle(const double& bearing) {
	const double to_rad = Cst::PI/180.;
	return (fmod(360.-bearing+90., 360.)*to_rad);
}

double IOUtils::angle_to_bearing(const double& angle) {
	const double to_deg = 180.0 / Cst::PI;
	return (fmod( 90.-angle*to_deg+360. , 360. ));
}

double IOUtils::bearing(const std::string& bearing_str)
{
	std::string s=bearing_str;
	trim(s);
	toUpper(s);

	if(s=="N") return 0.;
	if(s=="NNE") return 22.5;
	if(s=="NE") return 45.;
	if(s=="ENE") return 67.5;
	if(s=="E") return 90.;
	if(s=="ESE") return 112.5;
	if(s=="SE") return 135.;
	if(s=="SSE") return 157.5;
	if(s=="S") return 180.;
	if(s=="SSW") return 202.5;
	if(s=="SW") return 225.;
	if(s=="WSW") return 247.5;
	if(s=="W") return 270.;
	if(s=="WNW") return 292.5;
	if(s=="NW") return 315.;
	if(s=="NNW") return 337.5;

	//no match
	return nodata;
}

void IOUtils::stripComments(std::string& str)
{
	size_t found = str.find_first_of("#;");
	if (found != std::string::npos){
		str.erase(found); //rest of line disregarded
	}
}

void IOUtils::copy_file(const std::string& src, const std::string& dest)
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

std::string IOUtils::cleanPath(const std::string& in_path) {
	std::string out_path(in_path);

	size_t curr = out_path.find('\\', 0);
	while(curr!=std::string::npos){
		out_path.replace(curr, 1, "/");
		curr = out_path.find('\\', curr);
	}
	//out_path.replace(out_path.begin(), out_path.end(), '\\', '/');
	return out_path;
}

std::string IOUtils::getExtension(const std::string& filename)
{
	const std::string whitespaces(" \t\f\v\n\r");
	const size_t start_basename = filename.find_last_of("/\\"); //we will skip the path
	const size_t startpos = filename.find_last_of('.');
	if( startpos==std::string::npos ) return std::string("");
	if( start_basename!=std::string::npos && startpos<start_basename ) return std::string("");

	const size_t endpos = filename.find_last_not_of(whitespaces); // Find the first character position from reverse af

	return filename.substr(startpos+1, endpos-startpos);
}

std::string IOUtils::removeExtension(const std::string& filename)
{
	const std::string whitespaces(" \t\f\v\n\r");
	const size_t start_basename = filename.find_last_of("/\\"); //we will skip the path
	const size_t startpos = filename.find_last_of('.');
	if( startpos==std::string::npos ) return filename;
	if( start_basename!=std::string::npos && startpos<start_basename ) return filename;

	return filename.substr(0, startpos);
}

void IOUtils::trim(std::string& str)
{
	const std::string whitespaces(" \t\f\v\n\r");
	const size_t startpos = str.find_first_not_of(whitespaces); // Find the first character position after excluding leading blank spaces
	const size_t endpos = str.find_last_not_of(whitespaces); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( std::string::npos == startpos ) || ( std::string::npos == endpos)) {
		str = "";
	} else {
		str = str.substr( startpos, endpos-startpos+1 );
	}
}

void IOUtils::toUpper(std::string& str){
	for(size_t t=0; t<str.length(); t++) {
		str[t] = (char)toupper(str[t]);
	}
}

void IOUtils::toLower(std::string& str){
	for(size_t t=0; t<str.length(); t++) {
		str[t] = (char)tolower(str[t]);
	}
}

bool IOUtils::isNumeric(const std::string& str, const unsigned int& nBase)
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs
	std::istringstream iss(s);

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

bool IOUtils::readKeyValuePair(const std::string& in_line, const std::string& delimiter,
                               std::map<std::string,std::string>& out_map, const std::string& keyprefix, const bool& setToUpperCase)
{
	//size_t pos = in_line.find(delimiter); //first occurence of '='

	size_t pos = std::string::npos;
	if ((delimiter==" ") || (delimiter=="\t")) {
		pos = in_line.find_first_of(" \t", 0);
	} else {
		pos = in_line.find(delimiter); //first occurence of '='
	}


	if(pos != std::string::npos) { //ignore in_lines that are empty or without '='
		std::string key = in_line.substr(0, pos);
		std::string value = in_line.substr(pos + 1);

		IOUtils::trim(key);
		IOUtils::trim(value);
		//cerr << "key:" << key << " val:" << value << endl;

		if ((key == "") || (value=="")) {
			return false;
		}

		if (setToUpperCase)
			IOUtils::toUpper(key);


		out_map[keyprefix + key] = value;
	} else {
		return false;
		//cerr << "line:" << in_line << "delimiter" << endl;
	}

	return true;
}

bool IOUtils::validFileName(const std::string& filename)
{
	const size_t startpos = filename.find_first_not_of(" \t\n"); // Find the first character position after excluding leading blank spaces
	if((startpos!=0) || (filename==".") || (filename=="..")) {
		return false;
	}

	return true;
}

#ifdef _WIN32
bool IOUtils::fileExists(const std::string& filename)
{
	return ( GetFileAttributes( filename.c_str() ) != INVALID_FILE_ATTRIBUTES );
}

void IOUtils::readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	WIN32_FIND_DATA ffd;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	const size_t path_length = path.length();
	if (path_length > (MAX_PATH - 3)) {
		std::cerr << "Path " << path << "is too long (" << path_length << " characters)" << std::endl;
		throw FileAccessException("Error opening directory " + path, AT);
	}

	const std::string filepath = path+"\\"+pattern;
	hFind = FindFirstFile(filepath.c_str(), &ffd);
	if (INVALID_HANDLE_VALUE == hFind) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	do {
		if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			//this is a directory -> do nothing
		} else {
			std::string filename(ffd.cFileName);
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
bool IOUtils::fileExists(const std::string& filename)
{
	struct stat buffer ;

	if ((stat( filename.c_str(), &buffer))==0) {//File exists if stat returns 0
		return true ;
	}

	return false;
}

void IOUtils::readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	DIR *dp;
	struct dirent *dirp;

	if((dp  = opendir(path.c_str())) == NULL) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	while ((dirp = readdir(dp)) != NULL) {
		const std::string tmp(dirp->d_name);
		if( tmp.compare(".")!=0 && tmp.compare("..")!=0 ) { //skip "." and ".."
			if (pattern=="") {
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

std::string IOUtils::getLogName() {
	char *tmp;
	const size_t len=64;
	char logname[len]; //to keep the user name short enough

	if((tmp=getenv("USERNAME"))==NULL) { //Windows & Unix
		if((tmp=getenv("LOGNAME"))==NULL) { //Unix
			tmp=getenv("USER"); //Windows & Unix
		}
	}
	if(tmp==NULL) return std::string("N/A");

	strncpy(logname, tmp, len-1);
	return std::string(logname);
}

void IOUtils::readKeyValueHeader(std::map<std::string,std::string>& headermap,
                                 std::istream& fin,
                                 const size_t& linecount,
                                 const std::string& delimiter)
{
	size_t linenr = 0;
	std::string line="";

	//make a test for end of line encoding:
	char eol = IOUtils::getEoln(fin);

	for (size_t ii=0; ii< linecount; ii++){
		if (std::getline(fin, line, eol)) {
			linenr++;
			bool result = IOUtils::readKeyValuePair(line, delimiter, headermap);

			if (!result) { //  means if ((key == "") || (value==""))
				std::stringstream out;
				out << "Invalid key value pair in line: " << linenr << " of header";
				throw IOException(out.str(), AT);
			}
		} else {
			throw InvalidFormatException("Premature EOF while reading Header", AT);
		}
	}
}

char IOUtils::getEoln(std::istream& fin)
{
	std::streambuf* pbuf;
	char tmp = '0';
	int chars = 0;

	const std::streampos position = fin.tellg();

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
			pbuf->pubseekpos(position); //rewind
			fin.clear(); //reset eof flag, etc
			return tmp;
		}
	} while ((chars < 3000) && (!fin.eof()));

	pbuf = fin.rdbuf();
	pbuf->pubseekpos(position); //rewind
	fin.clear(); //reset eof flag, etc

	return '\n';
}

void IOUtils::skipLines(std::istream& fin, const size_t& nbLines, const char& eoln)
{
	std::string dummy;
	for (size_t ii=0; ii<nbLines; ii++) {
		if(!getline(fin, dummy, eoln)) {
			throw InvalidFormatException("Premature EOF while skipping lines", AT);
		}
	}
}

size_t IOUtils::readLineToVec(const std::string& line_in, std::vector<double>& vec_data)
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

size_t IOUtils::readLineToVec(const std::string& line_in, std::vector<std::string>& vecString)
{
	vecString.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input

	std::string tmp_string;
	while (!iss.eof()) {
		iss >> std::skipws >> tmp_string;

		if (tmp_string != "") {
			vecString.push_back(tmp_string);
		}
		tmp_string="";
	}

	return vecString.size();
}

size_t IOUtils::readLineToVec(const std::string& line_in, std::vector<std::string>& vecString, const char& delim)
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

template<> bool IOUtils::convertString<std::string>(std::string& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	(void)f;
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs

	t = s;
	return true;
}

template<> bool IOUtils::convertString<bool>(bool& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs

	if (toupper(s[0])=='T' || toupper(s[0])=='Y' ) {
		t = true;
	} else if (toupper(s[0])=='F' || toupper(s[0])=='N' ) {
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

	std::string::size_type pos = s.find_first_not_of(ALPHANUM);
	if (pos != std::string::npos) {
		std::string tmp = s.substr(pos);
		trim(tmp);
		if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
			return false;
		}
	}

	return true;
}

template<> bool IOUtils::convertString<unsigned int>(unsigned int& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs
	if (s.size() == 0) {
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
		std::string tmp="";
		getline(iss,  tmp); //get rest of line, if any
		trim(tmp);
		if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {
			//if line holds more than one value it's invalid
			return false;
		}
		return true;
	}
}

bool IOUtils::convertString(Date& t, const std::string& str, const double& time_zone, std::ios_base& (*f)(std::ios_base&))
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs

	(void)f;
	int year;
	unsigned int month, day, hour, minute, second;
	char rest[32] = "";

	//HACK: we read the seconds, but we ignore them...
	if (sscanf(s.c_str(), "%d-%u-%u %u:%u:%u%31s", &year, &month, &day, &hour, &minute, &second, rest) >= 6) {
		t.setDate(year, month, day, hour, minute, time_zone);
	} else if (sscanf(s.c_str(), "%d-%u-%uT%u:%u:%u%31s", &year, &month, &day, &hour, &minute, &second, rest) >= 6) {
		t.setDate(year, month, day, hour, minute, time_zone);
	} else if (sscanf(s.c_str(), "%d-%u-%u %u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
		t.setDate(year, month, day, hour, minute, time_zone);
	} else if (sscanf(s.c_str(), "%d-%u-%uT%u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
		t.setDate(year, month, day, hour, minute, time_zone);
	} else if (sscanf(s.c_str(), "%d-%u-%u%31s", &year, &month, &day, rest) >= 3) {
		t.setDate(year, month, day, (unsigned)0, (unsigned)0, time_zone);
	} else if (sscanf(s.c_str(), "%u:%u%31s", &hour, &minute, rest) >= 2) {
		t.setDate( ((double)hour)/24. + ((double)minute)/24./60. , time_zone);
	} else {
		//try to read purely numerical date, potentially surrounded by other chars
		const size_t in_len = str.length();
		const size_t beg = str.find_first_of(NUM);
		if(beg==npos || beg==in_len) return false;
		size_t end = str.find_first_not_of(NUM, beg+1);
		if(end==npos) end = in_len;

		const std::string datum = str.substr(beg, end-beg);
		const size_t d_len = datum.length();
		if(d_len<10 || d_len>14) return false;
		if( convertString(year,datum.substr(0,4))==false ) return false;
		if( convertString(month,datum.substr(4,2))==false ) return false;
		if( convertString(day,datum.substr(6,2))==false ) return false;
		if( convertString(hour,datum.substr(8,2))==false ) return false;
		if(d_len==10)
			minute=0;
		else {
			if(d_len>=12) {
				if( convertString(minute,datum.substr(10,2))==false ) return false;
			} else
				return false;
			if(d_len==12)
				second=0;
			else {
				if(d_len==14) {
					if( convertString(second,datum.substr(12,2))==false ) return false;
				} else
					return false;
			}
		}

		t.setDate( year, month, day, hour, minute, time_zone );
	}

	std::string tmp(rest);
	trim(tmp);
	if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
		return false;
	}

	return true;
}

template<> bool IOUtils::convertString<Coords>(Coords& t, const std::string& str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs

	(void)f;
	double lat, lon;
	try {
		Coords::parseLatLon(s,lat, lon);
	} catch(const IOException&) {
		return false;
	}
	t.setLatLon(lat, lon, nodata);

	return true;
}


void IOUtils::getProjectionParameters(const Config& cfg, std::string& coordin, std::string& coordinparam,
                                      std::string& coordout, std::string& coordoutparam) {
	cfg.getValue("COORDSYS", "Input", coordin);
	cfg.getValue("COORDPARAM", "Input", coordinparam, Config::nothrow);
	cfg.getValue("COORDSYS", "Output", coordout, Config::nothrow);
	cfg.getValue("COORDPARAM", "Output", coordoutparam, Config::nothrow);
}

void IOUtils::getTimeZoneParameters(const Config& cfg, double& tz_in, double& tz_out) {
	cfg.getValue("TIME_ZONE", "Input", tz_in, Config::nothrow);
	cfg.getValue("TIME_ZONE", "Output", tz_out, Config::nothrow);
}

size_t IOUtils::seek(const Date& soughtdate, const std::vector<MeteoData>& vecM, const bool& exactmatch){
	//returns index of element, if element does not exist it returns closest index after soughtdate
	//the element needs to be an exact hit or embedded between two measurments

	if (vecM.size() <= 0) {//no elements in buffer
		return IOUtils::npos;
	}

	//if we reach this point: at least one element in buffer
	if (vecM[0].date > soughtdate) {
		return IOUtils::npos;
	}

	if (vecM[vecM.size()-1].date < soughtdate) {//last element is earlier, return npos
		return IOUtils::npos;
	}

	if (vecM[0].date == soughtdate) {//closest element
		return 0;
	}

	//if we reach this point: the date is spanned by the buffer and there are at least two elements
	//HACK: would it be better to create a timeseries object and call vector's binary search on it?
	if (exactmatch){
		size_t first = 1, last = vecM.size()-1;

		//perform binary search
		while (first <= last) {
			size_t mid = (first + last) / 2;  // compute mid point
			if (soughtdate > vecM[mid].date)
				first = mid + 1;                   // repeat search in top half
			else if (soughtdate < vecM[mid].date)
				last = mid - 1;                    // repeat search in bottom half
			else
				return mid;                        // found it. return position
		}
	} else {
		size_t first = 0, last = vecM.size()-1;

		//perform binary search
		while (first <= last) {
			size_t mid = (first + last) / 2;  // compute mid point

			if (mid < (vecM.size()-1))
				if ((soughtdate > vecM[mid].date) && (soughtdate < vecM[mid+1].date))
					return mid+1;

			if (soughtdate > vecM[mid].date)
				first = mid + 1;                   // repeat search in top half
			else if (soughtdate < vecM[mid].date)
				last = mid - 1;                    // repeat search in bottom half
			else
				return mid;                        // found it. return position

		}
	}

	return IOUtils::npos;
}

std::string IOUtils::printFractionalDay(const double& fractional) {
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

void IOUtils::getArraySliceParams(const unsigned int& dimx, const unsigned int& nbworkers, const unsigned int &wk, unsigned int& startx, unsigned int& nx)
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

} //namespace
