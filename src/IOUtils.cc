#include "IOUtils.h"

bool IOUtils::checkEpsilonEquality(double val1, double val2, double epsilon)
{
	if (((val1-epsilon) < val2) && ((val1+epsilon) > val2)) {
		return true;
	}

	return false;
}

double IOUtils::pow2(const double val)
{
	return (val*val);
}

void IOUtils::trim(std::string& str)
{
	const std::string whitespaces (" \t\f\v\n\r");
	size_t startpos = str.find_first_not_of(whitespaces); // Find the first character position after excluding leading blank spaces  
	size_t endpos = str.find_last_not_of(whitespaces); // Find the first character position from reverse af  

	// if all spaces or empty return an empty string  
	if(( string::npos == startpos ) || ( string::npos == endpos)) {  
		str = "";  
	} else { 
		str = str.substr( startpos, endpos-startpos+1 );    
	}
}

void IOUtils::toUpper(std::string& str){
	for(unsigned int t=0; t<str.length(); t++) {
		str[t] = toupper(str[t]);
	}
}

bool IOUtils::readKeyValuePair(const std::string& in_line, const std::string& delimiter,
				std::map<std::string,std::string>& out_map, const std::string& keyprefix)
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

		out_map[keyprefix + key] = value;
	} else {
		return false;
		//cerr << "line:" << in_line << "delimiter" << endl;
	}

	return true;
}

bool IOUtils::fileExists(const std::string& filename)
{
	struct stat buffer ;

	if ((stat( filename.c_str(), &buffer))==0) {//File exists if stat returns 0
		return true ;
	}

	return false;
}

bool IOUtils::validFileName(const std::string& filename)
{
	size_t startpos = filename.find_first_not_of(" \t\n"); // Find the first character position after excluding leading blank spaces  
	if((startpos!=0) || (filename==".") || (filename=="..")) {
		return false;
	}

	return true;
}

void IOUtils::readKeyValueHeader(std::map<std::string,std::string>& headermap, 
				  std::istream& fin,
				  const unsigned int& linecount, 
				  const std::string& delimiter)
{
	int linenr = 0;
	std::string line="";

	//make a test for end of line encoding:
	char eol = IOUtils::getEoln(fin);

	for (unsigned int ii=0; ii< linecount; ii++){
		if (std::getline(fin, line, eol)) {
			//cout << line <<endl;
			linenr++;

			bool result = IOUtils::readKeyValuePair(line, delimiter, headermap);

			if (!result) { //  means if ((key == "") || (value==""))
				stringstream out;
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

	int position = fin.tellg();

	do {
		fin.get(tmp);
		chars++;

		if ((tmp == '\r') || (tmp == '\n')) {
			char peekc = tmp;
			//cout << (int)tmp << endl;
			while ((!fin.eof() && ((peekc=='\r') || (peekc=='\n')))) {
				tmp = peekc;
				fin.get(peekc);
				chars++;
			}
			//      cout << "Read: " << chars << "characters" << endl;
			pbuf = fin.rdbuf();
			pbuf->pubseekpos(position); //rewind
			//cout << fin.peek();
			return tmp;
		}
	} while ((chars < 3000) && (!fin.eof()));
	   
	pbuf = fin.rdbuf();
	pbuf->pubseekpos(position); //rewind

	return '\n';
}

void IOUtils::skipLines(std::istream& fin, unsigned int nbLines, char eoln)
{
	std::string dummy;
	for (unsigned int ii=0; ii<nbLines; ii++) {
		if(!getline(fin, dummy, eoln)) {
			throw InvalidFormatException("Premature EOF while skipping lines", AT);
		}
	}
}

unsigned int IOUtils::readLineToVec(const std::string& line_in, std::vector<std::string>& vecString)
{
	vecString.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input

	std::string tmp_string;
	while (!iss.eof()) {
		iss >> skipws >> tmp_string;

		if (tmp_string != "") {
			vecString.push_back(tmp_string);
		}
		tmp_string="";
	}

	return vecString.size();
}

unsigned int IOUtils::readLineToVec(const std::string& line_in, std::vector<std::string>& vecString, const char& delim)
{
	vecString.clear();
	std::string tmp_string;
	std::istringstream iss(line_in);

	while (getline(iss, tmp_string, delim)){
		vecString.push_back(tmp_string);
	}

	return vecString.size();
}


void IOUtils::readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	DIR *dp;
	struct dirent *dirp;

	if((dp  = opendir(path.c_str())) == NULL) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	while ((dirp = readdir(dp)) != NULL) {
		string tmp = string(dirp->d_name);
		if (pattern=="") {
			dirlist.push_back(tmp);
		} else {
			size_t pos = tmp.find(pattern);
			if (pos!=string::npos) {
				dirlist.push_back(tmp);
			}
		}
	}
	closedir(dp);
}

// generic template function convertString must be defined in the header

const char ALPHANUM[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

template<> bool IOUtils::convertString<std::string>(std::string& t, const std::string str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs

	t = s;
	return true;
}

template<> bool IOUtils::convertString<bool>(bool& t, const std::string str, std::ios_base& (*f)(std::ios_base&))
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
	if (pos != string::npos) {
		std::string tmp = s.substr(pos);
		trim(tmp);
		if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
			return false;
		}
	}

	return true;
}

template<> bool IOUtils::convertString<Date_IO>(Date_IO& t, const std::string str, std::ios_base& (*f)(std::ios_base&))
{
	std::string s = str;
	trim(s); //delete trailing and leading whitespaces and tabs

	(void)f;
	unsigned int year, month, day, hour, minute;
	char rest[32] = "";
	if (sscanf(s.c_str(), "%u-%u-%u %u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
		t.setDate(year, month, day, hour, minute);
	} else if (sscanf(s.c_str(), "%u-%u-%uT%u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
		t.setDate(year, month, day, hour, minute);
	} else if (sscanf(s.c_str(), "%u-%u-%u%31s", &year, &month, &day, rest) >= 3) {
		t.setDate(year, month, day);
	} else if (sscanf(s.c_str(), "%u:%u%31s", &hour, &minute, rest) >= 2) {
		t = Date_IO( ((double)hour)/24. + ((double)minute)/24./60. );
	} else {
		return false;
	}

	std::string tmp(rest);
	trim(tmp);
	if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
		return false;
	}

	return true;
}
