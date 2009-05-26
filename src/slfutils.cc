#include "slfutils.h"

bool slfutils::checkEpsilonEquality(double val1, double val2, double epsilon){
  if (((val1-epsilon) < val2) && ((val1+epsilon) > val2))
    return true;

  return false;
}

void slfutils::WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out){
//converts WGS84 coordinates (lat,long) to the Swiss coordinates. See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf
//The elevation is supposed to be above sea level, so it does not require any conversion
//lat and long must be decimal (and they will be converted to seconds)
	const double phi_p = (lat_in*3600. - 169028.66) / 10000.;
	const double lambda_p = (long_in*3600. - 26782.5) / 10000.;

	east_out = 600072.37
		+ 211455.93	* lambda_p
		- 10938.51	* lambda_p * phi_p
		- 0.36		* lambda_p * (phi_p*phi_p)
		- 44.54		* (lambda_p*lambda_p*lambda_p);

	north_out = 200147.07
		+ 308807.95	* phi_p
		+ 3745.25	* (lambda_p*lambda_p)
		+ 76.63		* (phi_p*phi_p)
		- 194.56	* (lambda_p*lambda_p) * phi_p
		+ 119.79	* (phi_p*phi_p*phi_p);

	/*// if necessary for the elevation, uncomment this block
	h_out = h_in - 49.55
		+ 2.73		* lambda_p
		+ 6.94		* phi_p;
	*/
}

void slfutils::CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out){
//converts Swiss coordinates to WGS84 coordinates (lat,long). See http://geomatics.ladetto.ch/ch1903_wgs84_de.pdf
//The elevation is supposed to be above sea level, so it does not require any conversion
//lat and long are decimal
	const double y_p = (east_in - 600000.) / 1000000.;
	const double x_p = (north_in - 200000.) / 1000000.;

	const double lambda_p = 2.6779094
		+ 4.728982	* y_p
		+ 0.791484	* y_p * x_p
		+ 0.1306	* y_p * (x_p*x_p)
		- 0.0436	* (y_p*y_p*y_p);

	const double phi_p = 16.9023892
		+ 3.238272	* x_p
		- 0.270978	* (y_p*y_p)
		- 0.002528	* (x_p*x_p)
		- 0.0447	* (y_p*y_p) * x_p
		- 0.0140	* (x_p*x_p*x_p);

	lat_out = phi_p * 100./36.;
	long_out = lambda_p * 100./36.;

	/*// if necessary for the elevation, uncomment this block
	h_out = h_in + 49.55
		- 12.60		* y_p
		- 6.94		* x_p;
	*/
}


void slfutils::trim(string& str){
  size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces  
  size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af  

  // if all spaces or empty return an empty string  
  if(( string::npos == startpos ) || ( string::npos == endpos)){  
    str = "";  
  } else  
    str = str.substr( startpos, endpos-startpos+1 );    
}

/**
  FUNCTION readKeyValuePair(const string& in_line, const string& delimiter, map<string,string>& out_map)
  - read a string line, parse it and save it into a map object, that is passed by reference
  - delimiter: string that separates key and value
  - return value: bool, depending on success of operation. true when line is empty
 */
bool slfutils::readKeyValuePair(const string& in_line, const string& delimiter, map<string,string>& out_map){

  //size_t pos = in_line.find(delimiter); //first occurence of '='

  size_t pos = string::npos;
  if ((delimiter==" ") || (delimiter=="\t")){
    pos = in_line.find_first_of(" \t", 0);
  } else {
    pos = in_line.find(delimiter); //first occurence of '='    
  }


  if(pos != string::npos) { //ignore in_lines that are empty or without '='
      string key = in_line.substr(0, pos);
      string value = in_line.substr(pos + 1);

      slfutils::trim(key);
      slfutils::trim(value);
      //cerr << "key:" << key << " val:" << value << endl;

      if ((key == "") || (value==""))
	return false;

      out_map[key] = value;
  } else {
    //cerr << "line:" << in_line << "delimiter" << endl;
  }
  
  return true;
}

bool slfutils::fileExists(const std::string& filename) {
  struct stat buffer ;

  if ((stat( filename.c_str(), &buffer))==0) //File exists if stat returns 0
    return true ;

  return false;
}

bool slfutils::validFileName(const std::string& filename) {
  size_t startpos = filename.find_first_not_of(" \t\n"); // Find the first character position after excluding leading blank spaces  
  if((startpos!=0) || (filename==".") || (filename==".."))
    return false;

  return true;
}

void slfutils::readKeyValueHeader(map<string,string>& headermap, 
				  std::istream& fin,
				  const unsigned int& linecount, 
				  const std::string& delimiter)
{
  int linenr = 0;
  std::string line="";

  //make a test for end of line encoding:
  char eol = slfutils::getEoln(fin);

  for (unsigned int ii=0; ii< linecount; ii++){
    if (std::getline(fin, line, eol)) {
      //cout << line <<endl;
      linenr++;

      bool result = slfutils::readKeyValuePair(line, delimiter, headermap);

      if (!result) { //  means if ((key == "") || (value==""))
	stringstream out;
	out << "Invalid key value pair in line: " << linenr << " of header";
       THROW SLFException(out.str(), AT);
      }
    } else {
      THROW InvalidFormatException("Premature EOF while reading Header", AT);
    }
  }
}

char slfutils::getEoln(std::istream& fin){
  streambuf* pbuf;
  char tmp = '0';
  int chars = 0;

  int position = fin.tellg();

  do {
    fin.get(tmp);
    chars++;

    if ((tmp == '\r') || (tmp == '\n')){
      char peekc = tmp;
      //cout << (int)tmp << endl;
      while ((!fin.eof() && ((peekc=='\r') || (peekc=='\n')))){
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

void slfutils::skipLines(std::istream& fin, unsigned int nbLines, char eoln){
  string dummy;
  for (unsigned int ii=0; ii<nbLines; ii++){
    if(!getline(fin, dummy, eoln))
      THROW InvalidFormatException("Premature EOF while skipping lines", AT);
  }
}

unsigned int slfutils::readLineToVec(const string& line_in, vector<string>& vecString){
  vecString.clear();
  std::istringstream iss(line_in); //construct inputstream with the string line as input

  string tmp_string;
  while (!iss.eof()){
    iss >> skipws >> tmp_string;
    
    if (tmp_string != ""){
      vecString.push_back(tmp_string);
    }
    tmp_string="";
  }

  return vecString.size();
}

void slfutils::readDirectory(const string& path, list<string>& dirlist, const string& pattern){
  DIR *dp;
  struct dirent *dirp;
  
  if((dp  = opendir(path.c_str())) == NULL) {
    THROW FileAccessException("Error opening directory " + path, AT);
  }

  while ((dirp = readdir(dp)) != NULL) {
    string tmp = string(dirp->d_name);
    if (pattern=="") {
      dirlist.push_back(tmp);
    } else {
      size_t pos = tmp.find(pattern);
      if (pos!=string::npos){
	dirlist.push_back(tmp);
      }
    }
  }
  closedir(dp);
}

// generic template function convertString must be defined in the header

const char ALPHANUM[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
template<> bool slfutils::convertString<bool>(bool& t, const std::string str, std::ios_base& (*f)(std::ios_base&)) {
  string s = str; 
  trim(s); //delete trailing and leading whitespaces and tabs
  
  if (toupper(s[0])=='T' ||
      toupper(s[0])=='Y' ) {
    t = true;
  } else
  if (toupper(s[0])=='F' ||
      toupper(s[0])=='N' ) {
    t = false;
  } else {
    std::istringstream iss(s);
    int i;
    iss >> f >> i; //Convert first part of stream with the formatter (e.g. std::dec, std::oct)
    if (iss.fail()) //Conversion failed
      return false;
    t = (i != 0);
  }
  
  unsigned int pos = s.find_first_not_of(ALPHANUM);
  if (pos != string::npos) {
    string tmp = s.substr(pos);
    trim(tmp);
    if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';')//if line holds more than one value it's invalid
      return false;
  }
  
  return true;
}

template<> bool slfutils::convertString<Date>(Date& t, const std::string str, std::ios_base& (*f)(std::ios_base&)) {
  string s = str; 
  trim(s); //delete trailing and leading whitespaces and tabs
  
  (void)f;
  unsigned int year, month, day, hour, minute;
  char rest[32] = "";
  if (sscanf(s.c_str(), "%u-%u-%u %u:%u%31s", &year, &month, &day, &hour, &minute, rest) >= 5) {
    t.setDate(year, month, day, hour, minute);
  } else if (sscanf(s.c_str(), "%u-%u-%u%31s", &year, &month, &day, rest) >= 3) {
    t.setDate(year, month, day);
  } else if (sscanf(s.c_str(), "%u:%u%31s", &hour, &minute, rest) >= 2) {
    t = Date( ((double)hour)/24. + ((double)minute)/24./60. );
  } else {
    return false;
  }
  
  string tmp(rest);
  trim(tmp);
  if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';')//if line holds more than one value it's invalid
    return false;

  return true;
}
