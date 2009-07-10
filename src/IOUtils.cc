#include "IOUtils.h"

#ifndef PI
	#define PI 3.141592653589
#endif

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

double IOUtils::normalizeBearing(double angle)
{
	if(angle<0.) angle = 360.0 + angle;
	if(angle>360.) angle = angle - 360.;
	return angle;
}

void IOUtils::WGS84_to_CH1903(const double& lat_in, const double& long_in, double& east_out, double& north_out)
{
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

void IOUtils::CH1903_to_WGS84(const double& east_in, const double& north_in, double& lat_out, double& long_out)
{
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
		- 22.64		* x_p;
	*/
}

void IOUtils::WGS84_to_local(const double& lat_ref, const double& lon_ref, const double& lat, const double& lon, double& easting, double& northing)
{
//	easting = VincentyDistance((lat_ref+lat)/2., lon_ref, (lat_ref+lat)/2., lon);
//	northing = VincentyDistance(lat_ref, (lon_ref+lon)/2., lat, (lon_ref+lon)/2.);
//	easting = cosineDistance((lat_ref+lat)/2., lon_ref, (lat_ref+lat)/2., lon);
//	northing = cosineDistance(lat_ref, (lon_ref+lon)/2., lat, (lon_ref+lon)/2.);
	double alpha;
	const double to_rad = PI / 180.0;
	const double distance = VincentyDistance(lat_ref, lon_ref, lat, lon, alpha);
	easting = -distance*sin(alpha*to_rad);
	northing = distance*cos(alpha*to_rad);
}

void IOUtils::local_to_WGS84(const double& lat_ref, const double& lon_ref, const double& easting, const double& northing, double& lat, double& lon)
{
	const double to_deg = 180.0 / PI;
	double bearing = atan2(northing, easting)*to_deg;
	const double distance = sqrt( pow2(easting) + pow2(northing) );

	bearing = normalizeBearing(bearing);
	VincentyInverse(lat_ref, lon_ref, distance, bearing, lat, lon);
}

double IOUtils::cosineDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2)
{
	const double Rearth = 6371.e3;
	const double to_rad = PI / 180.0;
	const double d = acos( 
		sin(lat1*to_rad) * sin(lat2*to_rad) 
		+ cos(lat1*to_rad) * cos(lat2*to_rad) * cos((lon2-lon1)*to_rad) 
		) * Rearth;
	return d;
}

double IOUtils::VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2)
{
	double alpha;
	return VincentyDistance(lat1, lon1, lat2, lon2, alpha);
}

double IOUtils::VincentyDistance(const double& lat1, const double& lon1, const double& lat2, const double& lon2, double& alpha)
{
	const double thresh = 1.e-12;	//convergence absolute threshold
	const int n_max = 100;		//maximum number of iterations
	const double a = 6378137.;	//major ellipsoid semi-axis, value for wgs84
	const double b = 6356752.3142;	//minor ellipsoid semi-axis, value for wgs84
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = PI / 180.0;
	
	const double L = (lon1 - lon2)*to_rad;
	const double U1 = atan( (1.-f)*tan(lat1*to_rad) );
	const double U2 = atan( (1.-f)*tan(lat2*to_rad) );

	double lambda = L, lambda_p=0., delta_sigma;
	double sin_sigma, cos_sigma, sigma, sin_alpha, cos_alpha2, cos_2sigma_m;
	double C, u2, A, B, s;
	int n=0;
	do {
		sin_sigma = sqrt( pow2(cos(U2)*sin(lambda)) + pow2(cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda)) );
		if(sin_sigma==0.) {
			//co-incident points
			return 0.;
		}
		cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
		sigma = atan2(sin_sigma,cos_sigma);
		sin_alpha = cos(U1)*cos(U2)*sin(lambda) / sin_sigma;
		cos_alpha2 = 1. - pow2(sin_alpha);
		if(lat1==0. && lat2==0.) {
			cos_2sigma_m = 0.;
		} else {
			cos_2sigma_m = cos_sigma - 2.*sin(U1)*sin(U2)/cos_alpha2;
		}
		C = f/16. * cos_alpha2*(4.+f*(4.-3.*cos_alpha2));
		lambda_p = lambda;
		lambda = L + (1.-C)*f*sin_alpha*( 
			sigma + C*sin_sigma*( cos_2sigma_m + C * cos_sigma * (-1.+2.*pow2(cos_2sigma_m)) ) 
			);
		n++;
	} while ( (n<n_max) && (fabs(lambda - lambda_p) > thresh) );
	
	if(n>n_max) {
		throw IOException("Distance calculation not converging", AT);
	}
	
	u2 = cos_alpha2 * (a*a - b*b) / (b*b);
	A = 1. + u2/16384. * ( 4096.+u2*(-768.+u2*(320.-175.*u2)) );
	B = u2/1024. * ( 256.+u2*(-128.+u2*(74.-47.*u2)) );
	delta_sigma = B*sin_sigma*( cos_2sigma_m+B/4.*( cos_sigma*(-1.+2.*pow2(cos_2sigma_m)) - B/6.*(cos_2sigma_m*(-3.+4.*pow2(sin_sigma))*(-3.+4.*pow2(cos_2sigma_m))) ) );

	s = b*A*(sigma - delta_sigma);	//distance between the two points
	
	//computation of the average forward bearing
	double alpha1 = atan2(cos(U2)*sin(lambda), cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda)) / to_rad; //forward azimuth
	double alpha2 = atan2(cos(U1)*sin(lambda), sin(U1)*cos(U2)-cos(U1)*sin(U2)*cos(lambda)) / to_rad; //reverse azimuth
	//converting reverse azimuth back to normal azimuth
	if(alpha2>180.) {
		alpha2 = alpha2 - 180.;
	} else {
		alpha2 = 180. - alpha2;
	}
	alpha1 = normalizeBearing(alpha1);
	alpha2 = normalizeBearing(alpha2);

	alpha = (alpha1+alpha2)/2.;
	return s;
}

void IOUtils::VincentyInverse(const double& lat_ref, const double& lon_ref, const double& distance, const double& bearing, double& lat, double& lon)
{
	const double thresh = 1.e-12;	//convergence absolute threshold
	const double a = 6378137.;	//major ellipsoid semi-axis, value for wgs84
	const double b = 6356752.3142;	//minor ellipsoid semi-axis, value for wgs84
	const double f = (a - b) / a;	//ellispoid flattening
	const double to_rad = PI / 180.0;
	const double to_deg = 180.0 / PI;

	const double alpha1 = bearing*to_rad;
	const double tanU1 = (1.-f)*tan(lat_ref*to_rad);
	const double cosU1 = 1./sqrt(1.+pow2(tanU1));
	const double sinU1 = tanU1*cosU1;
	const double sigma1 = atan2(tanU1,cos(alpha1));
	const double sinAlpha = cosU1*sin(alpha1);
	const double cos2alpha = 1. - pow2(sinAlpha);
	const double u2 = cos2alpha * (a*a - b*b) / (b*b);
	const double A = 1. + u2/16384. * (4096. + u2*(-768.+u2*(320.-175*u2)) );
	const double B = u2/1024. * (256. + u2*(-128.+u2*(74.-47.*u2)));
	
	double sigma = distance / (b*A);
	double sigma_p = 2.*PI;
	double cos2sigma_m;

	while (fabs(sigma - sigma_p) > thresh) {
		cos2sigma_m = cos( 2.*sigma1 + sigma );
		double delta_sigma = B*sin(sigma) * ( cos2sigma_m + B/4. * ( 
			cos(sigma)*(-1.+2.*pow2(cos2sigma_m)) 
			-B/6. * cos2sigma_m * (-3.+4.*pow2(sin(sigma))) * (-3.+4.*pow2(cos2sigma_m))
			) );
		sigma_p = sigma;
		sigma = distance / (b*A) + delta_sigma;
	}
	
	lat = atan2( sinU1*cos(sigma) + cosU1*cos(sigma)*cos(alpha1),
		     (1.-f) * sqrt( pow2(sinAlpha) + pow2(sinU1*sin(sigma) - cosU1*cos(sigma)*cos(alpha1)) )
		   );
	const double lambda = atan2( sin(sigma)*sin(alpha1), cosU1*cos(sigma) - sinU1*sin(sigma)*cos(alpha1) );
	const double C = f/16. * cos2alpha * (4.+f*(4.-3.*cos2alpha));
	const double L = lambda - (1.-C) * f * sinAlpha * (
				sigma + C * sin(sigma) * ( cos2sigma_m+C*cos(sigma) * (-1.+2.*pow2(cos2sigma_m)) )
				);

	lat = lat * to_deg;
	lon = lon_ref + (L*to_deg);
	//const double alpha2 = atan2( sinAlpha, -(sinU1*sin(sigma)-cosU1*cos(sigma)*cos(alpha1)) ); //reverse azimuth
}

void IOUtils::trim(string& str)
{
	const string whitespaces (" \t\f\v\n\r");
	size_t startpos = str.find_first_not_of(whitespaces); // Find the first character position after excluding leading blank spaces  
	size_t endpos = str.find_last_not_of(whitespaces); // Find the first character position from reverse af  

	// if all spaces or empty return an empty string  
	if(( string::npos == startpos ) || ( string::npos == endpos)) {  
		str = "";  
	} else { 
		str = str.substr( startpos, endpos-startpos+1 );    
	}
}

bool IOUtils::readKeyValuePair(const string& in_line, const string& delimiter, map<string,string>& out_map)
{
	//size_t pos = in_line.find(delimiter); //first occurence of '='

	size_t pos = string::npos;
	if ((delimiter==" ") || (delimiter=="\t")) {
		pos = in_line.find_first_of(" \t", 0);
	} else {
		pos = in_line.find(delimiter); //first occurence of '='    
	}


	if(pos != string::npos) { //ignore in_lines that are empty or without '='
		string key = in_line.substr(0, pos);
		string value = in_line.substr(pos + 1);

		IOUtils::trim(key);
		IOUtils::trim(value);
		//cerr << "key:" << key << " val:" << value << endl;

		if ((key == "") || (value=="")) {
			return false;
		}

		out_map[key] = value;
	} else {
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

void IOUtils::readKeyValueHeader(map<string,string>& headermap, 
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
	streambuf* pbuf;
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
	string dummy;
	for (unsigned int ii=0; ii<nbLines; ii++) {
		if(!getline(fin, dummy, eoln)) {
			throw InvalidFormatException("Premature EOF while skipping lines", AT);
		}
	}
}

unsigned int IOUtils::readLineToVec(const string& line_in, vector<string>& vecString)
{
	vecString.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input

	string tmp_string;
	while (!iss.eof()) {
		iss >> skipws >> tmp_string;

		if (tmp_string != "") {
			vecString.push_back(tmp_string);
		}
		tmp_string="";
	}

	return vecString.size();
}

void IOUtils::readDirectory(const string& path, list<string>& dirlist, const string& pattern)
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

template<> bool IOUtils::convertString<string>(string& t, const std::string str, std::ios_base& (*f)(std::ios_base&))
{
	string s = str; 
	trim(s); //delete trailing and leading whitespaces and tabs
	if (s.size() == 0) {
		t = "";
		return true;
	} else {
		std::istringstream iss(s);
		iss >> f >> t; //Convert first part of stream with the formatter (e.g. std::dec, std::oct)
	
		if (iss.fail()) {//Conversion failed
			return false;
		}
		string tmp="";
		getline(iss,  tmp); //get rest of line, if any
		trim(tmp);
		if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
			return false;
		}
		return true;
	}
}

template<> bool IOUtils::convertString<bool>(bool& t, const std::string str, std::ios_base& (*f)(std::ios_base&))
{
	string s = str; 
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
		string tmp = s.substr(pos);
		trim(tmp);
		if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
			return false;
		}
	}

	return true;
}

template<> bool IOUtils::convertString<Date_IO>(Date_IO& t, const std::string str, std::ios_base& (*f)(std::ios_base&))
{
	string s = str; 
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

	string tmp(rest);
	trim(tmp);
	if ((tmp.length() > 0) && tmp[0] != '#' && tmp[0] != ';') {//if line holds more than one value it's invalid
		return false;
	}

	return true;
}
