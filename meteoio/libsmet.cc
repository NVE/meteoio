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
#include <meteoio/libsmet.h>

using namespace std;

namespace smet {

const std::string SMETCommon::smet_version = "1.1";
set<string> SMETCommon::all_mandatory_header_keys = set<std::string>();
set<string> SMETCommon::all_optional_header_keys  = set<std::string>();
set<string> SMETCommon::all_decimal_header_values = set<std::string>();

SMETException::SMETException(const std::string& message, const std::string& position)
{
	if (position=="") {
		msg = "At unknown position: " + message;
	} else {
		msg = position + ": " + message;
	}
}

SMETException::~SMETException() throw() {}

const char* SMETException::what() const throw()
{
	return msg.c_str();
}

const bool SMETCommon::__init = SMETCommon::initStaticData();
bool SMETCommon::initStaticData()
{
	all_mandatory_header_keys.insert("station_id");
	all_mandatory_header_keys.insert("nodata");
	all_mandatory_header_keys.insert("fields");

	all_optional_header_keys.insert("latitude");
	all_optional_header_keys.insert("longitude");
	all_optional_header_keys.insert("altitude");
	all_optional_header_keys.insert("easting");
	all_optional_header_keys.insert("northing");
	all_optional_header_keys.insert("epsg");

	all_optional_header_keys.insert("station_name");
	all_optional_header_keys.insert("tz");
	all_optional_header_keys.insert("creation");
	all_optional_header_keys.insert("source");
	all_optional_header_keys.insert("units_offset");
	all_optional_header_keys.insert("units_multiplier");
	all_optional_header_keys.insert("comment");

	all_decimal_header_values.insert("nodata");
	all_decimal_header_values.insert("latitude");
	all_decimal_header_values.insert("longitude");
	all_decimal_header_values.insert("altitude");
	all_decimal_header_values.insert("easting");
	all_decimal_header_values.insert("northing");
	all_decimal_header_values.insert("epsg");
	all_decimal_header_values.insert("tz");

	return true;
}

double SMETCommon::convert_to_double(const std::string& in_string)
{
	istringstream ss(in_string);
	double value;
	if (!(ss >> value)) throw SMETException("Value '" + in_string + "' cannot be converted to double", SMET_AT);

	return value;
}

void SMETCommon::trim(std::string& str)
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

bool SMETCommon::readKeyValuePair(const std::string& in_line, const std::string& delimiter,
						    std::map<std::string,std::string>& out_map)
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

		SMETCommon::trim(key);
		SMETCommon::trim(value);
		//cerr << "key:" << key << " val:" << value << endl;

		if ((key == "") || (value=="")) {
			return false;
		}

		out_map[key] = value;
	} else {
		return false;
		//cerr << "line:" << in_line << "delimiter" << endl;
	}

	return true;
}

void SMETCommon::stripComments(std::string& str)
{
	size_t found = str.find_first_of("#;");
	if (found != std::string::npos){
		str.erase(found); //rest of line disregarded
	}
}

void SMETCommon::toUpper(std::string& str)
{
	for(size_t t=0; t<str.length(); t++) {
		str[t] = (char)toupper(str[t]);
	}
}

char SMETCommon::getEoln(std::istream& fin)
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
			//cout << (int)tmp << endl;
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

bool SMETCommon::is_decimal(const std::string& value)
{
	istringstream ss(value);
	double doublevalue;
	if (!(ss >> doublevalue)) return false;

	return true;
}

size_t SMETCommon::readLineToVec(const std::string& line_in, std::vector<std::string>& vec_string)
{
	vec_string.clear();
	std::istringstream iss(line_in); //construct inputstream with the string line as input

	std::string tmp_string;
	while (!iss.eof()) {
		iss >> std::skipws >> tmp_string;

		if (tmp_string != "") {
			vec_string.push_back(tmp_string);
		}
		tmp_string="";
	}

	return vec_string.size();
}

SMETWriter::SMETWriter(const std::string& in_filename, const SMETType& in_type, const bool& in_gzip) 
	: filename(in_filename), smet_type(in_type), gzip(in_gzip), nr_of_fields(0), 
	  julian_field(0), timestamp_field(0), nodata_value(-999.), nodata_string(""), 
	  location_in_header(false), location_in_data_wgs84(false), location_in_data_epsg(false),
	  timestamp_present(false), julian_present(false), file_is_binary(false),
	  location_wgs84(0), location_epsg(0)
{

}

SMETWriter::~SMETWriter()
{
	cleanup();
}

void SMETWriter::cleanup() throw()
{
	//clear ios flags
	fout << resetiosflags(ios_base::fixed | ios_base::left);

	if (fout.is_open()) //close fout if open
		fout.close();
}

void SMETWriter::set_header_value(const std::string& key, const double& value)
{
	//check if key is decimal, transform to string and add to header
	if (SMETCommon::all_decimal_header_values.find(key) != SMETCommon::all_decimal_header_values.end()){
		stringstream ss;
		if ((key == "latitude") || (key == "longitude") || (key == "easting") || (key == "northing")){
			ss << fixed << setprecision(6) << value;
		} else if (key == "altitude"){
			ss << fixed << setprecision(1) << value;
		} else if ((key == "epsg") || (key == "tz")){
			ss << fixed << setprecision(0) << value;
		} else {
			ss << value; //for nodata
		}

		set_header_value(key, ss.str());
	} else {
		throw SMETException("Trying to set a decimal value when a non-decimal is expected", SMET_AT);
	}
}

void SMETWriter::set_header_value(const std::string& key, const std::string& value)
{
	//check if header key/value pair is valid
	if (valid_header_pair(key, value)){
		header[key] = value;		
	} else {
		throw SMETException("Invalid, inconsistent or unknown key/value pair: " + key + " = " + value, SMET_AT);
	}
}

bool SMETWriter::valid_header_pair(const std::string& key, const std::string& value)
{
	bool key_ok = false;

	if (SMETCommon::all_mandatory_header_keys.find(key) != SMETCommon::all_mandatory_header_keys.end()){
		mandatory_header_keys.insert(key);
		key_ok = true;
	}

	//if (SMETCommon::all_optional_header_keys.find(key) != SMETCommon::all_optional_header_keys.end())
	key_ok = true; //it doesn't matter if key is mandatory or optional

	//nodata value needs extra treatment
	if (key == "nodata"){
		istringstream ss(value);
		if (!(ss >> nodata_value)) return false;
		nodata_string = value;
	}

	//Using WGS_84 coordinate system
	if (key == "latitude")  location_wgs84 |= 1;
	else if (key == "longitude") location_wgs84 |= 2;
	else if (key == "altitude")  location_wgs84 |= 4;

	//Using an EPSG coordinate system
	if (key == "easting")  location_epsg |= 1;
	else if (key == "northing") location_epsg |= 2;
	else if (key == "altitude") location_epsg |= 4;
	else if (key == "epsg")     location_epsg |= 8;

	//Now do some value checks
	if (key == "epsg"){
		istringstream ss(value);
		int intvalue;
		if (!(ss >> intvalue)) return false;
	}	

	if (SMETCommon::all_decimal_header_values.find(key) != SMETCommon::all_decimal_header_values.end()){
		//check if value is a double value
		if (!SMETCommon::is_decimal(value)) return false;
	}

	if ((location_epsg == 15) || (location_wgs84 == 7))
		location_in_header = true;

	//Rudimentary check on keys: fields, units_offset, units_multiplier
	if ((key == "fields") || (key == "units_offset") || (key == "units_multiplier"))
		key_ok = check_fields(key, value);

	return key_ok;
}

bool SMETWriter::check_fields(const std::string& key, const std::string& value)
{
	vector<string> tmp_vec;
	
	size_t counter = SMETCommon::readLineToVec(value, tmp_vec);

	//Firstly: Check if number of fields is consistent
	if (nr_of_fields != 0){
		if (counter != nr_of_fields) return false;
	} else {
		nr_of_fields = counter;
	}
	
	size_t count_wgs84 = 0, count_epsg = 0;
	if (key == "fields"){
		//set<string> fieldnames; //this will help us locate duplicate fields

		//check if location is in data and if timestamp is present
		for (size_t ii = 0; ii<tmp_vec.size(); ii++){
			if (tmp_vec[ii] == "latitude") count_wgs84++;
			else if (tmp_vec[ii] == "longitude") count_wgs84++;
			else if (tmp_vec[ii] == "easting") count_epsg++;
			else if (tmp_vec[ii] == "northing") count_epsg++;

			if (tmp_vec[ii] == "altitude") {
				count_wgs84++;
				count_epsg++;
			}

			if (tmp_vec[ii] == "timestamp"){
				if (timestamp_present) return false; //no duplicate timestamp field allowed
				timestamp_present = true;
				timestamp_field = ii;
			}

			if (tmp_vec[ii] == "julian"){
				if (julian_present) return false; //no duplicate julian field allowed
				julian_present = true;
				julian_field = ii;
			}
		}

		if (count_wgs84 == 3)
			location_in_data_wgs84 = true;

		if (count_epsg == 3)
			location_in_data_epsg = true;
	} else {
		//every value in units_offset and units_multiplier must be a decimal
		vector<string>::iterator it;
		for (it = tmp_vec.begin(); it != tmp_vec.end(); it++)
			if (!SMETCommon::is_decimal(*it)) return false;
	}

	return true;
}

void SMETWriter::write_signature()
{
	fout << "SMET " << SMETCommon::smet_version << " ";
	if (smet_type == ASCII)
		fout << "ASCII" << endl;
	else
		fout << "BINARY" << endl;
}

bool SMETWriter::valid_header()
{
	if (mandatory_header_keys.size() != 3)
		return false;

	if (location_in_data_epsg){ //EPSG code must be in the header anyway
		map<string,string>::iterator it = header.find("epsg");
		if (it == header.end()){
			return false;
		} else {
			return true;
		}
	}

	if (location_in_header || location_in_data_wgs84)
		return true;

	return false;
}

void SMETWriter::write(const std::vector<std::string>& vec_timestamp, const std::vector<double>& data)
{
	fout.open(filename.c_str());
	if (fout.fail())
		throw SMETException("Could not open file '" + filename + "' for writing", SMET_AT);

	write_header(); //Write the header info, always in ASCII format

	if (nr_of_fields == 0){
		cleanup();
		return;
	}

	if (!timestamp_present) 
		throw SMETException("No timestamp present, use write(const vector<double>& data)", SMET_AT);

	const size_t nr_of_data_fields = nr_of_fields - 1;

	size_t nr_of_lines = data.size() / (nr_of_fields-1);
	if ((nr_of_lines != vec_timestamp.size()) || ((data.size() % (nr_of_fields-1)) != 0))
		throw SMETException("Inconsistency between vec_timestamp and data detected, recheck your data", SMET_AT);

	std::vector<double> current_data;
	current_data.resize(nr_of_fields-1); //HACK: check if nr_of_fields>0
	check_formatting();

	if (smet_type == ASCII){
		for (size_t ii=0; ii<nr_of_lines; ii++){
			const size_t offset = ii*(nr_of_fields-1);
			if (data.size() != 0)
				copy(data.begin()+offset, data.begin()+offset+nr_of_data_fields, current_data.begin());
			write_data_line_ascii(vec_timestamp[ii], current_data);	
		}
	} else {
		throw SMETException("Cannot write a binary file with a timestamp, use julian instead", SMET_AT);
	}

	cleanup();
}

void SMETWriter::write(const std::vector<double>& data)
{
	fout.open(filename.c_str());
	if (fout.fail())
		throw SMETException("Could not open file '" + filename + "' for writing", SMET_AT);

	write_header(); //Write the header info, always in ASCII format

	if (nr_of_fields == 0){
		cleanup();
		return;
	}

	size_t nr_of_lines = data.size() / nr_of_fields;
	if ((data.size() % nr_of_fields) != 0)
		throw SMETException("Inconsistency between data and header fields detected, recheck your data", SMET_AT);

	std::vector<double> current_data;
	current_data.resize(nr_of_fields); 
	check_formatting();

	if (smet_type == ASCII){
		for (size_t ii=0; ii<nr_of_lines; ii++){
			if (data.size() != 0)
				copy(data.begin()+ii*nr_of_fields, data.begin()+ii*nr_of_fields+nr_of_fields, current_data.begin());
			write_data_line_ascii("0000-01-01T00:00", current_data); //dummy time
		}
	} else {
		for (size_t ii=0; ii<nr_of_lines; ii++){
			if (data.size() != 0)
				copy(data.begin()+ii*nr_of_fields, data.begin()+ii*nr_of_fields+nr_of_fields, current_data.begin());
			write_data_line_binary(current_data);
		}
		
		file_is_binary = false;
	}

	cleanup();
}

void SMETWriter::write_header()
{
	map<string,string>::iterator it;

	if (!valid_header())
		throw SMETException("The header data you supplied is not valid, file cannot be written", SMET_AT);

	write_signature();

	fout << "[HEADER]" << endl;
	fout << "station_id       = " << header["station_id"] << endl;
	
	it = header.find("station_name");
	if (it != header.end()) 
		fout << "station_name     = " << it->second << endl;

	if (location_in_header){
		if (location_wgs84 == 7){
			fout << "latitude         = " << header["latitude"]  << endl;
			fout << "longitude        = " << header["longitude"] << endl;
			fout << "altitude         = " << header["altitude"]  << endl;
		}

		if (location_epsg == 15){
			fout << "easting          = " << header["easting"]   << endl;
			fout << "northing         = " << header["northing"]  << endl;
			if (location_wgs84 != 7)
				fout << "altitude         = " << header["altitude"]  << endl;
			fout << "epsg             = " << header["epsg"]  << endl;
		}
	} else {
		if (location_in_data_epsg)
			fout << "epsg             = " << header["epsg"]      << endl;
	}

	fout << "nodata           = " << header["nodata"] << endl;

	//Optional header keys:
	it = header.find("tz");
	if (it != header.end()) 
		fout << "tz               = " << it->second << endl;

	it = header.find("creation");
	if (it != header.end()) 
		fout << "creation         = " << it->second << endl;

	it = header.find("source");
	if (it != header.end()) 
		fout << "source           = " << it->second << endl;

	it = header.find("units_offset");
	if (it != header.end()) 
		fout << "units_offset     = " << it->second << endl;

	it = header.find("units_multiplier");
	if (it != header.end()) 
		fout << "units_multiplier = " << it->second << endl;

	it = header.find("comment");
	if (it != header.end()) 
		fout << "comment          = " << it->second << endl;

	fout << "fields           = " << header["fields"] << endl;
	fout << "[DATA]" << endl;
}

void SMETWriter::write_data_line_binary(const std::vector<double>& data)
{
	const char eoln = '\n';

	if (!file_is_binary){ //open it as binary file
		fout.close();
		fout.open(filename.c_str(), ios::out | ios::app | ios::binary); //reopen as binary file
		if (fout.fail()) 
			throw SMETException("Could not access file " + filename, SMET_AT);

		file_is_binary = true;
	}

	float val = 0;
	for (size_t ii = 0; ii < data.size(); ii++){
		if (julian_present && (julian_field == ii)){
			double julian = data[ii];
			fout.write((char*)&julian, sizeof(double)); //the julian date is written in 64bit IEEE754 precision
		} else {
			val = (float)data[ii];
			fout.write((char*)&val, sizeof(float)); //normal data fields are written in 32bit IEEE754 precision
		}
	}

	fout.write((char*)&eoln, sizeof(char));
}

void SMETWriter::write_data_line_ascii(const std::string& timestamp, const std::vector<double>& data)
{
	fout.fill(' ');
	fout << right;
	fout << fixed;

	if ((data.size() == 0) && timestamp_present) fout << timestamp;

	for (size_t ii = 0; ii < data.size(); ii++){
		if (ii > 0) fout << " ";
		if (timestamp_present && (timestamp_field == ii))	fout << timestamp << " ";

		fout << setw(ascii_width[ii]) << setprecision(ascii_precision[ii]);
		if (data[ii] == nodata_value) fout << nodata_string; //to have a nicer representation
		else fout << data[ii];
	}
	fout << endl;
}

void SMETWriter::check_formatting()
{
	size_t nr_of_data_fields = nr_of_fields;
	if (timestamp_present) nr_of_data_fields--;

	if ((ascii_precision.size() != nr_of_data_fields) || (ascii_width.size() != nr_of_data_fields)){
		ascii_precision.resize(nr_of_data_fields, 3);
		ascii_width.resize(nr_of_data_fields, 8);
	}

	if (julian_present){
		ascii_precision[julian_field] = 8;
		ascii_width[julian_field] = 16;
	}
}

void SMETWriter::set_width(const std::vector<size_t>& vec_width)
{
	ascii_width = vec_width;
}

void SMETWriter::set_precision(const std::vector<size_t>& vec_precision)
{
	ascii_precision = vec_precision;
}


SMETReader::SMETReader(const std::string& in_fname) : filename(in_fname), nr_of_fields(0), timestamp_present(false), 
										    julian_present(false), timestamp_field(0), julian_field(0),
										    isAscii(true), location_wgs84(0), location_epsg(0),
										    location_data_wgs84(0), location_data_epsg(0), nodata_value(-999.),
										    timestamp_interval(false), julian_interval(false),
										    timestamp_start("-4714-11-24T00:00"),
										    timestamp_end("9999-12-31T00:00")
{
	std::ifstream fin; //Input file streams
	fin.clear();
	fin.open (filename.c_str(), ios::in);
	if (fin.fail())
		throw SMETException("Could not open file '" + filename + "' for reading", SMET_AT);

	try {
		eoln = SMETCommon::getEoln(fin); //get the end of line character for the file

		read_header(fin);
		process_header();
	} catch(...){
		cleanup(fin); //closes file
		throw;
	}
	cleanup(fin); //closes file
}

SMETReader::~SMETReader()
{

}

void SMETReader::cleanup(std::ifstream& fin) throw()
{
	if (fin.is_open()) //close fin if open
		fin.close();
}

void SMETReader::convert_to_MKSA(const bool& in_mksa)
{
	mksa = in_mksa;
}

std::string SMETReader::get_field_name(const size_t& nr_of_field)
{
	if (nr_of_field < nr_of_fields){
		return vec_fieldnames[nr_of_fields];
	} else {
		throw SMETException("The field you're trying to access is out of bounds", SMET_AT);
	}
}

size_t SMETReader::get_nr_of_fields() const
{
	return nr_of_fields;
}

void SMETReader::get_units_conversion(std::vector<double>& offset, std::vector<double>& multiplier) const
{
	multiplier = vec_multiplier;
	offset = vec_offset;
}

bool SMETReader::location_in_header(const LocationType& type) const
{
	if ((location_epsg == 15) && (type == EPSG)){
		return true;
	} else if ((location_wgs84 == 7) && (type == WGS84)){
		return true;
	}

	return false;
}

bool SMETReader::location_in_data(const LocationType& type) const
{
	if ((location_epsg == 8) && (location_data_epsg == 7) && (type == EPSG)){
		return true;
	} else if ((location_data_wgs84 == 7) && (type == WGS84)){
		return true;
	}

	return false;
}

void SMETReader::process_header()
{
	vector<string> tmp_vec;
	map<string,string>::iterator it;
	set<string> obligatory_keys;
	for (it = header.begin(); it != header.end(); it++){
		if (SMETCommon::all_mandatory_header_keys.find(it->first) != SMETCommon::all_mandatory_header_keys.end())
			obligatory_keys.insert(it->first);
		
		if (it->first == "fields"){
			SMETCommon::readLineToVec(it->second, tmp_vec);
			string newfields = "";
			if (tmp_vec.size() > 0){
				for (size_t ii=0; ii<tmp_vec.size(); ii++){
					if (tmp_vec[ii] == "timestamp"){
						timestamp_present = true;
						timestamp_field = ii;
					} else {
						if (nr_of_fields > 0) newfields += " ";
						vec_fieldnames.push_back(tmp_vec[ii]);
						newfields += tmp_vec[ii];
						nr_of_fields++;
					}

					if (tmp_vec[ii] == "julian"){
						julian_present = true;
						julian_field = ii;
					}
					
					//Using WGS_84 coordinate system
					if (tmp_vec[ii] == "latitude")  location_data_wgs84 |= 1;
					else if (tmp_vec[ii] == "longitude") location_data_wgs84 |= 2;
					else if (tmp_vec[ii] == "altitude")  location_data_wgs84 |= 4;

					//Using an EPSG coordinate system
					if (tmp_vec[ii] == "easting")  location_data_epsg |= 1;
					else if (tmp_vec[ii] == "northing") location_data_epsg |= 2;
					else if (tmp_vec[ii] == "altitude") location_data_epsg |= 4;
				}
				it->second = newfields;
			}
		}

		if (it->first == "nodata")
			nodata_value = SMETCommon::convert_to_double(it->second);

		//Using WGS_84 coordinate system
		if (it->first == "latitude")  location_wgs84 |= 1;
		else if (it->first == "longitude") location_wgs84 |= 2;
		else if (it->first == "altitude")  location_wgs84 |= 4;

		//Using an EPSG coordinate system
		if (it->first == "easting")  location_epsg |= 1;
		else if (it->first == "northing") location_epsg |= 2;
		else if (it->first == "altitude") location_epsg |= 4;
		else if (it->first == "epsg")     location_epsg |= 8;

		//Now do some value checks
		if (it->first == "epsg"){
			istringstream ss(it->second);
			int intvalue;
			if (!(ss >> intvalue)) 
				throw SMETException("In " + filename + ": EPSG code not an integer number", SMET_AT);
		}	
	}

	if (get_header_value("units_offset") != ""){
		vector<string> tmp_offset;
		string offsetstring = "";
		SMETCommon::readLineToVec(get_header_value("units_offset"), tmp_offset);

		for (size_t ii=0; ii<tmp_offset.size(); ii++){
			if (!timestamp_present || (ii != timestamp_field)){
				if (ii>0) offsetstring += " ";
				offsetstring += tmp_offset[ii];
				vec_offset.push_back(SMETCommon::convert_to_double(tmp_offset[ii]));
			}
		}
		header["units_offset"] = offsetstring;
	} else {
		vec_offset.resize(nr_of_fields, 0.0);
	}

	if (get_header_value("units_multiplier") != ""){
		vector<string> tmp_multiplier;
		string multiplierstring = "";
		SMETCommon::readLineToVec(get_header_value("units_multiplier"), tmp_multiplier);

		for (size_t ii=0; ii<tmp_multiplier.size(); ii++){
			if (!timestamp_present || (ii != timestamp_field)){
				if (ii>0) multiplierstring += " ";
				multiplierstring += tmp_multiplier[ii];
				vec_multiplier.push_back(SMETCommon::convert_to_double(tmp_multiplier[ii]));
			}
		}
		header["units_multiplier"] = multiplierstring;
	} else {
		vec_multiplier.resize(nr_of_fields, 1.0);
	}

	if (obligatory_keys.size() != SMETCommon::all_mandatory_header_keys.size())
		throw SMETException("Not a valid SMET file, mandatory information in header is missing", SMET_AT);

	if ((location_wgs84 == 7) || (location_epsg == 15) || (location_data_wgs84 == 7) 
	    || ((location_epsg == 8) && (location_data_epsg == 7))){
		//location info present somewhere
	} else {
		throw SMETException("Not a valid SMET file, mandatory location info is missing (header and data)", SMET_AT);
	}
}

void SMETReader::read_header(std::ifstream& fin)
{
	//1. Read signature
	std::string line = "";
	vector<string> tmpvec;
	getline(fin, line, eoln); //read complete signature line
	SMETCommon::stripComments(line);
	SMETCommon::readLineToVec(line, tmpvec);
	checkSignature(tmpvec, isAscii);

	//2. Read Header
	while (!fin.eof() && (fin.peek() != '[')) //skip lines until '[' is found
		getline(fin, line, eoln);

	getline(fin, line, eoln);
	SMETCommon::stripComments(line);
	SMETCommon::trim(line);
	SMETCommon::toUpper(line);

	if (line != "[HEADER]")
		throw;// InvalidFormatException("Section " + line + " in "+ filename + " invalid", AT);

	while (!fin.eof() && (fin.peek() != '[')){ //Read until next section
		getline(fin, line, eoln);

		SMETCommon::stripComments(line);
		SMETCommon::trim(line);

		if (line != "") {
			if (!SMETCommon::readKeyValuePair(line, "=", header))
				throw SMETException("Invalid key value pair in section [Header]: " + line, SMET_AT);
		}

	}

	//Read [DATA] section tag
	getline(fin, line, eoln);
	SMETCommon::stripComments(line);
	SMETCommon::trim(line);
	SMETCommon::toUpper(line);

	if (line != "[DATA]")
		throw SMETException("Section " + line + " in "+ filename + " invalid, expected [DATA]", SMET_AT);

	data_start_fpointer = fin.tellg();
}

void SMETReader::checkSignature(const std::vector<std::string>& vecSignature, bool& isAscii)
{
	if ((vecSignature.size() != 3) || (vecSignature[0] != "SMET"))
		throw SMETException("The signature of file " + filename + " is invalid", SMET_AT);

	std::string version = vecSignature[1];
	if ((version != "0.9") && (version != "0.95") && (version != "0.99") 
	    && (version != "1.0") && (version != SMETCommon::smet_version))
		throw SMETException("Unsupported file format version for file " + filename, SMET_AT);

	if(version=="0.9" || version=="0.95" || version=="0.99" || version=="1.0") {
		cout << "[W] SMET specification 1.1 changes the priorities of units_multiplier and units_offset. "
			<< "Please check/update your files and bring them to 1.1!!" << endl;
	}

	const std::string type = vecSignature[2];
	if (type == "ASCII")
		isAscii = true;
	else if (type == "BINARY")
		isAscii = false;
	else
		throw SMETException("The 3rd column of file " + filename + " must be either ASCII or BINARY", SMET_AT);
}

void SMETReader::read(const std::string& i_timestamp_start, const std::string& i_timestamp_end, 
                      std::vector<std::string>& vec_timestamp, std::vector<double>& vec_data)
{
	if (!timestamp_present){
		read(vec_timestamp, vec_data);
	} else {
		timestamp_interval = true;
		timestamp_start = i_timestamp_start;
		timestamp_end = i_timestamp_end;

		read(vec_timestamp, vec_data);

		timestamp_interval = false;
	}
}

void SMETReader::read(const double& i_julian_start, const double& i_julian_end, std::vector<double>& vec_data)
{
	if (!julian_present){
		read(vec_data);
	} else {
		julian_interval = true;
		julian_start = i_julian_start;
		julian_end = i_julian_end;

		read(vec_data);

		julian_interval = false;
	}
}

void SMETReader::read(std::vector<std::string>& vec_timestamp, std::vector<double>& vec_data)
{
	if (!timestamp_present) 
		throw SMETException("Requesting to read timestamp when there is none present", SMET_AT);

	ifstream fin;
	fin.clear();
	fin.open (filename.c_str(), ios::in);
	if (fin.fail())
		throw SMETException("Could not open file '" + filename + "' for reading", SMET_AT);

	try {
		fin.seekg(data_start_fpointer); //jump to data start position in the file
		if (timestamp_interval && timestamp_present){
			map<string, streampos>::const_iterator it = map_timestamp_streampos.find(timestamp_start);
			if ( it != map_timestamp_streampos.end())
				fin.seekg(it->second);
		} else if (julian_interval && julian_present){
			map<double, streampos>::const_iterator it = map_julian_streampos.find(julian_start);
			if ( it != map_julian_streampos.end())
				fin.seekg(it->second);
		}

		if (fin.fail() || fin.bad())
			fin.seekg(data_start_fpointer);

		if (isAscii)
			read_data_ascii(fin, vec_timestamp, vec_data);
		else
			throw SMETException("Binary SMET file has no field timestamp, only julian date", SMET_AT);
	} catch(...) {
		cleanup(fin);
		throw; 
	}

	cleanup(fin);
}

void SMETReader::read(std::vector<double>& vec_data)
{
	if (timestamp_present) 
		SMETException("Requesting not to read timestamp when there is one present", SMET_AT);
	
	vector<string> tmp_vec;

	ios_base::openmode mode = ios::in;
	if (!isAscii)
		mode = ios::in | ios::binary;

	ifstream fin;
	fin.open (filename.c_str(), ios::in);
	if (fin.fail()) 
		throw SMETException("Could not open file '" + filename + "' for reading", SMET_AT);
	
	try {
		fin.seekg(data_start_fpointer); //jump to data start position in the file

		if (isAscii)
			read_data_ascii(fin, tmp_vec, vec_data);
		else
			read_data_binary(fin, vec_data);
	} catch(...) {
		cleanup(fin);
		throw; 
	}

	cleanup(fin);
}

void SMETReader::read_data_ascii(std::ifstream& fin, std::vector<std::string>& vec_timestamp, std::vector<double>& vec_data)
{
	size_t nr_of_data_fields = nr_of_fields;
	if (timestamp_present) nr_of_data_fields++;

	vector<string> tmp_vec;
	string line;
	streampos current_fpointer = -1;

	while (!fin.eof()){
		line = "";
		streampos tmp_fpointer = fin.tellg();
		getline(fin, line, eoln);
		SMETCommon::stripComments(line);
		SMETCommon::trim(line);
		if (line == "") continue; //Pure comment lines and empty lines are ignored

		if (SMETCommon::readLineToVec(line, tmp_vec) == nr_of_data_fields){
			size_t shift = 0;
			if (julian_interval && julian_present){
				double current_julian = SMETCommon::convert_to_double(tmp_vec[julian_field]);
				if (current_julian < julian_start) 
					continue; //skip lines that don't hold the dates we're interested in
				else if (current_julian > julian_end)
					break; //skip the rest of the file
			}

			if (timestamp_interval && timestamp_present){
				const string& current_timestamp = tmp_vec[timestamp_field];
				if (current_timestamp < timestamp_start)
					continue; //skip lines that don't hold the dates we're interested in
				else if (current_timestamp > timestamp_end)
					break; //skip the rest of the file
			}

			for (size_t ii=0; ii<tmp_vec.size(); ii++){
				if (ii == timestamp_field){
					vec_timestamp.push_back(tmp_vec[ii]);
					shift = 1;
				} else {
					double tmp = SMETCommon::convert_to_double(tmp_vec[ii]);
					if ((mksa) && (tmp != nodata_value)){
						tmp *= vec_multiplier[ii-shift];
						tmp += vec_offset[ii-shift];
					}
					vec_data.push_back(tmp);
				}
			}
			current_fpointer = tmp_fpointer;
		}
	}

	if (current_fpointer != -1){
		if (timestamp_interval && timestamp_present)
			map_timestamp_streampos[timestamp_end] = current_fpointer;
		else if (julian_interval && julian_present)
			map_julian_streampos[julian_end] = current_fpointer;
	}
}

void SMETReader::read_data_binary(std::ifstream& fin, std::vector<double>& vec_data)
{
	while (!fin.eof()){
		double julian = -1.0;
		for (size_t ii=0; ii<nr_of_fields; ii++){
			if (julian_present && (ii == julian_field)){
				fin.read(reinterpret_cast < char * > (&julian), sizeof(double));
				vec_data.push_back(julian);
			} else {
				float tmpval;
				fin.read(reinterpret_cast < char * > (&tmpval), sizeof(float));
				double val = (double)tmpval;

				vec_data.push_back(val);
			}

			if (mksa){
				double& tmp = vec_data.back();
				if (tmp != nodata_value){
					tmp *= vec_multiplier[ii];
					tmp += vec_offset[ii];
				}
			}
		}
		
		if (julian_present && julian_interval){
			if ((julian < julian_start) || (julian > julian_end)){
				for (size_t ii=0; ii<nr_of_fields; ii++){
					vec_data.pop_back();
				}
			}
		}

		char c;
		fin.read(&c, sizeof(char));
		if (c != '\n')
			throw SMETException("Corrupted data in section [DATA]", SMET_AT);
	}
}

double SMETReader::get_header_doublevalue(const std::string& key) const
{
	map<string,string>::const_iterator it = header.find(key);
	if (it != header.end())
		return SMETCommon::convert_to_double(it->second);
	
	return nodata_value;
}

std::string SMETReader::get_header_value(const std::string& key) const
{
	map<string,string>::const_iterator it = header.find(key);
	if (it != header.end())
		return it->second;
	
	return "";
}

bool SMETReader::contains_timestamp() const
{
	return timestamp_present;
}

} //end namespace

