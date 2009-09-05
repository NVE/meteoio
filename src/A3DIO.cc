#include "A3DIO.h"

using namespace std;

//A3DIO::A3DIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

//Main constructor
A3DIO::A3DIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	//Nothing else so far
}

//Copy constructor
A3DIO::A3DIO(const A3DIO& aio) : IOInterface(NULL), cfg(aio.cfg)
{
	//Nothing else so far
}

A3DIO::A3DIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	//Nothing else so far
}

A3DIO::~A3DIO() throw()
{
	cleanup();
}

//Clone function
//A3DIO* A3DIO::clone() const { return new A3DIO(*this); }

void A3DIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
}

void A3DIO::read2DGrid(Grid2DObject&, const string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);	
}

void A3DIO::readDEM(DEMObject& dem_out)
{
	//Nothing so far
	(void)dem_out;
	throw IOException("Nothing implemented here", AT);	
}

void A3DIO::readLanduse(Grid2DObject& landuse_out)
{
	//Nothing so far
	(void)landuse_out;
	throw IOException("Nothing implemented here", AT);	
}

void A3DIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	//Nothing so far
	(void)date_in;
	(void)da_out;
	throw IOException("Nothing implemented here", AT);	
}

void A3DIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, vector< vector<MeteoData> >& vecMeteo)
{
	vector< vector<StationData> > vecStation;
	readMeteoData(dateStart, dateEnd, vecMeteo, vecStation);
}


void A3DIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
					 std::vector< std::vector<MeteoData> >& vecMeteo, 
					 std::vector< std::vector<StationData> >& vecStation,
					 const unsigned int& stationindex)
{
	//if dateStart and dateEnd are the same: return exact match for date
	//if dateStart > dateEnd: return first data set with date > dateStart
	//read in all data starting with dateStart until dateEnd
	//if there is no data at all then the vector will be empty, no exception will be thrown
	(void)stationindex;
	vecMeteo.clear();
	vecStation.clear();

	//first from meteo1d.txt
	read1DMeteo(dateStart, dateEnd, vecMeteo, vecStation);

	//then all corresponding data sets from the 2d meteo files
	//note: they have to be entirely corresponding (for every date)
	try {
		read2DMeteo(vecMeteo, vecStation);
	} catch(exception& e){
		cerr << "[I] No meteo2d data found or error while reading it, using only Meteo1D data: " 
			<< endl << "\t" << e.what() << endl;
	}
}

void A3DIO::convertUnits(MeteoData& meteo)
{
	//converts C to Kelvin, converts RH to [0,1]
	if(meteo.ta!=nodata) {
		meteo.ta=C_TO_K(meteo.ta);
	}
	if(meteo.tsg!=nodata) {
		meteo.tsg=C_TO_K(meteo.tsg);
	}
	if(meteo.rh!=nodata) {
		meteo.rh /= 100.;
	}
}

void A3DIO::read1DMeteo(const Date_IO& dateStart, const Date_IO& dateEnd, 
				    vector< vector<MeteoData> >& vecMeteo, vector< vector<StationData> >& vecStation)
{
	double latitude=IOUtils::nodata, longitude=IOUtils::nodata, 
		xcoord=IOUtils::nodata, ycoord=IOUtils::nodata, altitude=IOUtils::nodata;
	string tmp="", line="";
	Date_IO tmp_date;
	vector<string> tmpvec;
	map<string, string> header; // A map to save key value pairs of the file header
	MeteoData tmpdata;
	StationData sd;
	bool eofreached = false;

	cfg.getValue("METEOPATH", tmp); 
	tmp += "/meteo1d.txt";

	if (!IOUtils::fileExists(tmp)) {
		throw FileNotFoundException(tmp, AT);
	}

	fin.clear();
	fin.open (tmp.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(tmp,AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
   
	//Go through file, save key value pairs
	try {
		//Read in station meta data
		IOUtils::readKeyValueHeader(header, fin, 5, "="); //Read in 5 lines as header
		IOUtils::getValueForKey(header, "Latitude", latitude);
		IOUtils::getValueForKey(header, "Longitude", longitude);
		IOUtils::getValueForKey(header, "X_Coord", xcoord);
		IOUtils::getValueForKey(header, "Y_Coord", ycoord);
		IOUtils::getValueForKey(header, "Altitude", altitude);

		string coordsys="", coordparam="";
		try {
			cfg.getValue("COORDIN", coordsys);
			cfg.getValue("COORDPARAM", coordparam); 
		} catch(std::exception& e) {
			//problems while reading values for COORDIN or COORDPARAM
			cerr << "[E] reading configuration file: " << "\t" << e.what() << endl;
			throw;
		}
		MapProj mymapproj(coordsys, coordparam);

		//calculate coordinates if necessary
		if(latitude==IOUtils::nodata || longitude==IOUtils::nodata) {
			if(xcoord==IOUtils::nodata || ycoord==IOUtils::nodata) {
				throw InvalidFormatException("Too many nodata values for coordinates conversion in file " + tmp, AT);
			}
			mymapproj.convert_to_WGS84(xcoord, ycoord, latitude, longitude);
		} else if(xcoord!=IOUtils::nodata || ycoord!=IOUtils::nodata) {
			double tmp_lat, tmp_lon;
			mymapproj.convert_to_WGS84(xcoord, ycoord, tmp_lat, tmp_lon);
			if(!checkEpsilonEquality(latitude, tmp_lat, 1.e-4) || !checkEpsilonEquality(longitude, tmp_lon, 1.e-4)) {
				throw InvalidArgumentException("Latitude/longitude and Xcoord/Ycoord don't match in header of file " + tmp, AT);
			}
		}
		sd.setStationData(xcoord, ycoord, altitude, "", latitude, longitude);

		//Read one line, construct Date_IO object and see whether date is greater or equal than the date_in object
		IOUtils::skipLines(fin, 1, eoln); //skip rest of line
		
		//Loop going through the data sequentially until dateStart is found
		do {
			getline(fin, line, eoln); //read complete line
			eofreached = readMeteoDataLine(line, tmpdata, tmp);
			tmpdata.cleanData();
			convertUnits(tmpdata);

		} while((tmpdata.date < dateStart) && (!eofreached));

		if ((dateEnd < dateStart) && (!eofreached)){ //Special case
			vecMeteo.push_back( vector<MeteoData>() );
			vecStation.push_back( vector<StationData>() );
			vecMeteo[0].push_back(tmpdata);
			vecStation[0].push_back(sd);
		} else if ((tmpdata.date <= dateEnd)  && (!eofreached)) {
			vecMeteo.push_back( vector<MeteoData>() );
			vecStation.push_back( vector<StationData>() );
		}

		while ((tmpdata.date <= dateEnd)  && (!eofreached)) {
			//At this point tmpdata.date is >= dateStart
			vecMeteo[0].push_back(tmpdata);
			vecStation[0].push_back(sd);

			getline(fin, line, eoln); //read complete line
			eofreached = readMeteoDataLine(line, tmpdata, tmp);
			tmpdata.cleanData();
			convertUnits(tmpdata);
		}
		//cout << "Size of buffer: " << vecMeteo[0].size() << "   " << tmp_date.toString() << endl;
	} catch(...) {
		cout << "[e] " << AT << ": "<< endl;
		cleanup();
		throw;
	}

	cleanup();
}

bool A3DIO::readMeteoDataLine(std::string& line, MeteoData& tmpdata, string filename)
{
	Date_IO tmp_date;
	int tmp_ymdh[4];
	vector<string> tmpvec;
	double tmp_values[6];

	if (IOUtils::readLineToVec(line, tmpvec) != 10) {
		return true;
		//throw InvalidFormatException("Premature End of Line or no data for date " + date_in.toString() + " found in File " + filename, AT);
	}
      
	for (int ii=0; ii<4; ii++) {
		if (!IOUtils::convertString(tmp_ymdh[ii], tmpvec.at(ii), std::dec))
			throw InvalidFormatException(filename + ": " + line, AT);
	}

	tmp_date.setDate(tmp_ymdh[0],tmp_ymdh[1],tmp_ymdh[2],tmp_ymdh[3]);

	//Read rest of line with values ta, iswr, vw, rh, ea, nswc
  
	for (int ii=0; ii<6; ii++) { //go through the columns
		if (!IOUtils::convertString(tmp_values[ii], tmpvec.at(ii+4), std::dec)) {
			throw InvalidFormatException(filename + ": " + line, AT);
		}
	}
      
	tmpdata.setMeteoData(tmp_date, tmp_values[0], tmp_values[1], tmp_values[2], nodata, tmp_values[3], tmp_values[4], tmp_values[5], nodata, nodata, nodata);
	
	return false;
}


/*
  Preamble: Files are in METEOFILE directory. 4 types of files: 
  prec????.txt == nswc
  rh????.txt == rh
  ta????.txt == ta
  wspd????.txt == vw
  
  Remarks: The headers of the files may defer - for each unique 
  StationData one MeteoData and one StationData object will be created
*/
 void A3DIO::read2DMeteo(vector< vector<MeteoData> >& vecMeteo, vector< vector<StationData> >& vecStation)
{
	
	unsigned int stations=0, bufferindex=0;
	map<string, unsigned int> hashStations = map<string, unsigned int>();
	vector<string> filenames = vector<string>();

	//Requirement: meteo1D data must exist:
	if ((vecMeteo.size() == 0) || (vecMeteo[0].size() == 0))
		return;
	
	//1D and 2D data must correspond, that means that if there is 1D data
	//for a certain date (e.g. 1.1.2006) then 2D data must exist (prec2006.txt etc), 
	//otherwise throw FileNotFoundException
	constructMeteo2DFilenames(vecMeteo[0][0].date, filenames);

	stations = getNrOfStations(filenames, hashStations);
	cerr << "[I] Number of 2D meteo stations: " << stations << endl;

	if (stations < 1) {
		throw InvalidFormatException("No StationData found in 2D Meteo Files", AT); 
	} 

	vector<StationData> tmpvecS = vector<StationData>(stations); //stores unique stations

	try {
		read2DMeteoHeader(filenames[0], hashStations, tmpvecS);
		read2DMeteoHeader(filenames[1], hashStations, tmpvecS);
		read2DMeteoHeader(filenames[2], hashStations, tmpvecS);
		read2DMeteoHeader(filenames[3], hashStations, tmpvecS);
		if(IOUtils::fileExists(filenames[4])) { //for keeping dw optional
			read2DMeteoHeader(filenames[4], hashStations, tmpvecS);
		}

		//init vecStation with proper StationData, vecMeteo with nodata
		for (unsigned int jj=0; jj<tmpvecS.size(); jj++){
			vecMeteo.push_back( vector<MeteoData>() );
			vecStation.push_back( vector<StationData>() );
			for (unsigned int ii=0; ii<vecMeteo[0].size(); ii++){
				//NOTE: there needs to be the same amount of 1D and 2D data
				vecStation[jj+1].push_back(tmpvecS[jj]);
				vecMeteo[jj+1].push_back(MeteoData());
			}
		}

		do {
			unsigned int currentindex = bufferindex;
			read2DMeteoData(filenames[0], "nswc", hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[1], "rh", hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[2], "ta", hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[3], "vw", hashStations, vecMeteo, bufferindex);
			
			if(IOUtils::fileExists(filenames[4])) { //for keeping dw optional
				bufferindex = currentindex;
				read2DMeteoData(filenames[4], "dw", hashStations, vecMeteo, bufferindex);
			}
			//cerr << "bufferindex: " << bufferindex << "  Expected size()" << vecMeteo[0].size() << endl;

			if (bufferindex < (vecMeteo[0].size())) { //number of 1D meteo data
				//construct new filenames for the continued buffering
				constructMeteo2DFilenames(vecMeteo[0][bufferindex].date, filenames);
			}
		} while(bufferindex < (vecMeteo[0].size()));
	} catch(...) {
		//clear all 2D meteo data if error occurs
		if (vecMeteo.size() > 1)
			vecMeteo.erase(vecMeteo.begin()+1, vecMeteo.end());
		
		if (vecStation.size() > 1)
			vecStation.erase(vecStation.begin()+1, vecStation.end());

		cleanup();
		throw;
	}
  
	//clean data and convert the units
	for (unsigned int ii=1; ii<vecMeteo.size(); ii++) { //loop over all stations except 1D Meteo
		for (unsigned int jj=0; jj<vecMeteo[ii].size(); jj++){ //Meteo1D data already cleaned
			vecMeteo[ii][jj].cleanData();
			convertUnits(vecMeteo[ii][jj]);
		}
	}
}

void A3DIO::constructMeteo2DFilenames(const Date_IO& date_in, vector<string>& filenames)
{
	int year=0, dummy=0;
	string tmp;
	stringstream ss;

	filenames.clear();
 
	date_in.getDate(year, dummy, dummy);
	ss << year;

	cfg.getValue("METEOPATH", tmp); 

	string precFilename = tmp + "/prec" + ss.str() + ".txt";
	string rhFilename = tmp + "/rhum" + ss.str() + ".txt";
	string taFilename = tmp + "/tair" + ss.str() + ".txt";
	string wspdFilename = tmp + "/wspd" + ss.str() + ".txt";
	string wdirFilename = tmp + "/wdir" + ss.str() + ".txt";

	filenames.push_back(precFilename);
	filenames.push_back(rhFilename);
	filenames.push_back(taFilename);
	filenames.push_back(wspdFilename);
	filenames.push_back(wdirFilename);

	for (unsigned int ii=0; ii<filenames.size(); ii++) {
		if (!IOUtils::fileExists(filenames[ii]) && ii<4) { //for keeping dw optional
			throw FileNotFoundException(filenames[ii], AT);
		}
	}
}


unsigned int A3DIO::getNrOfStations(vector<string>& filenames, map<string, unsigned int>& hashStations)
{
	vector<string> tmpvec;
	string line_in="";
  
	for (unsigned int ii=0; ii<filenames.size(); ii++) {
		//cout << *it << endl;
		string filename = filenames[ii];

		fin.clear();
		fin.open (filename.c_str(), ifstream::in);
		if (fin.fail()) {
			if(ii==4) { //for keeping dw optional
				return (hashStations.size());
			}
			throw FileAccessException(filename, AT);
		}
  
		char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

		IOUtils::skipLines(fin, 4, eoln); 
		getline(fin, line_in, eoln); //5th line holds the names of the stations
		unsigned int cols = IOUtils::readLineToVec(line_in, tmpvec);
		if ( cols > 4) { // if there are any stations
			//check each station name and whether it's already hashed, otherwise: hash!
			for (unsigned int ii=4; ii<cols; ii++) {
				unsigned int tmp_int = hashStations.count(tmpvec.at(ii));
				if (tmp_int == 0) {
					hashStations[tmpvec.at(ii)] = hashStations.size();
				}
			} 
		}
		cleanup();      
	}  
  
	return (hashStations.size());
}

void A3DIO::read2DMeteoData(const string& filename, const string& parameter, 
					   map<string,unsigned int>& hashStations, 
					   vector< vector<MeteoData> >& vecM, unsigned int& bufferindex)
{
	
	string line_in = "";
	unsigned int columns;
	vector<string> tmpvec, vec_names;
	Date_IO tmp_date;
	int tmp_ymdh[4];

	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
    
	IOUtils::skipLines(fin, 4, eoln); //skip first 4 lines
	getline(fin, line_in, eoln); //line containing UNIQUE station names
	columns = IOUtils::readLineToVec(line_in, vec_names);
	if (columns < 4) {
		throw InvalidFormatException("Premature end of line in file " + filename, AT);
	}

	MeteoData& lastMeteoData = vecM[0][vecM[0].size()-1]; //last time stamp in buffer of 1D meteo

	do {
		getline(fin, line_in, eoln); 
		string tmpline = line_in;
		IOUtils::trim(tmpline);

		if (tmpline=="") {
			break;
		}

		if (IOUtils::readLineToVec(line_in, tmpvec)!=columns) { //Every station has to have its own column
			throw InvalidFormatException("Premature End of Line or no data for date " + vecM[0][bufferindex].date.toString() 
								    + " found in File " + filename, AT);
		}
    
		for (int ii=0; ii<4; ii++) {
			if (!IOUtils::convertString(tmp_ymdh[ii], tmpvec[ii], std::dec)) {
				throw InvalidFormatException("Check date columns in " + filename, AT);
			}
		}
		tmp_date.setDate(tmp_ymdh[0],tmp_ymdh[1],tmp_ymdh[2],tmp_ymdh[3]);

		MeteoData& currentMeteoData = vecM[0][bufferindex]; //1D Element to synchronize date
		if (tmp_date == currentMeteoData.date) {
			//Read in data
			for (unsigned int ii=4; ii<columns; ii++) {
				unsigned int stationnr = hashStations[vec_names.at(ii)]; 
				MeteoData& tmpmd = vecM[stationnr][bufferindex];
				tmpmd.date = tmp_date;

				if (parameter == "nswc") {
					if (!IOUtils::convertString(tmpmd.nswc, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For nswc value in " + filename + "  for date " + tmpmd.date.toString(), AT);
					}
	  
				} else if (parameter == "rh") {
					if (!IOUtils::convertString(tmpmd.rh, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For rh value in " + filename + "  for date " + tmpmd.date.toString(), AT);
					}
	  
				} else if (parameter == "ta") {
					if (!IOUtils::convertString(tmpmd.ta, tmpvec[ii], std::dec)) 
						throw ConversionFailedException("For ta value in " + filename + "  for date " + tmpmd.date.toString(), AT);
    
				} else if (parameter == "vw") {
					if (!IOUtils::convertString(tmpmd.vw, tmpvec[ii], std::dec)) { 
						throw ConversionFailedException("For vw value in " + filename + "  for date " + tmpmd.date.toString(), AT);
					}
				} else if (parameter == "dw") {
					if (!IOUtils::convertString(tmpmd.dw, tmpvec[ii], std::dec)) { 
						throw ConversionFailedException("For dw value in " + filename + "  for date " + tmpmd.date.toString(), AT);
					}
				}
			}

			bufferindex++;
		}
	} while((tmp_date<lastMeteoData.date) && (!fin.eof()));
	
	cleanup();
}

void A3DIO::read2DMeteoHeader(const string& filename, map<string,unsigned int>& hashStations, vector<StationData>& vecS)
{
	string line_in = "";
	unsigned int columns = 0;
	vector<string> vec_altitude, vec_xcoord, vec_ycoord, vec_names;

	//Build MapProj object to convert easting/northing values to lat/long in WGS84
	string coordsys="", coordparam="";
	try {
		cfg.getValue("COORDIN", coordsys);
		cfg.getValue("COORDPARAM", coordparam); 
	} catch(std::exception& e){
		//problems while reading values for COORDIN or COORDPARAM
	}
	MapProj mymapproj(coordsys, coordparam);

	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}
  
	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	IOUtils::skipLines(fin, 1, eoln);

	//Read all relevant lines in
	getline(fin, line_in, eoln); //Altitude
	columns = IOUtils::readLineToVec(line_in, vec_altitude);

	getline(fin, line_in, eoln); //xcoord
	if (IOUtils::readLineToVec(line_in, vec_xcoord) != columns) {
		throw InvalidFormatException("Column count doesn't match from line to line in " + filename, AT);
	}

	getline(fin, line_in, eoln); //ycoord
	if (IOUtils::readLineToVec(line_in, vec_ycoord) != columns) {
		throw InvalidFormatException("Column count doesn't match from line to line in " + filename, AT);
	}

	getline(fin, line_in, eoln); //names
	if (IOUtils::readLineToVec(line_in, vec_names) != columns) {
		throw InvalidFormatException("Column count doesn't match from line to line in " + filename, AT);
	}

	for (unsigned int ii=4; ii<columns; ii++) {
		unsigned int stationnr = hashStations[vec_names.at(ii)]; 
		if ((!IOUtils::convertString(vecS.at(stationnr-1).altitude, vec_altitude.at(ii), std::dec))
		    || (!IOUtils::convertString(vecS.at(stationnr-1).eastCoordinate, vec_xcoord.at(ii), std::dec))
		    || (!IOUtils::convertString(vecS.at(stationnr-1).northCoordinate, vec_ycoord.at(ii), std::dec))
		    || (!IOUtils::convertString(vecS.at(stationnr-1).stationName, vec_names.at(ii), std::dec))) {
			throw ConversionFailedException("Conversion of station description failed in " + filename, AT);  
		}

		//calculate coordinates if necessary
		if(vecS[stationnr-1].latitude==IOUtils::nodata || vecS[stationnr-1].longitude==IOUtils::nodata) {
			if(vecS[stationnr-1].eastCoordinate==IOUtils::nodata || vecS[stationnr-1].northCoordinate==IOUtils::nodata) {
				throw InvalidFormatException("Too many nodata values for coordinates conversion in file " + filename, AT);
			}
			mymapproj.convert_to_WGS84(vecS[stationnr-1].eastCoordinate, vecS[stationnr-1].northCoordinate, vecS[stationnr-1].latitude, vecS[stationnr-1].longitude);
		}
	}

	cleanup();
}

void A3DIO::readSpecialPoints(CSpecialPTSArray& pts)
{
	string filename="", line_in="";
	vector<string> tmpvec;
	vector< pair<int,int> > mypts;

	cfg.getValue("SPECIALPTSFILE", filename); // cout << tmp << endl;
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}

	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename,AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	while (!fin.eof()) {
		getline(fin, line_in, eoln); 

		if (IOUtils::readLineToVec(line_in, tmpvec)==2) { //Try to convert
			int x, y;
			if (!IOUtils::convertString(x, tmpvec.at(0), std::dec)) {
				throw ConversionFailedException("Conversion of a value failed in " + filename + " line: " + line_in, AT);        
			}

			if (!IOUtils::convertString(y, tmpvec.at(1), std::dec)) {
				throw ConversionFailedException("Conversion of a value failed in " + filename + " line: " + line_in, AT);        
			}

			pair<int,int> tmppair(x,y);
			mypts.push_back(tmppair);
		}
	}
	cleanup();

	//Now put everything into that legacy struct CSpecialPTSArray
	pts.resize(mypts.size());

	for (unsigned int jj=0; jj<mypts.size(); jj++) {
		pts[jj].ix = mypts.at(jj).first;
		pts[jj].iy = mypts.at(jj).second;
	}
}

void A3DIO::write2DGrid(const Grid2DObject&, const string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);	
}


/*extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}

	void* loadObject(const string& classname, const string& filename) {
		if(classname == "A3DIO") {
			cerr << "Creating handle to " << classname << endl;
			//return new A3DIO(deleteObject);
			return new A3DIO(deleteObject, filename);
		}
		cerr << "Could not load " << classname << endl;
		return NULL;
	}
}

*/
