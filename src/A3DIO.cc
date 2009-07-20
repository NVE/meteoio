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
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void A3DIO::get2DGridSize(int& nx, int& ny)
{
	string filename="";
	map<string, string> header; // A map to save key value pairs of the file header

	cfg.getValue("DEMFILE", filename); // cout << tmp << endl;

	if (!IOUtils::validFileName(filename)) {
		throw InvalidFileNameException(filename, AT);
	}
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}
  
	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);    
	}
   
	//Go through file, save key value pairs
	try {
		IOUtils::readKeyValueHeader(header, fin, 2, " ");
		IOUtils::getValueForKey(header, "ncols", nx);
		IOUtils::getValueForKey(header, "nrows", ny);
	} catch(std::exception& e) {
		cout << "[E] " << AT << ": "<< endl;
		nx=0;
		ny=0;
		cleanup();
		throw;
	}
	cleanup();
}

void A3DIO::read2DGrid(Grid2DObject& grid_out, const string& filename){

	int i_ncols, i_nrows;
	unsigned int ncols, nrows;
	double xllcorner, yllcorner, cellsize, nodata;
	double latitude, longitude;
	double tmp_val;
	vector<string> tmpvec;
	string line="";
	map<string, string> header; // A map to save key value pairs of the file header

	if (!IOUtils::validFileName(filename)) {
		throw InvalidFileNameException(filename, AT);
	}
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}
  
	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}
  
	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file
   
	//Go through file, save key value pairs
	try {
		IOUtils::readKeyValueHeader(header, fin, 6, " ");
		IOUtils::getValueForKey(header, "ncols", i_ncols);
		IOUtils::getValueForKey(header, "nrows", i_nrows);
		IOUtils::getValueForKey(header, "xllcorner", xllcorner);
		IOUtils::getValueForKey(header, "yllcorner", yllcorner);
		IOUtils::getValueForKey(header, "cellsize", cellsize);
		IOUtils::getValueForKey(header, "NODATA_value", nodata);

		if ((i_ncols==0) || (i_nrows==0)) {
			throw IOException("Number of rows or columns in 2D Grid given is zero, in file: " + filename, AT);
		}
		if((i_ncols<0) || (i_nrows<0)) {
			throw IOException("Number of rows or columns in 2D Grid read as \"nodata\", in file: " + filename, AT);
		}
		ncols = (unsigned int)i_ncols;
		nrows = (unsigned int)i_nrows;

		//compute WGS coordinates (considered as the true reference)
		//CH1903_to_WGS84(xllcorner, yllcorner, latitude, longitude);

		//HACK: check how we can input coordinates as WGS84 directly.
		//this HACK is a very cheap "extension" of the dem file format...
		if(xllcorner<360. || yllcorner<360.) {
			latitude = xllcorner;
			longitude = yllcorner;
		} else {
			CH1903_to_WGS84(xllcorner, yllcorner, latitude, longitude);
		}
    
		//Initialize the 2D grid
		grid_out.set(ncols, nrows, xllcorner, yllcorner, latitude, longitude, cellsize);
		
		//Read one line after the other and parse values into Grid2DObject
		for (unsigned int kk=nrows-1; (kk < nrows); kk--) {
			getline(fin, line, eoln); //read complete line
			//cout << "k:" << kk << "\n" << line << endl;

			if (IOUtils::readLineToVec(line, tmpvec) != ncols) {
				throw InvalidFormatException("Premature End " + filename, AT);
			}
			
			for (unsigned int ll=0; ll < ncols; ll++){
				if (!IOUtils::convertString(tmp_val, tmpvec.at(ll), std::dec)) {
					throw ConversionFailedException("For Grid2D value in line: " + line + " in file " + filename, AT);
				}
				
				if(tmp_val<=nodata) {
					//replace file's nodata by uniform, internal nodata
					grid_out.grid2D(ll, kk) = IOUtils::nodata;
				} else {
					grid_out.grid2D(ll, kk) = tmp_val;
				}
			}
		}
	} catch(std::exception& e) {
		cleanup();
		throw;
	}
	cleanup();
}

void A3DIO::readDEM(Grid2DObject& dem_out)
{
	string filename="";

	cfg.getValue("DEMFILE", filename); // cout << tmp << endl;

	read2DGrid(dem_out, filename);
}

void A3DIO::readLanduse(Grid2DObject& landuse_out)
{
	string filename="";

	cfg.getValue("LANDUSEFILE", filename); // cout << tmp << endl;

	read2DGrid(landuse_out, filename);
}

void A3DIO::readMeteoData(const Date_IO& date_in, vector<MeteoData>& meteo_out)
{
	vector<StationData> vecStation;
	readMeteoData(date_in, meteo_out, vecStation);
}

void A3DIO::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{
	//See whether data is already buffered
	//Filter
	//Resample if necessary
	//Filter resampled value
	//return that value
	/*
	  1. FilterFacade filterFacade("filters.ini");
	  2. filterFacade.getMinimalWindow(minNbPoints, minDeltaTime);
	  3. fiterFacade.doCheck(unfilteredMeteoBuffer, filteredMeteoBuffer,
	*/

	//FilterFacade filterFacade("filters.ini");

	unsigned int index = MeteoBuffer::npos;
	if (unfilteredMeteoBuffer.size()>0) {//check whether meteo data for the date exists in buffer
		index = unfilteredMeteoBuffer[0].seek(date_in);
	}
	//cerr << "Buffered MeteoBuffers: " << unfilteredMeteoBuffer.size() << "  Found data at index: " << index << endl;

	if (( index != MeteoBuffer::npos ) && (unfilteredMeteoBuffer[0].getMeteoData(index).date == date_in)) {//in the buffer, return from buffer
		cerr << "[I] Buffered data found for date: " << date_in.toString() << endl;

	} else if (index == MeteoBuffer::npos) { //not in buffer, rebuffering needed
		cerr << "[I] Data for date " << date_in.toString() << " not found in buffer, rebuffering" << endl;
		unfilteredMeteoBuffer.clear();
		filteredMeteoBuffer.clear();

		MeteoData md;
		read1DMeteo(date_in, md);
		index = unfilteredMeteoBuffer[0].seek(date_in);
		//cerr << "read 1D meteo" << endl;
    
		try {
			read2DMeteo(date_in, vecMeteo, vecStation);
		} catch(exception& e){
			cerr << "[I] No meteo2d data found or error while reading it, using only Meteo1D data: " << e.what() << endl;
      
			if (unfilteredMeteoBuffer.size()>1) {
				unfilteredMeteoBuffer.erase(unfilteredMeteoBuffer.begin()+1, unfilteredMeteoBuffer.end());
			}
		}
	}

	// RESAMPLING
	if (( index != MeteoBuffer::npos ) && (unfilteredMeteoBuffer[0].getMeteoData(index).date != date_in)) {//in the buffer, resampling needed 
		cerr << "[I] Resampling required for date: " << date_in.toString() << endl;
    
		Meteo1DResampler mresampler;
		for (unsigned int ii=0; ii<unfilteredMeteoBuffer.size(); ii++) {
			mresampler.resample(index, date_in, unfilteredMeteoBuffer[ii]);
		}
	}

	//FILTERING

	//fill the vectors with data, previously gathered into the MeteoBuffer
	vecMeteo.clear();
	vecStation.clear();
	for (unsigned int ii=0; ii<unfilteredMeteoBuffer.size(); ii++) {
		vecMeteo.push_back(unfilteredMeteoBuffer[ii].getMeteoData(index));
		vecStation.push_back(unfilteredMeteoBuffer[ii].getStationData(index));
	}
}


void A3DIO::read1DMeteo(const Date_IO& date_in, MeteoData& meteo_out)
{
	StationData sd_tmp;
	read1DMeteo(date_in, meteo_out, sd_tmp);
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

void A3DIO::read1DMeteo(const Date_IO& date_in, MeteoData& meteo_out, StationData& station_out)
{
	double latitude=nodata, longitude=nodata, xcoord=nodata, ycoord=nodata, altitude=nodata;
	string tmp="", line="";
	Date_IO tmp_date;
	vector<string> tmpvec;
	map<string, string> header; // A map to save key value pairs of the file header
	int parsedLines = 0;
	MeteoData tmpdata;

	unfilteredMeteoBuffer.push_back(MeteoBuffer(600));
	filteredMeteoBuffer.push_back(MeteoBuffer(600));

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
		IOUtils::readKeyValueHeader(header, fin, 5, "="); //Read in 5 lines as header
		IOUtils::getValueForKey(header, "Latitude", latitude);
		IOUtils::getValueForKey(header, "Longitude", longitude);
		IOUtils::getValueForKey(header, "X_Coord", xcoord);
		IOUtils::getValueForKey(header, "Y_Coord", ycoord);
		IOUtils::getValueForKey(header, "Altitude", altitude);

		//calculate coordinates if necessary
		if(latitude==IOUtils::nodata || longitude==IOUtils::nodata) {
			if(xcoord==IOUtils::nodata || ycoord==IOUtils::nodata) {
				throw InvalidFormatException("Too many nodata values for coordinates conversion in file " + tmp, AT);
			}
			CH1903_to_WGS84(xcoord, ycoord, latitude, longitude);
		} else if(xcoord!=IOUtils::nodata || ycoord!=IOUtils::nodata) {
			double tmp_lat, tmp_lon;
			CH1903_to_WGS84(xcoord, ycoord, tmp_lat, tmp_lon);
			if(!checkEpsilonEquality(latitude, tmp_lat, 1.e-4) || !checkEpsilonEquality(longitude, tmp_lon, 1.e-4)) {
				throw InvalidArgumentException("Latitude/longitude and Xcoord/Ycoord don't match in header of file " + tmp, AT);
			}
		}
		station_out.setStationData(xcoord, ycoord, altitude, "", latitude, longitude);

		//Read one line, construct Date_IO object and see whether date is greater or equal than the date_in object
		IOUtils::skipLines(fin, 1, eoln); //skip rest of line

		do {
			getline(fin, line, eoln); //read complete line
			parsedLines++;
			//cout << line << endl;

			readMeteoDataLine(line, date_in, tmpdata, tmp);
			tmpdata.cleanData();
			convertUnits(tmpdata);

			unfilteredMeteoBuffer[0].put(tmpdata, station_out);

			//} while(tmp_date<date_in);

		} while(tmpdata.date < date_in);

		meteo_out = tmpdata; //copy by value
		//meteo_out.cleanData();
		//convertUnits(meteo_out);

		//Lese weiter in den Buffer (bis zu 400 Werte)
		try {
			for (int ii=0; ii<400; ii++) {
				getline(fin, line, eoln); //read complete line
				parsedLines++;
				//cout << line << endl;
	
				readMeteoDataLine(line, date_in, tmpdata, tmp);
				tmpdata.cleanData();
				convertUnits(tmpdata);
	
				unfilteredMeteoBuffer[0].put(tmpdata, station_out);
			}
		} catch(...){;} //just continue
    
		//cout << tmp_date.toString() << endl;
	} catch(...) {
		cout << "[e] " << AT << ": "<< endl;
		cleanup();
		throw;
	}

	cleanup();
}

void A3DIO::readMeteoDataLine(std::string& line, const Date_IO& date_in, MeteoData& tmpdata, string filename)
{
	Date_IO tmp_date;
	int tmp_ymdh[4];
	vector<string> tmpvec;
	double tmp_values[6];

	if (IOUtils::readLineToVec(line, tmpvec) != 10) {
		throw InvalidFormatException("Premature End of Line or no data for date " + date_in.toString() + " found in File " + filename, AT);
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
}


void A3DIO::read2DMeteo(const Date_IO& date_in, vector<MeteoData>& meteo_out)
{
	vector<StationData> vecStation;
	read2DMeteo(date_in, meteo_out, vecStation);
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
void A3DIO::read2DMeteo(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{
	unsigned int stations=0, bufferindex=0;
	map<string, unsigned int> hashStations = map<string, unsigned int>();
	vector<string> filenames = vector<string>();

	vecMeteo.clear();
	vecStation.clear();
 
	constructMeteo2DFilenames(unfilteredMeteoBuffer[0].getMeteoData(0).date, filenames);
	stations = getNrOfStations(filenames, hashStations);
	cerr << "[I] Number of 2D meteo stations: " << stations << endl;
  
	if (stations < 1) {
		throw InvalidFormatException("No StationData found in 2D Meteo Files", AT); 
	} else {
		vecMeteo.clear();
		vecStation.clear();

		MeteoBuffer tmpbuffer(600, unfilteredMeteoBuffer[0].size());

		for (unsigned int ii=0; ii<stations; ii++) {
			StationData tmp_sd;
			MeteoData tmp_md;
      
			unfilteredMeteoBuffer.push_back(tmpbuffer);

			vecMeteo.push_back(tmp_md);
			vecStation.push_back(tmp_sd);
		}
	}

	try {
		read2DMeteoHeader(filenames[0], hashStations, vecStation);
		read2DMeteoHeader(filenames[1], hashStations, vecStation);
		read2DMeteoHeader(filenames[2], hashStations, vecStation);
		read2DMeteoHeader(filenames[3], hashStations, vecStation);
		if(IOUtils::fileExists(filenames[4])) { //for keeping dw optional
			read2DMeteoHeader(filenames[4], hashStations, vecStation);
		}

		for (unsigned int ii=1; ii<unfilteredMeteoBuffer.size(); ii++) {//loop through all MeteoBuffers
			for (unsigned int jj=0; jj<unfilteredMeteoBuffer[ii].size(); jj++) { //loop through all allocated StationData within one MeteoBuffer
				StationData& bufferedStation = unfilteredMeteoBuffer[ii].getStationData(jj);
				bufferedStation = vecStation[ii-1];
			}
		}

		do {

			unsigned int currentindex = bufferindex;
			read2DMeteoData(filenames[0], "nswc", date_in, hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[1], "rh", date_in, hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[2], "ta", date_in, hashStations, vecMeteo, bufferindex);
			bufferindex = currentindex;
			read2DMeteoData(filenames[3], "vw", date_in, hashStations, vecMeteo, bufferindex);
			
			if(IOUtils::fileExists(filenames[4])) { //for keeping dw optional
				bufferindex = currentindex;
				read2DMeteoData(filenames[4], "dw", date_in, hashStations, vecMeteo, bufferindex);
			}
			//cerr << "bufferindex: " << bufferindex << "  Expected size()" << unfilteredMeteoBuffer[0].size() << endl;

			if (bufferindex < (unfilteredMeteoBuffer[0].size()-1)) {
				//construct new filenames for the continued buffering
				constructMeteo2DFilenames(unfilteredMeteoBuffer[0].getMeteoData(bufferindex).date, filenames);
			}
		} while(bufferindex < (unfilteredMeteoBuffer[0].size()));
	} catch(...) {
		vecMeteo.clear();
		vecStation.clear();
    
		//TODO: cleanup of buffer

		cleanup();
		throw;
	}
  
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++) {
		vecMeteo[ii].cleanData();
		convertUnits(vecMeteo[ii]);
	}

	for (unsigned int ii=1; ii<unfilteredMeteoBuffer.size(); ii++) {//loop through all MeteoBuffers
		for (unsigned int jj=0; jj<unfilteredMeteoBuffer[ii].size(); jj++){ //loop through all allocated MeteoData within one MeteoBuffer
			MeteoData& bufferedMeteo = unfilteredMeteoBuffer[ii].getMeteoData(jj);
			bufferedMeteo.cleanData();
			convertUnits(bufferedMeteo);
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

void A3DIO::read2DMeteoData(const string& filename, const string& parameter, const Date_IO& date_in, 
						    map<string,unsigned int>& hashStations, vector<MeteoData>& vecM, unsigned int& bufferindex)
{
	string line_in = "";
	unsigned int columns;
	vector<string> tmpvec, vec_names;
	Date_IO tmp_date;
	int tmp_ymdh[4];
	//unsigned int bufferindex=0;
	bool old = true;

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

	MeteoData& lastMeteoData = unfilteredMeteoBuffer[0].getMeteoData(unfilteredMeteoBuffer[0].size()-1); //last time stamp in buffer of 1D meteo
	do {
		getline(fin, line_in, eoln); 
		string tmpline = line_in;
		IOUtils::trim(tmpline);

		if (tmpline=="") {
			break;
		}

		if (IOUtils::readLineToVec(line_in, tmpvec)!=columns) {
			throw InvalidFormatException("Premature End of Line or no data for date " + date_in.toString() 
								    + " found in File " + filename, AT);
		}
    
		for (int ii=0; ii<4; ii++) {
			if (!IOUtils::convertString(tmp_ymdh[ii], tmpvec[ii], std::dec)) {
				throw InvalidFormatException("Check date columns in " + filename, AT);
			}
		}
		tmp_date.setDate(tmp_ymdh[0],tmp_ymdh[1],tmp_ymdh[2],tmp_ymdh[3]);

		MeteoData& currentMeteoData = unfilteredMeteoBuffer[0].getMeteoData(bufferindex); //1D Element to synchronize date
		if (tmp_date == currentMeteoData.date) {
			//Read in data
			for (unsigned int ii=4; ii<columns; ii++) {
				unsigned int stationnr = hashStations[vec_names.at(ii)]; 
				MeteoData& tmpmd = unfilteredMeteoBuffer[stationnr].getMeteoData(bufferindex);
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
		//} while(tmp_date<date_in);
		if ((tmp_date >= date_in) && (old==true)) {
			old=false;
			//Write Data into MeteoData objects; which data is in the file determined by switch 'parameter'
			for (unsigned int ii=4; ii<columns; ii++) {
				unsigned int stationnr = hashStations[vec_names.at(ii)]; 
				vecM.at(stationnr-1).date = tmp_date;

				if (parameter == "nswc") {
					if (!IOUtils::convertString(vecM.at(stationnr-1).nswc, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For nswc value in " + filename + "  for date " + tmp_date.toString(), AT);
					}
    
				} else if (parameter == "rh") {
					if (!IOUtils::convertString(vecM.at(stationnr-1).rh, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For rh value in " + filename + "  for date " + tmp_date.toString(), AT);
					}
	  
				} else if (parameter == "ta") {
					if (!IOUtils::convertString(vecM.at(stationnr-1).ta, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For ta value in " + filename + "  for date " + tmp_date.toString(), AT);
					}
	  
				} else if (parameter == "vw") {
					if (!IOUtils::convertString(vecM.at(stationnr-1).vw, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For vw value in " + filename + "  for date " + tmp_date.toString(), AT);
					}
				} else if (parameter == "dw") {
					if (!IOUtils::convertString(vecM.at(stationnr-1).dw, tmpvec[ii], std::dec)) {
						throw ConversionFailedException("For dw value in " + filename + "  for date " + tmp_date.toString(), AT);
					}
				}
			}
		}

	} while((tmp_date<lastMeteoData.date) && (!fin.eof()));

	cleanup();
}

void A3DIO::read2DMeteoHeader(const string& filename, map<string,unsigned int>& hashStations, vector<StationData>& vecS)
{
	string line_in = "";
	unsigned int columns = 0;
	vector<string> vec_altitude, vec_xcoord, vec_ycoord, vec_names;

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
			CH1903_to_WGS84(vecS[stationnr-1].eastCoordinate, vecS[stationnr-1].northCoordinate, vecS[stationnr-1].latitude, vecS[stationnr-1].longitude);
		}
	}

	cleanup();
}

void A3DIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	int yyyy, mm, dd, hh;
	date_in.getDate(yyyy, mm, dd, hh);
	string filepath="";

	cfg.getValue("DAPATH", filepath); // cout << tmp << endl;
  
	stringstream ss;
	ss.fill('0');
	ss << filepath << "/" << setw(4) << yyyy << setw(2) << mm << setw(2) <<  dd << setw(2) <<  hh << ".sca";

	read2DGrid(da_out, ss.str());
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
	pts.SetSize(mypts.size());

	for (unsigned int jj=0; jj<mypts.size(); jj++) {
		pts[jj].ix = mypts.at(jj).first;
		pts[jj].iy = mypts.at(jj).second;
	}
}

void A3DIO::write2DGrid(const Grid2DObject& grid_in, const string& filename) {
  
	if (IOUtils::fileExists(filename)) {
		throw IOException("File " + filename + " already exists!", AT);
	}

	fout.open(filename.c_str());
	if (fout.fail()) {
		throw FileAccessException(filename, AT);
	}


	try {
		fout << "ncols \t\t" << grid_in.ncols << endl;
		fout << "nrows \t\t" << grid_in.nrows << endl;
		fout << "xllcorner \t" << grid_in.xllcorner << endl;
		fout << "yllcorner \t" << grid_in.yllcorner << endl;    
		fout << "cellsize \t" << grid_in.cellsize << endl;
		fout << "NODATA_value \t" << (int)(IOUtils::nodata) << endl;

		for (unsigned int kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
			for (unsigned int ll=0; ll < grid_in.ncols; ll++){
				fout << grid_in.grid2D(ll, kk) << "\t";
			}
			fout << endl;
		}
	} catch(...) {
		cout << "[E] " << AT << ": "<< endl;
		cleanup();
		throw;
	}

	cleanup();
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
