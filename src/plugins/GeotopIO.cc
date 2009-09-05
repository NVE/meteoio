#include "GeotopIO.h"

using namespace std;

GeotopIO::GeotopIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

GeotopIO::GeotopIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	//Nothing else so far
}

GeotopIO::GeotopIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	//Nothing else so far
}

GeotopIO::~GeotopIO() throw()
{
	cleanup();
}

void GeotopIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
}

void GeotopIO::read2DGrid(Grid2DObject&, const string& filename)
{
	//Nothing so far
	(void)filename;
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readDEM(DEMObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readLanduse(Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex)
{

	vector<string> tmpvec;
	string line="", filename="", str_stations="", path="", prefix="";
	unsigned int stations=0;

	(void)stationindex;
	vecMeteo.clear();
	vecStation.clear();

	cfg.getValue("NROFSTATIONS", str_stations);

	if (!IOUtils::convertString(stations, str_stations, std::dec)) {
		throw ConversionFailedException("Error while reading value for NROFSTATIONS", AT);
	}

	cfg.getValue("METEOPATH", path); // cout << tmp << endl;
	cfg.getValue("METEOPREFIX", prefix); // cout << tmp << endl;

	for (unsigned int ii=0; ii<stations; ii++) {
		vecMeteo.push_back( vector<MeteoData>() );
		vecStation.push_back( vector<StationData>() );

		stringstream ss;
		ss.fill('0');
		ss << path << "/" << prefix << setw(4) << (ii+1) << ".txt";

		filename = ss.str();
		//cout << ss.str() << endl;

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
			getline(fin, line, eoln); //read complete line meta information
			//cout << line << endl;

			unsigned int ncols = IOUtils::readLineToVec(line, tmpvec, ',');

			if (ncols == 0)
				throw InvalidFormatException("No meta data found in " + filename, AT);

			vector<double> tmpdata = vector<double>(ncols);
			vector<int> ymdh = vector<int>(4);
			while (!fin.eof()){
				getline(fin, line, eoln); //read complete line of data
				//cout << line << endl;

				MeteoData md;
				if (IOUtils::readLineToVec(line, tmpvec, ',') != ncols) {
					break;
					//throw InvalidFormatException("Premature End " + filename, AT);
				}

				if (!IOUtils::convertString(ymdh[0], "20" + tmpvec[0].substr(6,2), std::dec)) //day
					throw InvalidFormatException(filename + ": " + line, AT);
				if (!IOUtils::convertString(ymdh[1], tmpvec[0].substr(3,2), std::dec)) //month
					throw InvalidFormatException(filename + ": " + line, AT);
				if (!IOUtils::convertString(ymdh[2], tmpvec[0].substr(0,2), std::dec)) //year
					throw InvalidFormatException(filename + ": " + line, AT);
				if (!IOUtils::convertString(ymdh[3], tmpvec[0].substr(9,2), std::dec)) //hour
					throw InvalidFormatException(filename + ": " + line, AT);

				for (unsigned int jj=1; jj<ncols; jj++) {
					if (!IOUtils::convertString(tmpdata[jj], tmpvec.at(jj), std::dec))
						throw InvalidFormatException(filename + ": " + line, AT);
				}

				md.setMeteoData(Date_IO(ymdh[0],ymdh[1],ymdh[2],ymdh[3]), 
								 tmpdata[4], 
								 IOUtils::nodata, 
								 tmpdata[1], 
								 IOUtils::nodata, 
								 tmpdata[2], 
								 IOUtils::nodata, 
								 tmpdata[5], 
								 IOUtils::nodata, 
								 IOUtils::nodata, 
								 IOUtils::nodata, 
								 IOUtils::nodata);
				md.p = tmpdata[3]; //hPa
				//cout << md.toString() << endl;

				vecMeteo[ii].push_back(md);
				vecStation[ii].push_back(StationData());
			}
		} catch(std::exception& e) {
			cleanup();
			throw;
		}
		fin.close();
	}
}

void GeotopIO::readAssimilationData(const Date_IO&, Grid2DObject&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::readSpecialPoints(CSpecialPTSArray&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::write2DGrid(const Grid2DObject&, const string& name)
{
	//Nothing so far
	(void)name;
	throw IOException("Nothing implemented here", AT);
}

void GeotopIO::convertUnits(MeteoData& meteo)
{
	//converts C to Kelvin, converts lwr to ea, converts RH to [0,1]
	if(meteo.ta==nodata) {
		meteo.ta=nodata;
	} else {
		meteo.ta=C_TO_K(meteo.ta);
	}
	
	if(meteo.tsg==nodata) {
		meteo.tsg=nodata;
	} else {
		meteo.tsg=C_TO_K(meteo.tss);
	}
	
	if(meteo.tss==nodata) {
		meteo.tss=nodata;
	} else {
		meteo.tss=C_TO_K(meteo.tss);
	}

	if(meteo.rh==nodata) {
		meteo.rh=nodata;
	} else {
		meteo.rh /= 100.;
	}
}

extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
  
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "GeotopIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GeotopIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
