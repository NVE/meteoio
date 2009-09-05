#include "GrassIO.h"

using namespace std;

GrassIO::GrassIO(void (*delObj)(void*), const string& filename) : IOInterface(delObj), cfg(filename){}

GrassIO::GrassIO(const std::string& configfile) : IOInterface(NULL), cfg(configfile)
{
	//Nothing else so far
}

GrassIO::GrassIO(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader)
{
	//Nothing else so far
}

GrassIO::~GrassIO() throw()
{
	cleanup();
}

void GrassIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void GrassIO::read2DGrid(Grid2DObject& grid_out, const string& filename)
{

	int _nx, _ny;
	unsigned int ncols, nrows;
	double north, east, south, west, latitude, longitude;
	double tmp_val, xllcorner, yllcorner, cellsize;
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
		IOUtils::readKeyValueHeader(header, fin, 6, ":");
		IOUtils::getValueForKey(header, "cols",  _nx);
		IOUtils::getValueForKey(header, "rows",  _ny);
		IOUtils::getValueForKey(header, "north", north);
		IOUtils::getValueForKey(header, "east",  east);
		IOUtils::getValueForKey(header, "south", south);
		IOUtils::getValueForKey(header, "west",  west);

		if ((_nx==0) || (_ny==0)) {
			throw IOException("Number of rows or columns in 2D Grid given is zero, in file: " + filename, AT);
		}
		if((_nx<0) || (_ny<0)) {
			throw IOException("Number of rows or columns in 2D Grid read as \"nodata\", in file: " + filename, AT);
		}
		ncols = (unsigned int)_nx;
		nrows = (unsigned int)_ny;
		xllcorner = west;
		yllcorner = south;
		cellsize = (east - west) / (double)ncols;

		string coordsys="", coordparam="";
		try {
			cfg.getValue("COORDIN", coordsys);
			cfg.getValue("COORDPARAM", coordparam); 
		} catch(std::exception& e){
			//problems while reading values for COORDIN or COORDPARAM
		}

		//compute WGS coordinates (considered as the true reference)

		//HACK: check how we can input coordinates as WGS84 directly.
		//this HACK is a very cheap "extension" of the dem file format...
		if(xllcorner<360. || yllcorner<360.) {
			latitude = xllcorner;
			longitude = yllcorner;
		} else {
			MapProj mymapproj(coordsys, coordparam);
			mymapproj.convert_to_WGS84(xllcorner, yllcorner, latitude, longitude);
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
				if (tmpvec[ll] == "*"){
					tmp_val = IOUtils::nodata;
				} else {
					if (!IOUtils::convertString(tmp_val, tmpvec[ll], std::dec)) {
						throw ConversionFailedException("For Grid2D value in line: " + line + " in file " + filename, AT);
					}
				}
				
				if(tmp_val <= IOUtils::nodata) {
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

void GrassIO::readDEM(DEMObject& dem_out)
{
	string filename="";
	cfg.getValue("DEMFILE", filename); // cout << tmp << endl;
	read2DGrid(dem_out, filename);
}

void GrassIO::readLanduse(Grid2DObject& landuse_out)
{
	string filename="";
	cfg.getValue("LANDUSEFILE", filename); // cout << tmp << endl;
	read2DGrid(landuse_out, filename);
}

void GrassIO::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
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

void GrassIO::readMeteoData(const Date_IO&, const Date_IO&, 
					 std::vector< std::vector<MeteoData> >&, 
					 std::vector< std::vector<StationData> >&,
					 const unsigned int&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GrassIO::readSpecialPoints(CSpecialPTSArray&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GrassIO::write2DGrid(const Grid2DObject& grid_in, const string& name)
{  
	fout.open(name.c_str());
	if (fout.fail()) {
		throw FileAccessException(name, AT);
	}

	fout << setprecision(6) << fixed;

	try {
		fout << "north:" << (grid_in.yllcorner+grid_in.cellsize*grid_in.nrows) << endl;    
		fout << "south:" << grid_in.yllcorner << endl;    
		fout << "east:"  << (grid_in.xllcorner+grid_in.cellsize*grid_in.ncols)  << endl;
		fout << "west:"  << grid_in.xllcorner << endl;
		fout << "rows:"  << grid_in.nrows << endl;
		fout << "cols:"  << grid_in.ncols << endl;

		for (unsigned int kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
			unsigned int ll = 0;
			for (ll=0; ll < (grid_in.ncols-1); ll++){
				if (grid_in.grid2D(ll,kk) == IOUtils::nodata) {
					fout << "* ";
				} else {
					fout << grid_in.grid2D(ll, kk) << " ";
				}
			}

			//The last value in a line does not have a trailing " "
			if (grid_in.grid2D(ll,kk) == IOUtils::nodata) {
				fout << "*";
			} else {
				fout << grid_in.grid2D(ll, kk);
			}
			fout << endl;
		}
	} catch(std::exception& e) {
		cout << "[E] " << AT << ": " << e.what() << endl;
		cleanup();
		throw;
	}

	cleanup();
}

extern "C"
{
	void deleteObject(void* obj) {
		delete reinterpret_cast<PluginObject*>(obj);
	}
  
	void* loadObject(const string& classname, const string& filename) {
		if(classname == "GrassIO") {
			//cerr << "Creating dynamic handle for " << classname << endl;
			return new GrassIO(deleteObject, filename);
		}
		//cerr << "Could not load " << classname << endl;
		return NULL;
	}
}
