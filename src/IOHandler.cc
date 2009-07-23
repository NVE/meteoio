#ifdef _POPC_
#include "IOHandler.ph"
#else
#include "IOHandler.h"
#endif

using namespace std;


#ifdef _POPC_
IOHandler::IOHandler(const std::string& configfile) :  cfg(configfile), fileio(cfg){
#else
IOHandler::IOHandler(const std::string& configfile) : IOInterface(NULL), cfg(configfile), fileio(cfg)
{
#endif
	ascii_src = "FILE";
	boschung_src = "BOSCHUNG";
	imis_src = "IMIS";
	//load all dynamic plugins
	loadDynamicPlugins();
}

//Copy constructor
#ifdef _POPC_
//IOHandler::IOHandler(const IOHandler& aio) : cfg(aio.cfg), fileio(cfg), boschungio(cfg), imisio(cfg){
	//Nothing else so far
//}
#else
IOHandler::IOHandler(const IOHandler& aio) : IOInterface(NULL), cfg(aio.cfg), fileio(cfg)
{
	//Nothing else so far
	//TODO: Deal with the IOInterface* pointers, e.g. boschungio
}
#endif

#ifdef _POPC_
/*IOHandler::IOHandler(const ConfigReader& cfgreader) : cfg(cfgreader), fileio(cfg), boschungio(cfg), imisio(cfg){
  IOHandler::ascii_src = "FILE";
  IOHandler::boschung_src = "BOSCHUNG";
  IOHandler::imis_src = "IMIS";
  loadDynamicPlugins();
}*/
#else
IOHandler::IOHandler(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader), fileio(cfg)
{
	//Nothing else so far
	loadDynamicPlugins();
}
#endif

#ifdef _POPC_
IOHandler::~IOHandler(){
#else
IOHandler::~IOHandler() throw(){
#endif
	// Get rid of the object

	if (dynLibraryBoschung != NULL) {
		if (boschungio != NULL) {
			boschungio->deleteSelf();
			boschungio = NULL;
		}
		// Close the dynamic library
		delete dynLibraryBoschung;
	}

	if (dynLibraryImis != NULL) {
		if (imisio != NULL) {
			imisio->deleteSelf();
			imisio = NULL;
		}
		// Close the dynamic library
		delete dynLibraryImis;
	}

	cleanup();
}

//Clone function
//IOHandler* IOHandler::clone() const { return new IOHandler(*this); }
#ifdef _POPC_
void IOHandler::cleanup(){
#else
void IOHandler::cleanup() throw(){
#endif
}

void IOHandler::loadDynamicPlugins()
{
	string pluginpath = "";

	cout << "[i] " << AT << ": Loading dynamic plugins:" << endl;
	try {
		cfg.getValue("PLUGINPATH", pluginpath); 

		//BoschungIO dynamic library needs to be loaded
		cout << "\t" << "Trying to load libboschungio.so ... ";
		string filename = pluginpath + "/libboschungio.so";
		dynLibraryBoschung = DynamicLoader::loadObjectFile(filename, RTLD_NOW);

		if(dynLibraryBoschung == NULL) {
			cout << "failed\n\tCouldn't load the dynamic library " << filename << "\n\t" << dlerror() << endl;
		} else {
			boschungio = dynamic_cast<IOInterface*>(dynLibraryBoschung->newObject("BoschungIO", cfg.getFileName()));

			if(boschungio == NULL) {
				cout << "failed" << endl;
				//delete dynLibraryBoschung; This causes a segfault !!
			} else {
				cout << "success" << endl;
			}
		}
	} catch (exception& e) {
		if (dynLibraryBoschung != NULL)
		delete dynLibraryBoschung;
		cerr << "\t" << e.what() << endl;
	}

	try {
		cfg.getValue("PLUGINPATH", pluginpath); 

		//ImisIO dynamic library needs to be loaded
		cout << "\t" << "Trying to load libimisio.so ... ";
		string filename = pluginpath + "/libimisio.so";
		dynLibraryImis = DynamicLoader::loadObjectFile(filename, RTLD_NOW);

		if(dynLibraryImis == NULL) {
			cout << "failed\n\tCouldn't load the dynamic library " << filename << "\n\t" << dlerror() << endl;
		} else {
			imisio = dynamic_cast<IOInterface*>(dynLibraryImis->newObject("ImisIO", cfg.getFileName()));

			if(imisio == NULL) {
				cout << "failed" << endl;
				//delete dynLibraryImis; This causes a segfault !!
			} else {
				cout << "success" << endl;
			}
		}
	} catch (exception& e) {
		if (dynLibraryImis != NULL) {
			delete dynLibraryImis;
		}
		cerr << "\t" << e.what() << endl;
	}
}
 
void IOHandler::get2DGridSize(int& nx, int& ny)
{
	string demsource="";

	try {
		cfg.getValue("DEMSRC", demsource);

		if (demsource==ascii_src){
			//A3DIO fileio(cfg.getFileName());
			fileio.get2DGridSize(nx, ny);

		} else if (demsource==boschung_src){
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (demsource==imis_src){
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("DEMSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...){
		throw;
	}
}

void IOHandler::read2DGrid(Grid2DObject&, const string& filename)
{
	//Nothing so far
	(void)filename;
	throw IOException("Nothing implemented here", AT);
}

void IOHandler::readDEM(Grid2DObject& dem_out)
{
	string demsource="";

	try {
		cfg.getValue("DEMSRC", demsource); // cout << tmp << endl;

		if (demsource==ascii_src) {
			//A3DIO fileio(cfg.getFileName());
			fileio.readDEM(dem_out);

		} else if (demsource==boschung_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (demsource==imis_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("DEMSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}

void IOHandler::readLanduse(Grid2DObject& landuse_out)
{
	string landusesource="";

	try {
		cfg.getValue("LANDUSESRC", landusesource); 

		if (landusesource==ascii_src) {
			//A3DIO fileio(cfg.getFileName());
			fileio.readLanduse(landuse_out);

		} else if (landusesource==boschung_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (landusesource==imis_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("LANDUSESRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}


void IOHandler::readMeteoData(const Date_IO& date_in, vector<MeteoData>& meteo_out)
{
	vector<StationData> vecStation;
	readMeteoData(date_in, meteo_out, vecStation);
}

void IOHandler::readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
{
	//See whether data is already buffered
	//Filter
	//Resample if necessary
	//Filter resampled value
	//return that value
	string meteo2dsource="";

	try {
		cfg.getValue("METEOSRC", meteo2dsource); 

		if (meteo2dsource==ascii_src) {
			//A3DIO fileio(cfg.getFileName());
			fileio.readMeteoData(date_in, vecMeteo, vecStation);

		} else if (meteo2dsource==boschung_src) {
			if (boschungio != NULL) {
				boschungio->readMeteoData(date_in, vecMeteo, vecStation);
			} else {
				throw IOException("Requesting to read data with plugin libBoschungIO.so, but plugin is not loaded", AT);
			}
		} else if (meteo2dsource==imis_src) {
			if (imisio != NULL) {
				imisio->readMeteoData(date_in, vecMeteo, vecStation);
			} else {
				throw IOException("Requesting to read data with plugin libImisIO.so, but plugin is not loaded", AT);
			}
		} else {
			throw IOException("METEOSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}

}

void IOHandler::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	string dasrc="";

	try {
		cfg.getValue("DASRC", dasrc); 

		if (dasrc==ascii_src) {
			//A3DIO fileio(cfg.getFileName());
			fileio.readAssimilationData(date_in, da_out);

		} else if (dasrc==boschung_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (dasrc==imis_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("DASRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}

void IOHandler::readSpecialPoints(CSpecialPTSArray& pts) {
	string specialptssrc="";

	try {
		cfg.getValue("SPECIALPTSSRC", specialptssrc); 

		if (specialptssrc==ascii_src) {
			//A3DIO fileio(cfg.getFileName());
			fileio.readSpecialPoints(pts);

		} else if (specialptssrc==boschung_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (specialptssrc==imis_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("SPECIALPTSSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const string& name)
{
	string output="";

	try {
		cfg.getValue("OUTPUT", output); 

		if (output==ascii_src) {
			fileio.write2DGrid(grid_in, name);
		} else if (output==boschung_src) {
			throw IOException("Nothing implemented here", AT);

		} else if (output==imis_src) {
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("OUTPUT does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}

void IOHandler::write2DGrid(const Array2D<double>& grid_in, const double& xllcorner, const double& yllcorner, const double& cellsize, const string& name)
{
	string output="";

	try {
		cfg.getValue("OUTPUT", output); 

		if (output==ascii_src) {
			fileio.write2DGrid(grid_in, xllcorner, yllcorner, cellsize, name);
		} else if (output==boschung_src) {
			throw IOException("Nothing implemented here", AT);

		} else if (output==imis_src) {
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("OUTPUT does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}

