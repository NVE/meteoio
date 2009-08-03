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
	geotop_src = "GEOTOP";

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
	ascii_src = "FILE";
	boschung_src = "BOSCHUNG";
	imis_src = "IMIS";
	geotop_src = "GEOTOP";

	//Nothing else so far
	loadDynamicPlugins();
}
#endif

#ifdef _POPC_
IOHandler::~IOHandler(){
#else
IOHandler::~IOHandler() throw(){
#endif
	// Get rid of the objects
	deletePlugin(dynLibraryImis, imisio);
	deletePlugin(dynLibraryBoschung, boschungio);
	deletePlugin(dynLibraryGeoTOP, geotopio);

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
	cout << "[i] " << AT << ": Loading dynamic plugins:" << endl;

	loadPlugin("libboschungio.so", "BoschungIO", dynLibraryBoschung, boschungio);
	loadPlugin("libimisio.so", "ImisIO", dynLibraryImis, imisio);
	loadPlugin("libgeotopio.so", "GeotopIO", dynLibraryGeoTOP, geotopio);
}

void IOHandler::deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io) throw()
{
	if (dynLibrary != NULL) {
		if (io != NULL) {
			io->deleteSelf();
			io = NULL;
		}

		// Close the dynamic library
		delete dynLibrary;
	}
}

void IOHandler::loadPlugin(const string& libname, const string& classname, DynamicLibrary*& dynLibrary, IOInterface*& io)
{
	string pluginpath = "";

	try {
		cfg.getValue("PLUGINPATH", pluginpath); 

		//Which dynamic library needs to be loaded
		cout << "\t" << "Trying to load " << libname << " ... ";
		string filename = pluginpath + "/" + libname;
		dynLibrary = DynamicLoader::loadObjectFile(filename, RTLD_NOW);
		
		if(dynLibrary == NULL) {
			cout << "failed\n\tCouldn't load the dynamic library " << filename << "\n\t" << dlerror() << endl;
		} else {
			io = dynamic_cast<IOInterface*>((dynLibrary)->newObject(classname, cfg.getFileName()));

			if(io == NULL) {
				cout << "failed" << endl;
				//delete dynLibrary; This causes a segfault !!
			} else {
				cout << "success" << endl;
			}
			}
	} catch (exception& e) {
		if (dynLibrary != NULL)
		delete dynLibrary;
		cerr << "\t" << e.what() << endl;
	}
}

 
void IOHandler::read2DGrid(Grid2DObject& _grid, const string& _filename)
{
	string grid2Dsrc=""; //TODO: io.ini needs a field for this GRID2DSRC?

	try {
		cfg.getValue("GRID2DSRC", grid2Dsrc); // cout << tmp << endl;

		if (grid2Dsrc==ascii_src){
			//A3DIO fileio(cfg.getFileName());
			fileio.read2DGrid(_grid, _filename);

		} else if (grid2Dsrc==boschung_src){
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (grid2Dsrc==imis_src){
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("GRID2DSRC does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}
	} catch (...){
		throw;
	}
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


void IOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
						std::vector< std::vector<MeteoData> >& vecMeteo, 
						std::vector< std::vector<StationData> >& vecStation,
						unsigned int stationindex)
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
			fileio.readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);

		} else if (meteo2dsource==boschung_src) {
			if (boschungio != NULL) {
				boschungio->readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
			} else {
				throw IOException("Requesting to read data with plugin libBoschungIO.so, but plugin is not loaded", AT);
			}
		} else if (meteo2dsource==geotop_src) {
			if (geotopio != NULL) {
				geotopio->readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
			} else {
				throw IOException("Requesting to read data with plugin libgeotopio.so, but plugin is not loaded", AT);
			}
		} else if (meteo2dsource==imis_src) {
			if (imisio != NULL) {
				imisio->readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
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
			//A3DIO fileio(cfg.getFileName());
			fileio.write2DGrid(grid_in, name);
		} else if (output==boschung_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else if (output==imis_src) {
			//Nothing so far
			throw IOException("Nothing implemented here", AT);

		} else {
			throw IOException("OUTPUT does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
		}

	} catch (...) {
		throw;
	}
}


