#ifdef _POPC_
#include "IOHandler.ph"
#else
#include "IOHandler.h"
#endif

using namespace std;

void IOHandler::registerPlugins()
{
	mapPlugins["FILE"]      = IOPlugin("", "A3DIO", &fileio, NULL);
	mapPlugins["BOSCHUNG"]  = IOPlugin("libboschungio.so", "BoschungIO", NULL, NULL);
	mapPlugins["IMIS"]  = IOPlugin("libimisio.so", "ImisIO", NULL, NULL);
	mapPlugins["GEOTOP"]  = IOPlugin("libgeotopio.so", "GeotopIO", NULL, NULL);
	mapPlugins["GSN"]  = IOPlugin("libgsnio.so", "GSNIO", NULL, NULL);
	mapPlugins["ESRI"]  = IOPlugin("libesriio.so", "EsriIO", NULL, NULL);
}

#ifdef _POPC_
IOHandler::IOHandler(const std::string& configfile) :  cfg(configfile), fileio(cfg){
#else
IOHandler::IOHandler(const std::string& configfile) : IOInterface(NULL), cfg(configfile), fileio(cfg)
{
#endif
	registerPlugins();
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
IOHandler::IOHandler(const ConfigReader& cfgreader) : cfg(cfgreader), fileio(cfg)
#else
IOHandler::IOHandler(const ConfigReader& cfgreader) : IOInterface(NULL), cfg(cfgreader), fileio(cfg)
#endif
{
	registerPlugins();	
}

#ifdef _POPC_
IOHandler::~IOHandler(){
#else
IOHandler::~IOHandler() throw(){
#endif
	// Get rid of the objects
	for (mapit = mapPlugins.begin(); mapit!=mapPlugins.end(); mapit++){
		deletePlugin((mapit->second).dynLibrary, (mapit->second).io);
	}
}
#ifdef _POPC_
void IOHandler::deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io)
#else
void IOHandler::deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io) throw()
#endif
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
	cout << "[i] " << AT << ": Loading dynamic plugin: " << libname << endl;
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

IOInterface* IOHandler::getPlugin(const std::string& cfgvalue)
{
	string op_src="";
	cfg.getValue(cfgvalue, op_src); 

	mapit = mapPlugins.find(op_src);		
	if (mapit == mapPlugins.end())
		throw IOException(cfgvalue + "does not seem to be valid descriptor in file " + cfg.getFileName(), AT);
	
	if ((mapit->second).io == NULL){
		loadPlugin((mapit->second).libname, (mapit->second).classname, (mapit->second).dynLibrary, (mapit->second).io);
	}			
	
	if ((mapit->second).io == NULL)	
		throw IOException("Requesting to read/write data with plugin for " + cfgvalue + ", but plugin is not loaded", AT);			

	return (mapit->second).io;
}

void IOHandler::read2DGrid(Grid2DObject& _grid, const string& _filename)
{
	IOInterface *plugin = getPlugin("GRID2DSRC");
	plugin->read2DGrid(_grid, _filename);
}

void IOHandler::readDEM(Grid2DObject& dem_out)
{
	IOInterface *plugin = getPlugin("DEMSRC");
	plugin->readDEM(dem_out);
}

void IOHandler::readLanduse(Grid2DObject& landuse_out)
{
	IOInterface *plugin = getPlugin("LANDUSESRC");
	plugin->readLanduse(landuse_out);
}

void IOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
						std::vector< std::vector<MeteoData> >& vecMeteo, 
						std::vector< std::vector<StationData> >& vecStation,
						const unsigned int& stationindex)
{
	IOInterface *plugin = getPlugin("METEOSRC");
	plugin->readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
}

void IOHandler::readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out)
{
	IOInterface *plugin = getPlugin("DASRC");
	plugin->readAssimilationData(date_in, da_out);
}

void IOHandler::readSpecialPoints(CSpecialPTSArray& pts) {
	IOInterface *plugin = getPlugin("SPECIALPTSSRC");
	plugin->readSpecialPoints(pts);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const string& name)
{
	IOInterface *plugin = getPlugin("OUTPUT");
	plugin->write2DGrid(grid_in, name);
}


